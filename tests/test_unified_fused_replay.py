"""Replica fused replay should select work once and preserve completed products."""
import json

import numpy as np
import pandas as pd
import pytest
from validation_fixture import make_project

from mdscope.pipeline import fused_all, runner
from mdscope.registry import fused_output


def _completed(tmp_path):
    rep, cfg = make_project(tmp_path, 2)
    cfg["progress"]["enabled"] = False
    root = tmp_path / "results"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    return rep, cfg, root, state, root / rep.system / rep.label


def test_plan_retains_named_instances_and_only_pending_products(tmp_path):
    rep, cfg = make_project(tmp_path, 2)
    repout = tmp_path / "results" / rep.system / rep.label
    expected, pending = runner._plan_fused_replay(cfg, repout)
    assert pending == expected
    assert {
        entry.run_id for entry in pending
        if entry.name == "residue_membrane_distance"
    } == {"A", "B"}
    existing = fused_output(repout, "residue_membrane_distance", run_id="A")
    existing.path.parent.mkdir(parents=True)
    existing.path.write_bytes(b"persisted-primary")
    _, selected = runner._plan_fused_replay(cfg, repout)
    assert existing not in selected
    assert fused_output(repout, "residue_membrane_distance", run_id="B") in selected


def test_completed_stream_replays_three_missing_families_in_one_pass(tmp_path, monkeypatch):
    rep, cfg, root, state, repout = _completed(tmp_path)
    targets = (
        fused_output(repout, "bilayer_normal"),
        fused_output(repout, "solvation"),
        fused_output(repout, "residue_membrane_distance", run_id="B"),
    )
    saved = [entry.path.read_bytes() for entry in targets]
    completed_sibling = fused_output(
        repout, "residue_membrane_distance", run_id="A"
    ).path
    sibling_bytes = completed_sibling.read_bytes()
    for entry in targets:
        entry.path.unlink()

    opened = []
    processed = []
    original = fused_all.mda.Universe
    original_track = fused_all.track

    def counting_universe(*args, **kwargs):
        opened.append(args)
        return original(*args, **kwargs)

    def counting_frames(frames, *args, **kwargs):
        for ts in original_track(frames, *args, **kwargs):
            processed.append(int(ts.frame))
            yield ts

    monkeypatch.setattr(fused_all.mda, "Universe", counting_universe)
    monkeypatch.setattr(fused_all, "track", counting_frames)
    runner.run_replicate(rep, cfg, root, state)
    assert len(opened) == 1, "compatible recovery must share one trajectory"
    assert processed == list(range(0, 8, 2)), "recovery must traverse the expected slice exactly once"
    assert completed_sibling.read_bytes() == sibling_bytes
    for entry, before in zip(targets, saved):
        assert entry.path.is_file()
        if entry.path.suffix == ".parquet":
            import pyarrow.parquet as pq
            # Compare saved values, not parquet serialization metadata.
            from io import BytesIO
            pd.testing.assert_frame_equal(
                pq.read_table(BytesIO(before)).to_pandas(),
                pd.read_parquet(entry.path),
                check_exact=False, rtol=1e-12, atol=1e-12,
            )
        else:
            from io import BytesIO
            np.testing.assert_array_equal(
                np.load(BytesIO(before), allow_pickle=False),
                np.load(entry.path, allow_pickle=False),
            )
    stages = json.loads((repout / "stage_timings.json").read_text())["stages"]
    assert [entry["stage"] for entry in stages] == ["fused-recovery", "summary"]


def test_complete_noop_does_not_reopen_trajectory(tmp_path, monkeypatch):
    rep, cfg, root, state, repout = _completed(tmp_path)

    def forbidden(*args, **kwargs):
        raise AssertionError("no-op run reopened a trajectory")

    verified = []
    original = runner._validate_analysis_output

    def counted(entry, *args, **kwargs):
        verified.append(entry.checkpoint_key)
        return original(entry, *args, **kwargs)

    monkeypatch.setattr(fused_all.mda, "Universe", forbidden)
    monkeypatch.setattr(runner, "_validate_analysis_output", counted)
    runner.run_replicate(rep, cfg, root, state)
    # Existing products require preflight and final audit, not a third read.
    from mdscope.registry import analysis_output_plan
    assert len(verified) == len(analysis_output_plan(cfg, repout)) * 2


@pytest.mark.parametrize("damage,match", [
    ("timestamp", "Trajectory times"),
    ("length", "frame count"),
])
def test_replay_rejects_mismatched_time_grid_before_publishing(
        tmp_path, monkeypatch, damage, match):
    rep, cfg, root, state, repout = _completed(tmp_path)
    target = fused_output(repout, "bilayer_normal").path
    sibling = fused_output(repout, "residue_membrane_distance", run_id="A").path
    sibling_bytes = sibling.read_bytes()
    target.unlink()
    summary = repout / "summary" / "timeseries.parquet"
    df = pd.read_parquet(summary)
    if damage == "timestamp":
        df.loc[df.index[-1], "time_ps"] += 1.
    else:
        df = df.iloc[:-1].copy()
    df.to_parquet(summary, index=False)
    marker = state.marker(runner._replicate_state_prefix(rep) + "-bilayer_normal")
    marker_bytes = marker.read_bytes()
    with pytest.raises(ValueError, match=match):
        runner.run_replicate(rep, cfg, root, state)
    assert not target.exists()
    assert sibling.read_bytes() == sibling_bytes
    assert marker.read_bytes() == marker_bytes
    assert not list(repout.glob(".fused-only-*")), "failed replay leaked staging"


@pytest.mark.parametrize("raises", [False, True])
def test_distance_finalize_preserves_frame_zero_metadata_and_restores_frame(
        tmp_path, monkeypatch, raises):
    """Metadata selection uses the old frame-zero semantics, including on error."""
    from types import SimpleNamespace

    from mdscope.pipeline import fused_distance

    class Timeline:
        frame = 5

        def __getitem__(self, index):
            self.frame = index
            return index

    observer = fused_distance.FusedDistances.__new__(fused_distance.FusedDistances)
    observer.u = SimpleNamespace(trajectory=Timeline())
    observer.cfg = {}
    observer.repout = tmp_path
    observer.cache = tmp_path / "cache"
    observer.runs = {"near": {"spec": {}}}
    observer.frames = 4
    observer.compute_seconds = 1.
    called = []

    def finalize_run(*args, **kwargs):
        called.append(kwargs["_universe"].trajectory.frame)
        if raises:
            raise RuntimeError("injected finalization error")
        return {"frames": 4}

    monkeypatch.setattr(fused_distance, "run_residue_membrane_distance", finalize_run)
    if raises:
        with pytest.raises(RuntimeError, match="injected finalization error"):
            observer.finalize(SimpleNamespace())
    else:
        assert observer.finalize(SimpleNamespace()) == {"near": {"frames": 4}}
    assert called == [0]
    assert observer.u.trajectory.frame == 5


def test_observer_only_mixed_families_replay_once_and_then_noop(tmp_path, monkeypatch):
    """No primary stream: still one fused pass, then no re-open on restart."""
    rep, cfg = make_project(tmp_path, 2)
    for name in (
        "alignment", "rmsd", "rmsf", "radius_of_gyration",
        "density_1d", "membrane", "salt_bridges", "ligand_contacts",
    ):
        cfg["analyses"][name]["enabled"] = False
    cfg["cache"]["aligned_coordinates"] = False
    cfg["progress"]["enabled"] = False
    root = tmp_path / "results"
    state = runner.RunState(root)

    def unexpected_primary(*args, **kwargs):
        raise AssertionError("observer-only task invoked primary stream")

    original_universe = fused_all.mda.Universe
    opened = []

    def counting_universe(*args, **kwargs):
        opened.append(args)
        return original_universe(*args, **kwargs)

    monkeypatch.setattr(runner, "run_streaming_pass", unexpected_primary)
    monkeypatch.setattr(fused_all.mda, "Universe", counting_universe)
    runner.run_replicate(rep, cfg, root, state)
    assert len(opened) == 1, "observer-only families must share one reader"
    repout = root / rep.system / rep.label
    for name in ("bilayer_normal", "solvation", "sdc4_mechanisms"):
        assert fused_output(repout, name).path.is_file()
    assert fused_output(
        repout, "residue_membrane_distance", run_id="A"
    ).path.is_file()
    runner.run_replicate(rep, cfg, root, state)
    assert len(opened) == 1, "complete observer-only restart must not reopen"
