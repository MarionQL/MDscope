"""Mixed raw/fused/specialized observers survive a checkpoint interruption."""
from copy import deepcopy
from types import SimpleNamespace
from dataclasses import replace
import hashlib
import shutil

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from mdscope.config import DEFAULTS
from mdscope.pipeline import streaming, runner
from validation_fixture import make_project


def _final_tables(root):
    return {str(p.relative_to(root)): pd.read_parquet(p) for p in root.rglob("*.parquet")
            if "cache" not in p.parts and "raw" not in p.parts}


@pytest.mark.parametrize("chunk,failure", [(1, "during_frame"), (2, "after_write"), (2, "before_checkpoint")])
def test_mixed_restart_and_offline_finalization(tmp_path, monkeypatch, chunk, failure):
    rep, cfg = make_project(tmp_path, chunk)
    cfg["analyses"]["sdc4_mechanisms"]["residue_occlusion"] = True
    cfg["statistics"]["families"]["profiles"] = True
    clean, resumed = tmp_path/"clean", tmp_path/"resumed"
    reference_cfg = deepcopy(cfg)
    reference_cfg["trajectory"]["chunk_frames"] = 3
    streaming.run_streaming_pass(rep, reference_cfg, clean, "reference")
    calls = 0
    target = {"during_frame": "compute_rg", "after_write": "_flush",
              "before_checkpoint": "_checkpoint_payload"}[failure]
    original = getattr(streaming, target)
    def interrupted(*a, **k):
        nonlocal calls
        calls += 1
        if failure == "after_write":
            result = original(*a, **k)
        if calls == (3 if failure == "during_frame" else 2):
            raise RuntimeError("injected mixed-observer interruption")
        return result if failure == "after_write" else original(*a, **k)
    monkeypatch.setattr(streaming, target, interrupted)
    with pytest.raises(RuntimeError, match="injected"):
        streaming.run_streaming_pass(rep, cfg, resumed, "restart")
    monkeypatch.setattr(streaming, target, original)
    streaming.run_streaming_pass(rep, cfg, resumed, "restart")

    # Summary consolidation and the actual project finalizer must be offline.
    rep.topology.unlink(); rep.trajectory.unlink()
    def forbidden(*a, **k):
        raise AssertionError("offline finalization attempted to reopen a trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    for root in (clean, resumed):
        runner.summarize_replicate(root/"test/rep1", cfg)
        runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=root,
                                state=SimpleNamespace(config_hash="mixed"), phase="statistics")
    expected, actual = _final_tables(clean), _final_tables(resumed)
    assert expected.keys() == actual.keys()
    assert any("density_1d/protein" in path for path in actual)
    assert any("density_1d/lipid" in path for path in actual)
    assert "test/rep1/extensions/sdc4_mechanisms/residue_occlusion.parquet" in actual
    profiles = actual["project_summary/statistics/profiles/replicate_profiles.parquet"]
    assert set(profiles.loc[profiles.profile.str.startswith("residue_occlusion__"), "profile"]) == {
        "residue_occlusion__protein", "residue_occlusion__hs",
        "residue_occlusion__lipid", "residue_occlusion__any",
    }
    assert any("solvation" in path for path in actual)
    for path in expected:
        if path.endswith("extensions/solvation/shell_stability.parquet"):
            # This diagnostic intentionally treats processing chunks as
            # time windows, so its sample count changes with chunk_frames.
            assert set(actual[path].columns) == set(expected[path].columns)
            assert actual[path].chunks_total.eq((4+chunk-1)//chunk).all()
            assert expected[path].chunks_total.eq(2).all()
            continue
        pd.testing.assert_frame_equal(actual[path], expected[path], check_exact=False, rtol=1e-9, atol=1e-9)
    for chain in ("A", "B"):
        relative = f"test/rep1/profiles/residue_membrane_distance/{chain}/distance_matrix_A.npy"
        np.testing.assert_array_equal(np.load(clean/relative), np.load(resumed/relative))
    table = actual["test/rep1/summary/timeseries.parquet"]
    assert len(table) == 4 and table.time_ps.is_unique and table.time_ps.is_monotonic_increasing
    np.testing.assert_allclose(table.time_ps, [0., 2., 4., 6.], rtol=1e-6, atol=1e-6)
    stats = actual["project_summary/statistics/time_series/replicate_statistics.parquet"]
    assert set(stats.replicate) == {"rep1"}
    relative = "test/rep1/extensions/sdc4_mechanisms/lipid_density_maps.npz"
    with np.load(clean/relative) as left, np.load(resumed/relative) as right:
        assert int(right["nframes"]) == 4
        assert left.files == right.files
        for key in left.files:
            np.testing.assert_array_equal(left[key], right[key])


@pytest.mark.parametrize("times", [[0.,0.,1.], [0.,2.,1.], [0.,np.inf], [0.,np.nan]])
def test_invalid_persisted_times_fail_before_archiving_existing_results(tmp_path, times):
    folder = tmp_path/"test/rep1/summary"
    folder.mkdir(parents=True)
    pd.DataFrame({"time_ps": times, "metric": np.ones(len(times))}).to_parquet(folder/"timeseries.parquet")
    previous = tmp_path/"project_summary/keep.txt"
    previous.parent.mkdir(); previous.write_text("existing valid results")
    cfg = {"systems": {"test": {"topology": "absent.pdb", "replicates": [{"label": "rep1", "trajectory": "absent.dcd"}]}}}
    with pytest.raises(RuntimeError, match="test/rep1: invalid time_ps"):
        runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=tmp_path,
                                state=SimpleNamespace(config_hash="bad"), phase="statistics")
    assert previous.read_text() == "existing valid results"
    assert not (tmp_path/".archives").exists()


def _bytes_snapshot(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


@pytest.mark.parametrize("relative,damage", [
    ("cache/fused_membrane/footprint_series/part-000000.parquet", "delete"),
    ("cache/sdc4_mechanisms/part-000000.parquet", "drop_frame"),
    ("cache/sdc4_mechanisms/part-000000.parquet", "duplicate_frame"),
    ("cache/sdc4_mechanisms/map-000000.npz", "corrupt"),
    ("profiles/density_1d/protein/frames/protein__part-000000.npz", "delete"),
    ("cache/fused_distance/A/part-000000.npz", "delete"),
    ("cache/solvation/counts-000000.parquet", "corrupt"),
])
def test_damaged_committed_observer_chunk_fails_before_rewriting_or_pruning(tmp_path, relative, damage):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path/"results"
    streaming.run_streaming_pass(rep, cfg, root, "mixed-integrity")
    repout = root/"test/rep1"
    path = repout/relative
    assert path.exists()
    if damage == "delete":
        path.unlink()
    elif damage == "corrupt":
        path.write_bytes(b"corrupt committed array")
    else:
        data = pd.read_parquet(path)
        data = data.iloc[:1] if damage == "drop_frame" else pd.concat([data.iloc[:1]]*2)
        data.to_parquet(path, index=False)
    # This is uncommitted and must survive an unsafe restart for diagnosis.
    orphan = repout/"cache/sdc4_mechanisms/part-999999.parquet"
    orphan.write_bytes(b"uncommitted, preserve on failure")
    before = _bytes_snapshot(repout)
    with pytest.raises(RuntimeError, match="committed artifact"):
        streaming.run_streaming_pass(rep, cfg, root, "mixed-integrity")
    assert _bytes_snapshot(repout) == before


def test_checkpoint_copied_to_other_replica_is_rejected(tmp_path):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path/"results"
    streaming.run_streaming_pass(rep, cfg, root, "same-project-hash")
    destination = root/"test/rep2"
    shutil.copytree(root/"test/rep1", destination, copy_function=shutil.copy2)
    before = _bytes_snapshot(destination)
    with pytest.raises(RuntimeError, match="replicate identity"):
        streaming.run_streaming_pass(replace(rep, label="rep2"), cfg, root, "same-project-hash")
    assert _bytes_snapshot(destination) == before


@pytest.mark.parametrize("damage", ["duplicate_chunk", "missing_checkpoint", "legacy", "counter_rollback",
                                  "fractional_counter", "accumulator_count", "group_count", "density_count", "counter_range"])
def test_checkpoint_integrity_fails_without_rewriting_existing_outputs(tmp_path, damage):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path/"results"
    streaming.run_streaming_pass(rep, cfg, root, "integrity")
    repout = root/"test/rep1"
    checkpoint = repout/"cache/stream_checkpoint.npz"
    if damage == "duplicate_chunk":
        source = repout/"cache/sdc4_mechanisms/part-000000.parquet"
        shutil.copy2(source, source.with_name("part-000000-copy.parquet"))
    elif damage == "missing_checkpoint":
        checkpoint.unlink()
        # Specialized parts alone must not be pruned as a fresh run.
        for p in (repout/"raw").glob("*part-*.parquet"):
            p.unlink()
    else:
        with np.load(checkpoint, allow_pickle=False) as cp:
            payload = {key: cp[key] for key in cp.files}
        if damage == "legacy":
            del payload["artifact_inventory_version"]
        elif damage == "fractional_counter":
            payload["processed"] = np.asarray(4.25)
        elif damage == "accumulator_count":
            payload["w_n"] = np.asarray(3)
        elif damage == "group_count":
            payload["g_A_n"] = np.asarray(3)
        elif damage == "density_count":
            payload["density_protein_0_frames"] = np.asarray(3)
        elif damage == "counter_range":
            payload["next_part"] = np.asarray(10**12)
        else:
            payload.update(processed=np.asarray(2), next_part=np.asarray(1))
        np.savez_compressed(checkpoint, **payload)
    before = _bytes_snapshot(repout)
    with pytest.raises(RuntimeError, match="artifact|Legacy|counter|accumulator"):
        streaming.run_streaming_pass(rep, cfg, root, "integrity")
    assert _bytes_snapshot(repout) == before


def test_offline_finalization_retains_unequal_replica_lengths_and_missing_values(tmp_path, monkeypatch):
    cfg = deepcopy(DEFAULTS)
    # This fixture only persists a scalar table. Its unrelated observers are disabled.
    for spec in cfg["analyses"].values():
        spec["enabled"] = False
    cfg["systems"] = {}
    for system, lengths, means in (("a", [3, 5], [0., 4.]), ("b", [4, 6], [2., 10.])):
        cfg["systems"][system] = {"topology": "unavailable.pdb", "replicates": []}
        for i, (length, mean) in enumerate(zip(lengths, means), 1):
            label = f"rep{i}"
            cfg["systems"][system]["replicates"].append({"label": label, "trajectory": "unavailable.dcd"})
            folder = tmp_path/system/label/"summary"
            folder.mkdir(parents=True)
            pd.DataFrame({
                "time_ps": np.arange(length)*10.,
                "value_A": np.arange(length)-(length-1)/2+mean,
                "unavailable_species": np.nan,
            }).to_parquet(folder/"timeseries.parquet")
    cfg["statistics"].update(bootstrap_draws=100, families={"profiles": False, "interactions": False, "solvation_rdfs": False})
    def forbidden(*a, **k):
        raise AssertionError("offline finalization reopened a trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=tmp_path,
                            state=SimpleNamespace(config_hash="unequal"), phase="statistics")
    folder = tmp_path/"project_summary/statistics/time_series"
    rep = pd.read_parquet(folder/"replicate_statistics.parquet")
    assert rep.set_index(["system", "replicate"]).n.to_dict() == {
        ("a", "rep1"): 3, ("a", "rep2"): 5, ("b", "rep1"): 4, ("b", "rep2"): 6}
    system = pd.read_parquet(folder/"system_statistics.parquet")
    assert system.set_index("system")["mean"].to_dict() == {"a": 2., "b": 6.}
    assert system.n_replicates.eq(2).all()
    assert set(system.metric) == {"value_A"}  # missing values did not become physical zeros
    for system in ("a", "b"):
        for label in ("rep1", "rep2"):
            assert pd.read_parquet(tmp_path/system/label/"summary/timeseries.parquet").unavailable_species.isna().all()


@pytest.mark.parametrize("damage", ["replicate_swap", "dropped_summary_frame"])
def test_finalization_rejects_checkpoint_summary_identity_or_frame_mismatch(tmp_path, damage):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path/"results"
    streaming.run_streaming_pass(rep, cfg, root, "summary-integrity")
    repout = root/"test/rep1"
    runner.summarize_replicate(repout, cfg)
    if damage == "replicate_swap":
        shutil.copytree(repout, root/"test/rep2", copy_function=shutil.copy2)
        cfg["systems"]["test"]["replicates"].append({"label": "rep2", "trajectory": str(rep.trajectory)})
        diagnostic = "replicate identity"
    else:
        path = repout/"summary/timeseries.parquet"
        pd.read_parquet(path).iloc[:-1].to_parquet(path, index=False)
        diagnostic = "frame count"
    previous = root/"project_summary/keep.txt"
    previous.parent.mkdir(); previous.write_text("previous results")
    before = _bytes_snapshot(root)
    with pytest.raises(RuntimeError, match=diagnostic):
        runner.finalize_project(tmp_path/"unused.yaml", cfg=cfg, root=root,
                                state=SimpleNamespace(config_hash="summary-integrity"), phase="statistics")
    assert _bytes_snapshot(root) == before


def test_completed_stream_rejects_missing_committed_chunk_without_reading_trajectory(tmp_path, monkeypatch):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    chunk = root / "test/rep1/raw/timeseries-part-000000.parquet"
    chunk.unlink()
    before = _bytes_snapshot(root)

    def forbidden(*args, **kwargs):
        raise AssertionError("completed stream reopened the trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    with pytest.raises(RuntimeError, match="missing committed timeseries parts"):
        runner.run_replicate(rep, cfg, root, state)
    assert _bytes_snapshot(root) == before


def test_offline_finalization_rejects_missing_enabled_fused_product(tmp_path):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    (root / "test/rep1/extensions/solvation/counts.parquet").unlink()
    before = _bytes_snapshot(root)
    with pytest.raises(RuntimeError, match="invalid required solvation artifact"):
        runner.finalize_project(tmp_path / "unused.yaml", cfg=cfg, root=root, state=state, phase="statistics")
    assert _bytes_snapshot(root) == before

@pytest.mark.parametrize("legacy_inventory", [False, True])
def test_older_completed_replica_is_reused_without_archiving_or_reading_trajectory(tmp_path, monkeypatch, legacy_inventory):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    repout = root / "test/rep1"
    (repout / ".config_hash").write_text("old-code-and-config")
    marker_path = state.marker(runner._replicate_state_prefix(rep) + "-stream")
    marker_path.write_text('{"status":"complete","config_hash":"old-code-and-config"}')
    checkpoint = repout / "cache/stream_checkpoint.npz"
    with np.load(checkpoint, allow_pickle=False) as data:
        payload = dict(data)
    if legacy_inventory:
        for key in ("artifact_inventory", "artifact_inventory_version", "replicate_identity"):
            del payload[key]
    np.savez_compressed(checkpoint, **payload, config_hash=np.asarray("old-code-and-config"))
    before = _bytes_snapshot(repout)

    def forbidden(*args, **kwargs):
        raise AssertionError("reuse opened the trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    runner.run_replicate(rep, cfg, root, state)
    assert _bytes_snapshot(repout) == before
    assert not (root / ".archives/replicates").exists()


def test_run_adds_missing_occlusion_and_reuses_completed_replica(tmp_path, monkeypatch):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    repout = root / "test/rep1"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    cfg["analyses"]["sdc4_mechanisms"]["residue_occlusion"] = True
    before = _bytes_snapshot(repout)
    with pytest.raises(RuntimeError, match="residue_occlusion.parquet"):
        runner.finalize_project(tmp_path / "unused.yaml", cfg=cfg, root=root,
                                state=state, phase="statistics")
    assert _bytes_snapshot(repout) == before

    def forbidden(*args, **kwargs):
        raise AssertionError("already completed work was called")

    monkeypatch.setattr(runner, "run_streaming_pass", forbidden)
    monkeypatch.setattr(runner, "run_sdc4_mechanisms", forbidden)
    monkeypatch.setattr(runner, "summarize_replicate", forbidden)
    runner.run_replicate(rep, cfg, root, state)

    path = repout / "extensions/sdc4_mechanisms/residue_occlusion.parquet"
    table = pd.read_parquet(path)
    assert table.time_ps.nunique() == 4
    assert table.groupby("time_ps").size().eq(24).all()
    after = _bytes_snapshot(repout)
    assert after.keys() - before.keys() == {str(path.relative_to(repout))}
    assert all(after[key] == value for key, value in before.items() if key != "stage_timings.json")

    monkeypatch.setattr(mda, "Universe", forbidden)
    runner.run_replicate(rep, cfg, root, state)
    assert _bytes_snapshot(repout) == after


def test_run_adds_membrane_density_and_keeps_the_completed_stream(tmp_path, monkeypatch):
    from mdscope.observables.density_1d import temporal_density_profiles

    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    repout = root / "test/rep1"
    state = runner.RunState(root)
    runner.run_replicate(rep, cfg, root, state)
    before = _bytes_snapshot(repout)
    cfg["analyses"]["density_1d"]["runs"]["membrane"] = {
        "membrane_center_selection": "resname POPC and name P",
        "half_width_A": 20.0,
        "bin_width_A": 1.0,
        "selections": {"phosphate": "resname POPC and name P"},
    }

    with pytest.raises(RuntimeError, match="missing membrane-relative density"):
        runner.finalize_project(tmp_path / "unused.yaml", cfg=cfg, root=root, state=state,
                                phase="statistics")
    assert _bytes_snapshot(repout) == before

    def forbidden(*args, **kwargs):
        raise AssertionError("completed observers were called again")

    monkeypatch.setattr(runner, "run_streaming_pass", forbidden)
    monkeypatch.setattr(runner, "summarize_replicate", forbidden)
    runner.run_replicate(rep, cfg, root, state)
    profiles = temporal_density_profiles(repout, cfg)
    relative = profiles["density_1d__membrane__phosphate"]
    np.testing.assert_allclose(relative["times"], [0., 2., 4., 6.], rtol=1e-6)
    assert relative["values"].shape == (4, 40)
    after = _bytes_snapshot(repout)
    assert all(after[key] == value for key, value in before.items() if key != "stage_timings.json")
    monkeypatch.setattr(mda, "Universe", forbidden)
    runner.run_replicate(rep, cfg, root, state)
    assert _bytes_snapshot(repout) == after


def test_completed_run_adds_only_footprint_and_rebuilds_summary(tmp_path, monkeypatch):
    rep, cfg = make_project(tmp_path, 2)
    root = tmp_path / "results"
    repout = root / "test/rep1"
    state = runner.RunState(root)
    cfg["analyses"]["lipid_enrichment"]["enabled"] = False
    runner.run_replicate(rep, cfg, root, state)
    before = _bytes_snapshot(repout)
    cfg["analyses"]["lipid_enrichment"]["enabled"] = True

    def forbidden(*args, **kwargs):
        raise AssertionError("completed streaming pass ran again")

    monkeypatch.setattr(runner, "run_streaming_pass", forbidden)
    runner.run_replicate(rep, cfg, root, state)
    new_path = repout / "extensions/lipid_footprint/timeseries.parquet"
    footprint = pd.read_parquet(new_path)
    summary = pd.read_parquet(repout / "summary/timeseries.parquet")
    pd.testing.assert_frame_equal(summary[footprint.columns], footprint)
    normal = "extensions/bilayer_normal/timeseries.parquet"
    assert _bytes_snapshot(repout)[normal] == before[normal]

    monkeypatch.setattr(mda, "Universe", forbidden)
    runner.run_replicate(rep, cfg, root, state)
