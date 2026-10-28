"""Optional second trajectory for raw MSD: independent coordinates and strict sync."""
from pathlib import Path
from types import SimpleNamespace

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.config import iter_replicates
from mdscope.analysis.core import trajectory_sources
from mdscope.analysis.membrane import membrane_metrics
from mdscope.analysis.msd import initialize_msd


def _universe(translation):
    u = mda.Universe.empty(
        4, n_residues=4, atom_resindex=np.arange(4), trajectory=True,
    )
    u.add_TopologyAttr("names", ["P"] * 4)
    u.add_TopologyAttr("resnames", ["POPC"] * 4)
    u.add_TopologyAttr("resids", [1, 2, 3, 4])
    u.add_TopologyAttr("masses", [31.] * 4)
    base = np.asarray([[4, 1, 8], [7, 1, 8], [4, 3, 2], [7, 3, 2]], np.float32)
    frames = np.asarray([base, base + [translation, 0., 0.]], np.float32)
    boxes = np.tile([10., 10., 10., 90., 90., 90.], (2, 1))
    u.load_new(frames, format=MemoryReader, dimensions=boxes)
    return u


def _membrane():
    state = SimpleNamespace(
        enabled=True, msd=True, apl=False, thickness=False, ai=2,
        upper_resix=np.asarray([0, 1]), lower_resix=np.asarray([2, 3]),
    )
    return state


def test_independent_msd_reader_uses_original_frames(monkeypatch):
    centered_main = _universe(0.)
    uncentered_raw = _universe(2.)
    state = _membrane()
    spec = {"headgroup_selection": "all", "msd": {"enabled": True}}
    initialize_msd(centered_main, state, spec)
    rep = SimpleNamespace(topology=Path("unused.tpr"), msd_trajectory=Path("raw.xtc"))
    monkeypatch.setattr(trajectory_sources.mda, "Universe", lambda *args: uncentered_raw)
    raw = trajectory_sources.open_msd_universe(
        rep, {"analyses": {"membrane": spec}}, state, centered_main, 0,
    )
    assert state.msd_atoms.universe is raw
    values = []
    for _ in trajectory_sources.paired_frames(
        centered_main.trajectory[:], raw.trajectory[:]
    ):
        values.append(membrane_metrics(centered_main, state, raw))
    np.testing.assert_allclose(values[0]["membrane_msd_lateral_A2"], 0.)
    np.testing.assert_allclose(values[1]["membrane_msd_lateral_A2"], 4.)
    np.testing.assert_allclose(centered_main.atoms.positions[:, 0], [4, 7, 4, 7])


def test_msd_secondary_stream_rejects_missing_frames_and_timestamp_mismatch():
    ts = lambda t: SimpleNamespace(time=t)
    with pytest.raises(ValueError, match="sampled frame counts"):
        list(trajectory_sources.paired_frames([ts(0), ts(1)], [ts(0)]))
    with pytest.raises(ValueError, match="time synchronized"):
        list(trajectory_sources.paired_frames([ts(0), ts(1)], [ts(0), ts(2)]))


def test_replicate_raw_msd_trajectory_path_resolution(tmp_path):
    cfg = {"systems": {"membrane": {
        "topology": "topol.tpr",
        "replicates": [{
            "trajectory": "processed.xtc",
            "msd_trajectory": "unmodified.xtc",
            "label": "rep1",
        }],
    }}}
    paths = list(iter_replicates(cfg, tmp_path / "cfg.yaml"))
    assert len(paths) == 1
    assert paths[0].trajectory == tmp_path / "processed.xtc"
    assert paths[0].msd_trajectory == tmp_path / "unmodified.xtc"


def test_raw_reader_rejects_mismatched_frame_count(monkeypatch):
    main = _universe(0.)
    raw = _universe(2.)
    # A reader that has fewer frames is never silently aligned by index.
    raw.load_new(np.zeros((1, 4, 3), np.float32),
                 format=MemoryReader,
                 dimensions=np.asarray([[10., 10., 10., 90., 90., 90.]]))
    state = _membrane()
    spec = {"headgroup_selection": "all", "msd": {"enabled": True}}
    initialize_msd(main, state, spec)
    monkeypatch.setattr(trajectory_sources.mda, "Universe", lambda *args: raw)
    rep = SimpleNamespace(topology=Path("unused.tpr"), msd_trajectory=Path("raw.xtc"))
    with pytest.raises(ValueError, match="different frame counts"):
        trajectory_sources.open_msd_universe(
            rep, {"analyses": {"membrane": spec}}, state, main, 0,
        )
