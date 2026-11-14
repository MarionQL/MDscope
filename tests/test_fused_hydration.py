"""Issue #155: fused HS, protein and membrane hydration with unique-water unions."""
from copy import deepcopy

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.mechanisms.hs_hydration import OUTPUT, temporal_hs_hydration_profiles
from mdscope.pipeline.fused_extensions import FusedExtensions, enabled_names
from mdscope.registry import fused_output, required_artifacts


def _universe():
    # The two waters are separate residues; one straddles a periodic boundary.
    u = mda.Universe.empty(9, n_residues=9,
                           atom_resindex=np.arange(9), trajectory=True)
    u.add_TopologyAttr("names", ["OG", "O", "O1", "O2", "O3", "O4", "O5", "OW", "OW"])
    u.add_TopologyAttr("resnames",
                       ["SER", "ALA", "BXYL", "BGAL", "BGAL",
                        "POPC", "POPE", "TIP3", "TIP3"])
    u.add_TopologyAttr("resids", [39, 40, 99, 41, 15, 5, 6, 100, 101])
    u.add_TopologyAttr("chainIDs", ["PROA"] * 9)
    u.add_bonds([(0, 2), (2, 3), (3, 4)])
    base = np.array([[1., 0., 0.], [4., 0., 0.], [1., 0., 0.],
                     [2., 0., 0.], [3., 0., 0.], [2., 0., 0.],
                     [20., 0., 0.], [1.2, 0., 0.], [29.9, 0., 0.]],
                    dtype=np.float32)
    second = base.copy()
    second[-1] = [23., 0., 0.]
    boxes = np.tile([30., 30., 30., 90., 90., 90.], (2, 1))
    u.load_new(np.stack([base, second]), format=MemoryReader, dimensions=boxes)
    return u


def _cfg(*, general=True):
    return {
        "selections": {"protein": "resname SER ALA"},
        "trajectory": {"start": 0, "stop": 2, "stride": 1, "chunk_frames": 1},
        "time": {"origin_ps": 10., "dt_ps": 5.},
        "analyses": {
            "sdc4_mechanisms": {
                "hs_sugar_water_cutoff_A": 1.4,
                "water_oxygen_selection": "name OW",
                "hs_auto": {"glycan_selection": "resname BXYL BGAL",
                            "attachment_selection": "name OG"},
            },
            "hydration": {
                "enabled": general, "cutoff_A": 1.4,
                "water_oxygen_selection": "name OW",
                "protein_selection": "resname SER ALA",
                "membrane_selection": "resname POPC POPE",
            },
        },
    }


def _run(root, cfg, *, resume=False):
    universe = _universe()
    fused = FusedExtensions(universe, cfg, root)
    if resume:
        fused.update(universe.trajectory[0], 0)
        fused.flush(0)
        fused.update(universe.trajectory[1], 1)
        fused.flush(1)  # written after the last durable checkpoint
        universe = _universe()
        fused = FusedExtensions(universe, cfg, root)
        fused.prune_from(1)
        start = 1
    else:
        fused.prune_from(0)
        start = 0
    for i in range(start, 2):
        fused.update(universe.trajectory[i], i)
        fused.flush(i)
    fused.finalize()


def test_hydration_fused_unique_water_aggregation_and_time_grid(tmp_path):
    cfg = _cfg()
    assert enabled_names(cfg) == ("hydration",)
    _run(tmp_path, cfg)
    sugar = pd.read_parquet(tmp_path / OUTPUT)
    assert sugar.sugar_resid.tolist() == [99, 41, 15] * 2
    assert sugar.water_count.tolist() == [2, 1, 0, 1, 1, 0]
    assert sugar.time_ps.tolist() == [10.] * 3 + [15.] * 3

    counts = pd.read_parquet(tmp_path / "extensions/hydration/counts.parquet")
    aggregate = pd.read_parquet(tmp_path / "extensions/hydration/aggregate_counts.parquet")
    assert len(counts) == 2 * (3 + 2 + 2)
    for kind, level, group, expected in [
        ("hs", "chain", "PROA_S39", [2, 1]),
        ("hs", "total", "all", [2, 1]),
        ("protein", "total", "all", [2, 1]),
        ("membrane", "total", "all", [1, 1]),
        ("membrane", "lipid_type", "POPE", [0, 0]),
    ]:
        rows = aggregate.loc[
            aggregate.compartment.eq(kind) & aggregate.level.eq(level)
            & aggregate.group_id.eq(group)]
        assert rows.water_count.tolist() == expected
    # A shared water is counted once per sugar, but once in a chain union.
    assert sugar.water_count.iloc[:3].sum() > aggregate.loc[
        aggregate.compartment.eq("hs") & aggregate.level.eq("total"),
        "water_count"].iloc[0]
    meta = pd.read_parquet(tmp_path / "extensions/hydration/target_metadata.parquet")
    assert len(meta) == 7
    assert set(meta.compartment) == {"hs", "protein", "membrane"}
    times = np.array([10., 15.])
    profiles = temporal_hs_hydration_profiles(tmp_path, cfg, times)
    np.testing.assert_array_equal(
        profiles["hs_water_by_position__PROA_S39"]["values"],
        [[2., 1., 0.], [1., 1., 0.]])

    from mdscope.statistics.time_series import collect_time_series
    series = collect_time_series(tmp_path)
    assert series[("hydration/hs/chain/PROA_S39", "unique_water_count")].value.tolist() == [2, 1]
    assert series[("hydration/protein/total/all", "unique_water_count")].value.tolist() == [2, 1]
    assert series[("hydration/membrane/total/all", "unique_water_count")].value.tolist() == [1, 1]
    assert np.isfinite(series[("hydration/protein/total/all",
                               "unique_water_per_site")].value).all()

    entry = fused_output(tmp_path, "hydration")
    assert all(path.is_file() for path in required_artifacts(entry, cfg))


def test_hydration_restart_replays_only_uncommitted_chunks(tmp_path):
    cfg = _cfg()
    clean = tmp_path / "clean"
    resumed = tmp_path / "resumed"
    _run(clean, cfg)
    _run(resumed, cfg, resume=True)
    for file in [OUTPUT,
                 "extensions/hydration/counts.parquet",
                 "extensions/hydration/aggregate_counts.parquet",
                 "extensions/hydration/target_metadata.parquet"]:
        pd.testing.assert_frame_equal(
            pd.read_parquet(clean / file), pd.read_parquet(resumed / file))


def test_legacy_hs_only_does_not_require_general_solvation(tmp_path):
    cfg = _cfg(general=False)
    _run(tmp_path, cfg)
    assert len(pd.read_parquet(tmp_path / OUTPUT)) == 6
    counts = pd.read_parquet(tmp_path / "extensions/hydration/counts.parquet")
    assert set(counts.compartment) == {"hs"}


def test_hydration_water_selection_and_box_contract(tmp_path):
    cfg = deepcopy(_cfg())
    u = _universe()
    cfg["analyses"]["hydration"]["water_oxygen_selection"] = "name O1"
    # Different water oxygen populations must not be silently substituted.
    import pytest
    with pytest.raises(ValueError, match="water selections differ"):
        FusedExtensions(u, cfg, tmp_path / "mismatch")
    cfg["analyses"]["sdc4_mechanisms"]["water_oxygen_selection"] = "name O1"
    with pytest.raises(ValueError, match="targets overlap water selection"):
        FusedExtensions(u, cfg, tmp_path / "overlap")
    cfg = _cfg()
    cfg["analyses"]["hydration"]["membrane_selection"] = "resname POPC BGAL"
    with pytest.raises(ValueError, match="target selections overlap"):
        FusedExtensions(u, cfg, tmp_path / "membrane_overlap")
    cfg = _cfg()
    observer = FusedExtensions(u, cfg, tmp_path / "bad_box")
    ts = u.trajectory[0]
    ts.dimensions = [0., 30., 30., 90., 90., 90.]
    with pytest.raises(ValueError, match="periodic box"):
        observer.update(ts, 0)
