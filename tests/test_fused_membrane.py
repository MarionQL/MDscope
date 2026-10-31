"""Shared-trajectory observer tests including crash/replay semantics."""
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.fused_membrane import FusedMembrane, enabled_names
from mdscope.analysis.core.lipid_geometry import anchor_groups


def _universe():
    u = mda.Universe.empty(
        10, n_residues=10, atom_resindex=np.arange(10), trajectory=True
    )
    u.add_TopologyAttr("names", ["P"] * 8 + ["CA", "CA"])
    u.add_TopologyAttr("resnames", ["POPC"] * 8 + ["ALA"] * 2)
    u.add_TopologyAttr("resids", list(range(1, 11)))
    x = np.array([10., 25., 10., 25., 10., 25., 10., 25.])
    y = np.array([10., 10., 25., 25., 10., 10., 25., 25.])
    z = np.array([20., 20., 20., 20., 60., 60., 60., 60.])
    first = np.column_stack([x, y, z])
    protein = np.array([[99., 15., 40.], [1., 15., 40.]])
    coordinates = []
    for i in range(5):
        frame = np.vstack([first + [i * .1, 0, 0], protein])
        frame[8, 0] -= 100. if i % 2 else 0.
        coordinates.append(frame)
    boxes = np.tile([100., 100., 100., 90., 90., 90.], (5, 1))
    u.load_new(np.asarray(coordinates, dtype=np.float32),
               format=MemoryReader, dimensions=boxes)
    return u


def _cfg():
    selections = {"POPC": "resname POPC and name P"}
    return {
        "trajectory": {"start": 0, "stop": 5, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": None, "origin_ps": 0.},
        "selections": {"protein": "name CA"},
        "analyses": {
            "membrane_order": {"enabled": False},
            "bilayer_normal": {
                "enabled": True, "axis": "z",
                "species_selections": selections,
                "neighbor_radius_A": 30., "min_neighbors": 3,
                "save_lipid_angles": True,
            },
            "lipid_enrichment": {
                "enabled": True, "axis": "z",
                "species_selections": selections,
                "protein_selection": "name CA",
                "inner_radius_A": 0., "outer_radius_A": 30.,
            },
        },
    }


def _execute(u, cfg, root, resume=False):
    u.trajectory[0]
    observer = FusedMembrane(u, cfg, root)
    if resume:
        for k in range(3):
            observer.update(u.trajectory[k], k)
            if k == 1:
                observer.flush(0)
        observer.flush(1)  # simulate uncommitted part written before crash
        del observer
        u.trajectory[0]
        observer = FusedMembrane(u, cfg, root)
        observer.prune_from(1)
        start = 2
        next_part = 1
    else:
        observer.prune_from(0)
        start = 0
        next_part = 0
    for k in range(start, 5):
        observer.update(u.trajectory[k], k)
        if (k + 1) % 2 == 0:
            observer.flush(next_part)
            next_part += 1
    if observer.frame_count % 2:
        observer.flush(next_part)
    return observer.finalize()


def test_shared_anchor_and_checkpoint_replay(tmp_path):
    cfg = _cfg()
    assert enabled_names(cfg) == ("bilayer_normal", "lipid_enrichment")
    root_a, root_b = tmp_path / "clean", tmp_path / "resumed"
    result_a = _execute(_universe(), cfg, root_a)
    result_b = _execute(_universe(), cfg, root_b, resume=True)
    assert result_a["bilayer_normal"]["frames"] == 5
    assert result_b["lipid_enrichment"]["frames"] == 5
    for path in (
        "extensions/bilayer_normal/timeseries.parquet",
        "extensions/lipid_enrichment/timeseries.parquet",
        "extensions/lipid_enrichment/counts.parquet",
        "extensions/lipid_enrichment/selection_audit.parquet",
    ):
        pd.testing.assert_frame_equal(
            pd.read_parquet(root_a / path), pd.read_parquet(root_b / path),
            check_dtype=False,
        )
    for root in (root_a, root_b):
        timing = json.loads((root / "fused_membrane_timings.json").read_text())
        assert timing["shared_anchor_mappings"] == 1
        angles = sorted((root / "extensions/bilayer_normal/per_lipid").glob("*.parquet"))
        assert len(angles) == 3
        assert sum(len(pd.read_parquet(path)) for path in angles) == 5 * 8
