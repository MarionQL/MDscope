"""Fused observer tests: MDAnalysis detector parity and chunk recovery."""
import json

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from mdscope.analysis.pipeline.fused_extensions import FusedExtensions
from mdscope.analysis.interactions.hydrogen_bonds import write_hydrogen_bonds


def _universe():
    u = mda.Universe.empty(
        8, n_residues=4, atom_resindex=[0, 0, 0, 1, 1, 1, 2, 3],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "H", "O", "N", "H", "O", "P", "P"])
    u.add_TopologyAttr("resnames", ["ALA", "ALA", "POPC", "POPE"])
    u.add_TopologyAttr("resids", [1, 2, 3, 4])
    u.add_TopologyAttr("bonds", [(0, 1), (3, 4)])
    base = np.array([
        [1., 1., 1.], [2., 1., 1.], [11.6, 1., 1.],
        [9., 1., 1.], [10., 1., 1.], [3.6, 1., 1.],
        [3., 3., 1.], [7., 3., 1.],
    ], dtype=np.float32)
    frames = np.array([base + [0., .05 * i, 0.] for i in range(5)])
    boxes = np.tile([50., 50., 50., 90., 90., 90.], (5, 1))
    u.load_new(frames, format=MemoryReader, dimensions=boxes)
    return u


def _cfg():
    return {
        "selections": {"protein": "resname ALA"},
        "trajectory": {"start": 0, "stop": 5, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": 5., "origin_ps": 10.},
        "analyses": {
            "hydrogen_bonds": {
                "enabled": True, "mode": "protein_protein",
                "donors_sel": "resname ALA and name N",
                "hydrogens_sel": "resname ALA and name H",
                "acceptors_sel": "resname ALA and name O",
                "d_h_cutoff_A": 1.2, "d_a_cutoff_A": 3.,
                "angle_cutoff_deg": 150., "update_selections": False,
            },
            "lipid_hydrogen_bonds": {"enabled": False},
            "dssp": {"enabled": False},
            "sasa": {"enabled": False},
            "local_environment": {
                "enabled": True, "protein_selection": "resname ALA",
                "lipid_selection": "resname POPC POPE",
                "motifs": {}, "contact_cutoff_A": 5.,
            },
        },
    }


def _run(root, resume=False):
    u, cfg = _universe(), _cfg()
    observer = FusedExtensions(u, cfg, root)
    if resume:
        for k in range(3):
            observer.update(u.trajectory[k], k)
            if k == 1:
                observer.flush(0)
        observer.flush(1)  # simulate a chunk written just before a crash
        u, cfg = _universe(), _cfg()
        observer = FusedExtensions(u, cfg, root)
        observer.prune_from(1)
        first, part = 2, 1
    else:
        observer.prune_from(0)
        first, part = 0, 0
    for k in range(first, 5):
        observer.update(u.trajectory[k], k)
        if (k + 1) % 2 == 0:
            observer.flush(part)
            part += 1
    if observer.frames % 2:
        observer.flush(part)
    observer.finalize()
    return u, cfg


def test_hbond_and_local_observers_restart_parity(tmp_path):
    clean, resumed = tmp_path / "clean", tmp_path / "resumed"
    u, cfg = _run(clean)
    _run(resumed, resume=True)
    for path in (
        "interactions/hydrogen_bonds/events.parquet",
        "interactions/hydrogen_bonds/count_timeseries.parquet",
        "extensions/local_environment/timeseries.parquet",
        "extensions/local_environment/residue_contacts.parquet",
    ):
        pd.testing.assert_frame_equal(
            pd.read_parquet(clean / path),
            pd.read_parquet(resumed / path),
            check_dtype=False,
        )
    # Independently traverse the same MemoryReader with the stock
    # MDAnalysis HydrogenBondAnalysis.run implementation.
    s = cfg["analyses"]["hydrogen_bonds"]
    hba = HydrogenBondAnalysis(
        u, donors_sel=s["donors_sel"], hydrogens_sel=s["hydrogens_sel"],
        acceptors_sel=s["acceptors_sel"],
        between=["resname ALA", "resname ALA"],
        d_h_cutoff=s["d_h_cutoff_A"], d_a_cutoff=s["d_a_cutoff_A"],
        d_h_a_angle_cutoff=s["angle_cutoff_deg"], update_selections=False,
    )
    hba.run(start=0, stop=5, step=1)
    expected = tmp_path / "reference"
    write_hydrogen_bonds(u, cfg, expected, hba.results.hbonds, hba.frames, hba.times)
    for name in ("events.parquet", "count_timeseries.parquet"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(clean / "interactions/hydrogen_bonds" / name),
            pd.read_parquet(expected / "interactions/hydrogen_bonds" / name),
            check_dtype=False, check_exact=False, rtol=1e-6, atol=1e-6,
        )
    assert len(pd.read_parquet(clean / "extensions/local_environment/timeseries.parquet")) == 5
    timing = json.loads((clean / "fused_extensions_timings.json").read_text())
    assert set(timing["enabled"]) == {"hydrogen_bonds", "local_environment"}
