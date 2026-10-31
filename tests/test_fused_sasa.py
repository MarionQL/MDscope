"""SASA streaming serialization parity with a deterministic FreeSASA test double.

The physical SASA algorithm is FreeSASA's responsibility. This test verifies
identical atom coordinates, per-residue reduction, time handling and output
schemas between the independent and fused frame iterators.
"""
import sys
import types

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.pipeline.fused_extensions import FusedExtensions
from mdscope.analysis.observables.sasa import _run_sasa


class FakeStructure:
    def __init__(self):
        self.coordinates = []

    def addAtom(self, name, resname, resid, segid, x, y, z):
        self.coordinates.append((x, y, z))


class FakeParameters:
    def setAlgorithm(self, value): pass
    def setNPoints(self, value): pass
    def setNSlices(self, value): pass
    def setProbeRadius(self, value): pass
    def setNThreads(self, value): pass


class FakeResult:
    def __init__(self, structure):
        self.values = np.linalg.norm(structure.coordinates, axis=1)

    def atomArea(self, index):
        return self.values[index]


def test_sasa_fused_vs_independent(monkeypatch, tmp_path):
    fake = types.SimpleNamespace(
        Structure=FakeStructure, Parameters=FakeParameters,
        ShrakeRupley=1, LeeRichards=2,
        calc=lambda structure, params: FakeResult(structure),
    )
    monkeypatch.setitem(sys.modules, "freesasa", fake)
    u = mda.Universe.empty(
        6, n_residues=2, atom_resindex=[0, 0, 0, 1, 1, 1],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["N", "CA", "C"] * 2)
    u.add_TopologyAttr("resnames", ["ALA", "ALA"])
    u.add_TopologyAttr("resids", [1, 2])
    u.add_TopologyAttr("segids", ["SYSTEM"])
    first = np.arange(18, dtype=np.float32).reshape(6, 3) + 1
    frames = np.stack([first + i for i in range(4)])
    u.load_new(frames, format=MemoryReader,
               dimensions=np.tile([100., 100., 100., 90., 90., 90.], (4, 1)))
    cfg = {
        "trajectory": {"start": 0, "stop": 4, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": 2., "origin_ps": 0.},
        "selections": {"protein": "resname ALA"},
        "analyses": {
            "sasa": {"enabled": True, "selection": "resname ALA"},
            "dssp": {"enabled": False},
            "hydrogen_bonds": {"enabled": False},
            "lipid_hydrogen_bonds": {"enabled": False},
            "local_environment": {"enabled": False},
        },
        "progress": {"enabled": False},
    }
    independent = tmp_path / "independent"
    _run_sasa(u, cfg, independent)
    fused = tmp_path / "fused"
    u.trajectory[0]
    observer = FusedExtensions(u, cfg, fused)
    for i in range(4):
        observer.update(u.trajectory[i], i)
        if (i + 1) % 2 == 0:
            observer.flush(i // 2)
    observer.finalize()
    pd.testing.assert_frame_equal(
        pd.read_parquet(independent / "profiles/sasa/residue_sasa.parquet"),
        pd.read_parquet(fused / "profiles/sasa/residue_sasa.parquet"),
        check_dtype=False, rtol=1e-12, atol=1e-12,
    )
