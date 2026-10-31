"""Parity of DSSP's native per-frame implementation against .run()."""
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.dssp import _run_dssp
from mdscope.analysis.pipeline.fused_extensions import FusedExtensions


def test_dssp_fused_vs_independent(tmp_path):
    nres = 7
    u = mda.Universe.empty(
        4 * nres, n_residues=nres,
        atom_resindex=np.repeat(np.arange(nres), 4), trajectory=True
    )
    u.add_TopologyAttr("names", ["N", "CA", "C", "O"] * nres)
    u.add_TopologyAttr("resnames", ["ALA"] * nres)
    u.add_TopologyAttr("resids", list(range(1, nres + 1)))
    coords = []
    for i in range(nres):
        # Synthetic extended backbone with realistic bond lengths.
        x = i * 3.8
        coords.extend([
            [x, 0., 0.], [x + 1.46, 0.5, 0.],
            [x + 2.90, 0., 0.], [x + 3.35, -1.05, 0.],
        ])
    coords = np.asarray(coords, dtype=np.float32)
    frames = np.stack([coords, coords + [0., .01, 0.], coords + [0., .02, 0.]])
    u.load_new(frames, format=MemoryReader,
               dimensions=np.tile([100., 100., 100., 90., 90., 90.], (3, 1)))
    cfg = {
        "trajectory": {"start": 0, "stop": 3, "stride": 1, "chunk_frames": 2},
        "time": {"dt_ps": 10., "origin_ps": 0.},
        "selections": {"protein": "resname ALA"},
        "analyses": {
            "dssp": {"enabled": True, "selection": "resname ALA", "guess_hydrogens": True},
            "sasa": {"enabled": False},
            "hydrogen_bonds": {"enabled": False},
            "lipid_hydrogen_bonds": {"enabled": False},
            "local_environment": {"enabled": False},
        },
        "progress": {"enabled": False},
    }
    expected = tmp_path / "independent"
    _run_dssp(u, cfg, expected)
    u.trajectory[0]
    actual = tmp_path / "fused"
    observer = FusedExtensions(u, cfg, actual)
    for i in range(3):
        observer.update(u.trajectory[i], i)
        if i == 1:
            observer.flush(0)
    observer.flush(1)
    observer.finalize()
    pd.testing.assert_frame_equal(
        pd.read_parquet(expected / "profiles/dssp/assignments.parquet"),
        pd.read_parquet(actual / "profiles/dssp/assignments.parquet"),
        check_dtype=False,
    )
