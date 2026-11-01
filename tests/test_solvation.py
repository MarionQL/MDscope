"""Independent tiny-system references for solvation counting, PBC and kinetics."""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.observables.solvation import (SolvationObserver, contact_episodes,
                                           rdf_from_counts)


def _universe():
    u = mda.Universe.empty(5, n_residues=4, atom_resindex=[0, 0, 1, 2, 3],
                           trajectory=True)
    u.add_TopologyAttr("names", ["OE1", "OE2", "OW", "OW", "SOD"])
    u.add_TopologyAttr("resnames", ["GLU", "SOL", "SOL", "SOD"])
    u.add_TopologyAttr("resids", [7, 1, 2, 3])
    xyz = np.array([[1, 1, 1], [1.5, 1, 1], [3, 1, 1],
                    [19.5, 1, 1], [3, 1, 2]], np.float32)
    shifted = xyz.copy()
    shifted[2] = [12, 1, 1]
    u.load_new(np.array([xyz, shifted]), format=MemoryReader,
               dimensions=np.tile([20, 20, 20, 90, 90, 90], (2, 1)))
    return u


def test_bulk_rdf_formula():
    edges = np.array([0., 1., 2.])
    expected = np.array([3., 21.])
    np.testing.assert_allclose(rdf_from_counts(expected, edges, 3 / (4*np.pi/3)),
                               [1, 1])


def test_contact_episodes_censor_and_gaps():
    events = [(0, "water", "glu", 4), (1, "water", "glu", 4),
              (3, "water", "glu", 4)]
    out = contact_episodes(events, [(i, i*5.) for i in range(4)])
    assert len(out) == 2
    assert out.duration_ps.tolist() == [5., 0.]
    assert out.left_censored.tolist() == [True, False]
    assert out.right_censored.tolist() == [False, True]


def test_periodic_unique_water_count_and_restart_histograms(tmp_path):
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "kinetics",
        "targets": {"protein": {"carboxylate": "resname GLU"}},
        "water_selection": "resname SOL and name OW",
        "ion_selection": "resname SOD",
        "search_radius_A": 5., "rdf_bin_width_A": 0.5,
        "shell_cutoffs_A": {"water": 3.5, "ions": 4.}
    }}}
    u = _universe()
    obs = SolvationObserver(u, cfg, tmp_path)
    obs.update(u.trajectory[0], 0, 0.)
    obs.flush(0)
    obs.update(u.trajectory[1], 1, 5.)
    obs.flush(1)
    obs.finalize()
    root = tmp_path / "extensions" / "solvation"
    counts = pd.read_parquet(root / "counts.parquet")
    water = counts[counts.solvent == "water"]
    assert water.n_solvent.tolist() == [2, 1]
    ions = counts[counts.solvent == "ions"]
    assert ions.n_solvent.tolist() == [1, 1]
    assert len(pd.read_parquet(root / "rdf.parquet")) == 20
    events = pd.read_parquet(root / "events.parquet")
    assert len(events[events.solvent == "water"]) == 3
    assert (root / "episodes.parquet").exists()
    assert (root / "metadata.json").exists()


def test_local_grid_uses_same_water_candidates(tmp_path):
    u = _universe()
    cfg = {"analyses": {"solvation": {
        "enabled": True, "preset": "standard",
        "targets": {"protein": {"oxygen": "resname GLU and name OE1 OE2"}},
        "water_selection": "resname SOL and name OW",
        "grid_frames": {"protein/oxygen/0": "index 0 1 4"},
        "grid_width_A": 10., "grid_spacing_A": 1.,
        "search_radius_A": 8.,
    }}}
    observer = SolvationObserver(u, cfg, tmp_path)
    observer.update(u.trajectory[0], 0, 0.)
    observer.flush(0)
    observer.finalize()
    density = pd.read_parquet(tmp_path/"extensions/solvation/water_density_3d.parquet")
    assert density["count"].sum() == 2
    assert np.all(np.isfinite(density["g_bulk"]))
