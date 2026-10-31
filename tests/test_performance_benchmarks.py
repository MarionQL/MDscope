"""Numerically checked microbenchmarks; timings are reported, not flaky pass gates."""
from time import perf_counter

import MDAnalysis as mda
import numpy as np

from mdscope.analysis.observables.msd import _residue_centers
from mdscope.analysis.bilayer_normal import fit_local_normals


def _time(fn, repeat=4):
    durations = []
    for _ in range(repeat):
        start = perf_counter()
        fn()
        durations.append(perf_counter() - start)
    return min(durations)


def test_batch_lipid_centers_are_equivalent_and_report_speedup():
    n_lipids, atoms_per_lipid = 2048, 4
    u = mda.Universe.empty(
        n_lipids * atoms_per_lipid,
        n_residues=n_lipids,
        atom_resindex=np.repeat(np.arange(n_lipids), atoms_per_lipid),
        trajectory=True,
    )
    u.add_TopologyAttr("masses", np.full(len(u.atoms), 12.0))
    u.atoms.positions = np.random.default_rng(2468).uniform(
        5, 40, (len(u.atoms), 3)
    ).astype(np.float32)
    def baseline():
        return np.stack([res.atoms.center_of_mass() for res in u.residues])
    def batched():
        return _residue_centers(u.atoms)
    np.testing.assert_allclose(batched(), baseline(), rtol=1e-6, atol=1e-5)
    slow_s = _time(baseline)
    fast_s = _time(batched)
    print(
        f"BENCHMARK lipid_COM n={n_lipids} reference={slow_s:.6f}s "
        f"batched={fast_s:.6f}s measured_speedup={slow_s / max(fast_s, 1e-12):.2f}x",
        flush=True,
    )


def test_local_normals_large_sparse_neighbor_smoke():
    rng = np.random.default_rng(37)
    n = 1024
    xy = rng.uniform(0., 100., (n, 2))
    z = np.where(np.arange(n) < n // 2, 25., 75.)
    xyz = np.column_stack((xy, z)).astype(np.float32)
    upper = np.arange(n) >= n // 2
    box = np.asarray([100., 100., 100., 90., 90., 90.])
    def compute():
        return fit_local_normals(xyz, upper, box, radius=15., min_neighbors=5)
    angles, neighbors = compute()
    assert angles.shape == (n,) and neighbors.shape == (n,)
    assert np.nanmax(np.abs(angles)) < 1e-4  # planar leaflets
    seconds = _time(compute, repeat=3)
    print(
        f"BENCHMARK local_normals n={n} radius=15A time={seconds:.6f}s",
        flush=True,
    )
