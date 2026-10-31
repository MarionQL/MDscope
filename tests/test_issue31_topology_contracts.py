"""Topology edge cases without real trajectories for scientific geometry."""
import numpy as np
import pytest
import MDAnalysis as mda
from mdscope.observables.msd import _residue_centers
from mdscope.core.atom_selection import heavy_atom_mask


def _atoms():
    u = mda.Universe.empty(
        n_atoms=5, n_residues=3, atom_resindex=[2, 0, 2, 1, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["C1", "H1", "C2", "O1", "H2"])
    u.atoms.positions = np.array([
        [1., 2., 3.], [10., 0., 0.], [3., 2., 3.],
        [5., 5., 5.], [12., 0., 0.],
    ])
    return u


def test_residue_centers_reordered_topology_indices_and_missing_masses():
    u = _atoms()
    atoms = u.atoms
    ids, inverse = np.unique(atoms.resindices, return_inverse=True)
    geometric = np.column_stack([
        np.bincount(inverse, weights=atoms.positions[:, i], minlength=len(ids))
        / np.bincount(inverse, minlength=len(ids)) for i in range(3)
    ])
    np.testing.assert_array_equal(ids, [0, 1, 2])
    np.testing.assert_allclose(geometric, [[11., 0., 0.], [5., 5., 5.], [2., 2., 3.]])
    # MSD's fallback is center of geometry, but preserves MDAnalysis
    # residue-group order, not arbitrary selected-atom order.
    np.testing.assert_allclose(_residue_centers(atoms), geometric)


def test_heavy_selection_policies_remain_distinct_without_elements():
    u = _atoms()
    assert heavy_atom_mask(u.atoms, policy="extended").tolist() == [
        True, False, True, True, False
    ]
    assert heavy_atom_mask(u.atoms, policy="mechanistic").tolist() == [
        True, False, True, True, False
    ]
