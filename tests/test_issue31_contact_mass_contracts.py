"""Independent sparse contact and mass-selection equivalence checks."""
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance

from mdscope.analysis.biophysical_common import atom_contacts, residue_pairs
from mdscope.analysis.extended.common import pairs
from mdscope.analysis.mechanisms.mechanistic_common import contact_count
from mdscope.analysis.observables.msd import _residue_centers


def _universe():
    u = mda.Universe.empty(5, n_residues=3, atom_resindex=[0, 0, 1, 2, 2], trajectory=True)
    u.add_TopologyAttr("names", ["C1", "O1", "C2", "N1", "H1"])
    u.atoms.positions = np.array([[0.5, 1, 1], [1.5, 1, 1], [9.5, 1, 1],
                                  [19, 19, 19], [19.5, 19, 19]], float)
    u.dimensions = [10, 20, 20, 90, 90, 90]
    return u


def test_sparse_contacts_and_residue_pairs_against_independent_reference():
    u = _universe()
    left, right = u.atoms[:2], u.atoms[2:]
    expected = capped_distance(left.positions, right.positions,
                               max_cutoff=2.1, box=u.dimensions,
                               return_distances=False)
    np.testing.assert_array_equal(atom_contacts(left, right, 2.1, u.dimensions), expected)
    assert contact_count(left, right, u.dimensions, cutoff=2.1) == len(expected)
    np.testing.assert_array_equal(
        residue_pairs(left, right, 2.1, u.dimensions),
        np.unique(np.stack([left.resindices[expected[:, 0]],
                            right.resindices[expected[:, 1]]], axis=1), axis=0),
    )
    extended_expected = capped_distance(left.positions.astype(np.float32),
                                        right.positions.astype(np.float32),
                                        max_cutoff=2.1, box=u.dimensions,
                                        return_distances=False)
    np.testing.assert_array_equal(pairs(left, right, 2.1, u.dimensions), extended_expected)


def test_sparse_contacts_empty_input_and_self_exclusion():
    u = _universe()
    assert pairs(u.atoms[[]], u.atoms, 2.0, u.dimensions).shape == (0, 2)
    results = pairs(u.atoms[:2], u.atoms[:2], 2.0, u.dimensions, exclude_self=True)
    assert all(i != j for i, j in results)


def test_msd_residue_center_mass_weighting_and_fallback():
    u = _universe()
    # No mass metadata: legacy MSD falls back to residue center of geometry.
    geometric = u.atoms.center_of_geometry(compound="residues")
    np.testing.assert_allclose(_residue_centers(u.atoms), geometric)
    u.add_TopologyAttr("masses", [12., 16., 12., 14., 1.])
    mass_centers = u.atoms.center_of_mass(compound="residues")
    np.testing.assert_allclose(_residue_centers(u.atoms), mass_centers)
    assert not np.allclose(mass_centers[0], geometric[0])
