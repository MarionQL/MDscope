"""Independent sparse contact and mass-selection equivalence checks."""
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance

from mdscope.analysis.biophysical_common import atom_contacts, residue_pairs, collapse_atom_pairs_to_residues
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


def test_sparse_contacts_empty_input_and_self_pairs():
    u = _universe()
    assert atom_contacts(u.atoms[[]], u.atoms, 2.0, u.dimensions).shape == (0, 2)
    contacts = atom_contacts(u.atoms[:2], u.atoms[:2], 2.0, u.dimensions)
    assert contacts.shape[1] == 2


def test_msd_residue_center_mass_weighting_and_fallback():
    u = _universe()
    # No mass metadata: legacy MSD falls back to residue center of geometry.
    geometric = u.atoms.center_of_geometry(compound="residues")
    np.testing.assert_allclose(_residue_centers(u.atoms), geometric)
    u.add_TopologyAttr("masses", [12., 16., 12., 14., 1.])
    mass_centers = u.atoms.center_of_mass(compound="residues")
    np.testing.assert_allclose(_residue_centers(u.atoms), mass_centers)
    assert not np.allclose(mass_centers[0], geometric[0])


def test_residue_pair_kernel_with_repeated_residue_numbers_across_chains():
    # Different topological residues can share the same human-readable resid.
    pairs = np.array([[0, 0], [1, 0], [2, 1], [2, 1]], dtype=int)
    left_resindices = np.array([7, 7, 11], dtype=int)
    right_resindices = np.array([3, 9], dtype=int)
    expected = np.array([[7, 3], [11, 9]], dtype=int)
    np.testing.assert_array_equal(
        collapse_atom_pairs_to_residues(pairs, left_resindices, right_resindices),
        expected,
    )
    assert collapse_atom_pairs_to_residues(
        np.empty((0, 2), dtype=int), left_resindices, right_resindices
    ).shape == (0, 2)


def test_contact_selection_contract_for_protein_membrane_and_membrane_only():
    u = mda.Universe.empty(
        8, n_residues=4, atom_resindex=[0, 0, 1, 1, 2, 2, 3, 3],
        n_segments=2, residue_segindex=[0, 0, 1, 1], trajectory=True,
    )
    u.add_TopologyAttr("names", ["C1", "O1"] * 4)
    u.add_TopologyAttr("resnames", ["ALA", "POPC", "POPE", "CHL1"])
    u.add_TopologyAttr("resids", [1, 1, 1, 1])
    u.add_TopologyAttr("segids", ["PROT", "MEMB"])
    u.atoms.positions = np.array([
        [0.5, 2, 2], [1.0, 2, 2], [9.5, 2, 2], [9., 2, 2],
        [8.8, 2, 2], [7.7, 2, 2], [19, 15, 15], [19.5, 15, 15],
    ], float)
    u.dimensions = [10., 20., 20., 90., 90., 90.]
    protein = u.select_atoms("resname ALA")
    membrane = u.select_atoms("resname POPC POPE CHL1")
    expected = capped_distance(
        protein.positions, membrane.positions,
        max_cutoff=2.5, box=u.dimensions, return_distances=False,
    )
    np.testing.assert_array_equal(
        residue_pairs(protein, membrane, 2.5, u.dimensions),
        collapse_atom_pairs_to_residues(expected, protein.resindices, membrane.resindices),
    )
    lipids = u.select_atoms("resname POPC")
    others = u.select_atoms("resname POPE CHL1")
    expected_lipids = capped_distance(
        lipids.positions, others.positions,
        max_cutoff=2.5, box=u.dimensions, return_distances=False,
    )
    np.testing.assert_array_equal(
        residue_pairs(lipids, others, 2.5, u.dimensions),
        collapse_atom_pairs_to_residues(expected_lipids, lipids.resindices, others.resindices),
    )
