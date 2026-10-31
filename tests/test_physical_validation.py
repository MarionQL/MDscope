"""Scientific-regression checks for PBC, leaflet tracking and glycolipid identity."""
import numpy as np
import MDAnalysis as mda

from mdscope.analysis.core.lipid_geometry import (
    anchor_groups, aligned_leaflets, classify_leaflets, periodic_center,
)
from mdscope.analysis.observables.lipid_enrichment import enrichment_frame
from mdscope.analysis.observables.bilayer_normal import fit_local_normals
from mdscope.analysis.lipid_hydrogen_bonds import lipid_molecule_map, tabulate_events
from mdscope.analysis.interaction_tables import INTERACTION_KEY_COLUMNS


BOX = np.array([100., 100., 100., 90., 90., 90.])


def test_protein_periodic_center_invariant_to_independent_wrapping():
    protein = np.array([[99., 30., 50.], [1., 31., 50.], [98., 29., 51.]])
    baseline = periodic_center(protein, BOX)
    shifted = protein.copy()
    shifted[0, 0] -= 100.
    shifted[1, 0] += 100.
    shifted[2, 1] += 100.
    np.testing.assert_allclose(periodic_center(shifted, BOX), baseline, atol=1e-10)
    # Naive center incorrectly places the protein at x~66 rather than x~99.
    assert abs(protein.mean(axis=0)[0] - baseline[0]) > 25


def test_enrichment_invariant_to_lipid_and_protein_wrapping():
    species = np.array(["POPC", "POPE", "POPC", "POPE", "POPC", "POPE"])
    upper = np.array([1, 1, 1, 0, 0, 0], bool)
    anchors = np.array([
        [98., 30., 97.], [2., 30., 98.], [60., 30., 96.],
        [99., 31., 7.], [1., 30., 8.], [60., 30., 9.],
    ])
    protein = np.array([[99., 30., 50.], [1., 30., 50.]])
    def compute(points, group):
        result = enrichment_frame(species, upper, points,
                                  periodic_center(group, BOX), BOX, 0., 8., 2)
        return [(r["leaflet"], r["species"], r["count_local"], r["n_local"]) for r in result]
    expected = compute(anchors, protein)
    anchors2 = anchors.copy()
    anchors2[[0, 3], 0] -= 100
    protein2 = protein.copy()
    protein2[0, 0] -= 100
    assert compute(anchors2, protein2) == expected


def test_leaflet_tracking_invariant_to_z_wrap_and_label_reversal():
    first = np.array([94., 96., 98., 95., 4., 6., 8., 5.])
    reference = classify_leaflets(first, 100.)
    after = first.copy() + 102.
    after[[0, 3, 5]] -= 100.
    labels, changes = aligned_leaflets(after, 100., reference)
    np.testing.assert_array_equal(labels, reference)
    assert changes == 0


def test_local_normals_invariant_to_wrap():
    x, y = np.meshgrid(np.arange(4)*5.+2, np.arange(4)*5.+2)
    z = .12*x + .05*y + 20
    plane = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    points = np.vstack([plane, plane + [0, 0, 40.]])
    upper = np.array([False]*len(plane) + [True]*len(plane))
    angles, _ = fit_local_normals(points, upper, BOX, radius=12, min_neighbors=5)
    shifted = points.copy()
    shifted[::3, 0] += 100
    shifted[1::3, 1] -= 100
    shifted[2::3, 2] += 100
    angles2, _ = fit_local_normals(shifted, upper, BOX, radius=12, min_neighbors=5)
    np.testing.assert_allclose(angles2, angles, equal_nan=True, atol=1e-5)


def test_anchor_groups_vectorized_mapping():
    u = mda.Universe.empty(8, n_residues=8, atom_resindex=np.arange(8), trajectory=True)
    u.add_TopologyAttr("names", ["P"]*8)
    u.add_TopologyAttr("resnames", ["POPC"]*4+["POPE"]*4)
    u.add_TopologyAttr("resids", list(range(1, 9)))
    u.dimensions = BOX
    u.atoms.positions = np.array([
        [1, 1, 95], [2, 2, 96], [3, 3, 7], [4, 4, 8],
        [5, 5, 95], [6, 6, 96], [7, 7, 7], [8, 8, 8],
    ], float)
    species, resids, tracker, upper = anchor_groups(
        u, {"species_selections": {"POPC": "resname POPC and name P",
                                    "POPE": "resname POPE and name P"}})
    assert len(tracker) == 8
    assert len(np.unique(resids)) == 8
    np.testing.assert_allclose(tracker.positions(), u.atoms.positions)
    assert set(species) == {"POPC", "POPE"}
    assert upper.sum() == 4


def test_glycolipid_cross_residue_intramolecular_hbonds_excluded():
    # GM1a comprises two bonded residues, as does BGAL. The chainID
    # identifies the species but NOT the individual molecule.
    u = mda.Universe.empty(
        8, n_residues=4, atom_resindex=[0,0,1,1,2,2,3,3], trajectory=True
    )
    u.add_TopologyAttr("names", ["O","H","O","H"]*2)
    u.add_TopologyAttr("resnames", ["CER160", "BGAL", "CER160", "BGAL"])
    u.add_TopologyAttr("resids", [1,2,3,4])
    u.add_TopologyAttr("chainIDs", ["GM1a"]*4 + ["BGAL"]*4)
    u.add_TopologyAttr("bonds", [(0,1),(1,2),(2,3),(4,5),(5,6),(6,7)])
    ids, species = lipid_molecule_map(u, u.atoms)
    assert ids[0] == ids[3] and ids[4] == ids[7]
    assert ids[0] != ids[4]
    assert species[2] == "GM1a" and species[6] == "BGAL"
    raw = np.array([
        [0, 1, 0, 2, 2.7, 160.],  # same GM1a, different residues: exclude
        [0, 1, 0, 6, 2.7, 160.],  # GM1a to BGAL: keep
    ])
    events = tabulate_events(
        raw, u.atoms, [0], [100.], molecule_ids=ids, atom_species=species,
    )
    assert len(events) == 1
    assert events.iloc[0]["species_pair"] == "BGAL__GM1a"
    assert events.iloc[0]["donor_molecule_id"] != events.iloc[0]["acceptor_molecule_id"]


def test_lipid_hbond_lifetime_identity_includes_bonded_molecule_ids():
    keys = INTERACTION_KEY_COLUMNS["lipid_hydrogen_bonds"]
    assert "donor_molecule_id" in keys
    assert "acceptor_molecule_id" in keys
