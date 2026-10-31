"""Multi-residue glycolipid and motif-boundary regressions for issue #31."""
from types import SimpleNamespace

import MDAnalysis as mda
import numpy as np
import pytest

from mdscope.core.lipid_geometry import anchor_groups
from mdscope.biophysical_common import sequence_motifs
from mdscope.mechanisms.mechanistic_common import sequence_matches


def _glycolipid_universe():
    # Four molecules, each with a ceramide anchor residue and a carbohydrate
    # residue.  Carbohydrates must NOT be counted as additional lipids.
    u = mda.Universe.empty(
        12, n_residues=8,
        atom_resindex=[0, 0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7],
        trajectory=True,
    )
    u.add_TopologyAttr("resnames", ["CER", "SUG", "CER", "SUG",
                                    "CER", "SUG", "CER", "SUG"])
    u.add_TopologyAttr("resids", [100, 101, 200, 201, 300, 301, 400, 401])
    u.add_TopologyAttr("names", ["C1", "C2", "O1"] * 4)
    u.atoms.positions = np.array([
        [29.5, 1, 5], [0.5, 1, 5], [4, 3, 5],
        [5, 5, 6], [7, 5, 6], [8, 8, 6],
        [10, 10, 25], [12, 10, 25], [12, 12, 25],
        [20, 20, 26], [22, 20, 26], [23, 20, 26],
    ], float)
    u.dimensions = [30, 30, 30, 90, 90, 90]
    return u


def test_multiresidue_glycolipids_count_once_and_cross_pbc():
    u = _glycolipid_universe()
    names, indices, tracker, leaf = anchor_groups(
        u, {"species_selections": {"glycolipid": "resname CER"},
            "axis": "z"},
    )
    assert names.tolist() == ["glycolipid"] * 4
    np.testing.assert_array_equal(indices, [0, 2, 4, 6])
    assert len(tracker) == 4
    assert len(leaf) == 4 and leaf.sum() == 2
    positions = tracker.positions()
    # First ceramide C1/C2 straddles x=0. Its center should be at
    # the periodic boundary, not the misleading wrapped midpoint x=15.
    assert np.isclose(positions[0, 0] % 30, 0, atol=1e-5)
    wrapped_anchors = u.select_atoms("resname CER")
    np.testing.assert_array_equal(np.unique(wrapped_anchors.resindices), [0, 2, 4, 6])
    # A naive, nonperiodic mean of the first two anchor atoms is x=15.
    assert np.isclose(wrapped_anchors.positions[:2].mean(axis=0)[0], 15)
    # The different values are intentional: extended centers use wrapped COG.


def test_multiresidue_selection_overlap_is_rejected():
    u = _glycolipid_universe()
    with pytest.raises(ValueError, match="overlap"):
        anchor_groups(
            u, {"species_selections": {"first": "resname CER",
                                        "second": "resid 100"}},
        )


def test_motif_search_does_not_bridge_two_chains():
    # A G at the end of chain A and a G at the beginning of chain B
    # must not be detected as an adjacent GG motif.
    residues = [
        SimpleNamespace(chainID="A", segid="A", resname="ALA", resid=10, ix=0),
        SimpleNamespace(chainID="A", segid="A", resname="GLY", resid=20, ix=1),
        SimpleNamespace(chainID="B", segid="B", resname="GLY", resid=90, ix=2),
        SimpleNamespace(chainID="B", segid="B", resname="ALA", resid=100, ix=3),
    ]
    assert sequence_motifs(residues, {"glycine_pair": "GG"}) == []
    # A positive match within B retains actual noncontiguous topology IDs.
    residues.append(SimpleNamespace(chainID="B", segid="B",
                                    resname="GLY", resid=140, ix=4))
    result = sequence_motifs(residues, {"gag": "GAG"})
    assert len(result) == 1
    assert result[0]["chain"] == "B"
    assert result[0]["resindices"] == [2, 3, 4]


def test_mechanistic_overlapping_motif_matches_are_intentionally_retained():
    atoms = [
        SimpleNamespace(chainID="A", segid="A", resname="GLY", resindex=i)
        for i in range(3)
    ]
    atoms += [SimpleNamespace(chainID="B", segid="B", resname="GLY",
                              resindex=3)]
    hits = sequence_matches(atoms, {"pair": "GG"})
    assert set(hits) == {"pair_A_1", "pair_A_2"}
    np.testing.assert_array_equal(hits["pair_A_1"], [0, 1])
    np.testing.assert_array_equal(hits["pair_A_2"], [1, 2])
