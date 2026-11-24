"""Regression coverage for shared topology-aware sequence motif matching (#167)."""
from dataclasses import dataclass

import MDAnalysis as mda
import pytest

from mdscope.core.sequence import sequence_matches


@dataclass
class CA:
    resindex: int
    resid: int
    resname: str
    chainID: str = "A"
    segid: str = "PROA"
    segindex: int = 0
    name: str = "CA"


def test_overlapping_matches_have_stable_topology_indices():
    ca = [CA(i, 10 + i, aa) for i, aa in enumerate(("ALA", "ALA", "ALA"))]
    matches = sequence_matches(ca, {"di": "AA"})
    assert {k: v.tolist() for k, v in matches.items()} == {
        "di_A_1": [0, 1], "di_A_2": [1, 2],
    }
    nonoverlap = sequence_matches(ca, {"di": "AA"}, allow_overlaps=False)
    assert [v.tolist() for v in nonoverlap.values()] == [[0, 1]]


def test_identical_resids_on_different_chains_do_not_mix():
    ca = [CA(0, 10, "ALA"), CA(1, 11, "GLY"),
          CA(2, 10, "ALA", chainID="B", segid="PROB", segindex=1),
          CA(3, 11, "GLY", chainID="B", segid="PROB", segindex=1)]
    found = sequence_matches(ca, {"pair": "AG"})
    assert {k: v.tolist() for k, v in found.items()} == {
        "pair_A_1": [0, 1], "pair_B_1": [2, 3],
    }


def test_selection_gaps_cannot_create_artificial_motifs():
    ca = [CA(0, 10, "ALA"), CA(2, 12, "GLY")]
    assert sequence_matches(ca, {"pair": "AG"}) == {}


def test_repeated_chain_labels_across_segments_have_distinct_ids():
    ca = [CA(0, 1, "ALA"), CA(1, 2, "GLY"),
          CA(2, 1, "ALA", segid="SECOND", segindex=2),
          CA(3, 2, "GLY", segid="SECOND", segindex=2)]
    found = sequence_matches(ca, {"pair": "AG"})
    assert set(found) == {"pair_A_seg0_1", "pair_A_seg2_1"}


def test_modified_or_unknown_residue_does_not_match_wildcard():
    ca = [CA(0, 1, "ASN"), CA(1, 2, "MSE"), CA(2, 3, "PRO")]
    assert sequence_matches(ca, {"nxp": "N.P"}) == {}


def test_missing_motif_and_invalid_topology_selection():
    assert sequence_matches([CA(0, 1, "ALA")], {"pair": "AG"}) == {}
    with pytest.raises(ValueError, match="Duplicate CA"):
        sequence_matches([CA(0, 1, "ALA"), CA(0, 1, "ALA")], {"a": "A"})
    with pytest.raises(ValueError, match="CA-only"):
        sequence_matches([CA(0, 1, "ALA", name="CB")], {"a": "A"})
    with pytest.raises(ValueError, match="Invalid motif regex"):
        sequence_matches([CA(0, 1, "ALA")], {"broken": "["})


@pytest.mark.parametrize("chains,segments,expected", [
    (None, None, ["?:ALA1", "?:ALA1"]),
    (None, ["PROA", "PROB"], ["PROA:ALA1", "PROB:ALA1"]),
    (["", "  "], ["PROA", "PROB"], ["PROA:ALA1", "PROB:ALA1"]),
    (["A", "B"], ["same", "same"], ["A:ALA1", "B:ALA1"]),
    (["A", "B"], None, ["A:ALA1", "B:ALA1"]),
])
def test_atom_and_residue_labels_share_topology_identity(chains, segments, expected):
    from mdscope.core.sequence import residue_label
    u = mda.Universe.empty(2, n_residues=2, n_segments=2,
                           atom_resindex=[0, 1], residue_segindex=[0, 1])
    u.add_TopologyAttr("resnames", ["ALA", "ALA"])
    u.add_TopologyAttr("resids", [1, 1])
    if chains is not None:
        u.add_TopologyAttr("chainIDs", chains)
    if segments is not None:
        u.add_TopologyAttr("segids", segments)
    assert [residue_label(a) for a in u.atoms] == expected
    assert [residue_label(r) for r in u.residues] == expected
    assert list(u.residues.ix) == [0, 1]
