"""local_environment motifs use the shared matcher and keep legacy identities (#167)."""
import MDAnalysis as mda
import numpy as np
import pytest

from mdscope.config import resolve_config
from mdscope.mechanisms.local_environment import motif_groups

# Segments SA, SB, SC (chainIDs A, B, C), resids restart at 1 in each.
SEGMENTS = {"SA": ["ASN", "HIS", "ILE", "PRO", "GLY"],
            "SB": ["ASN", "ALA", "ILE", "PRO"],
            "SC": ["ASN", "HIS", "ALA", "ILE", "PRO"]}


def _universe():
    resnames = [name for names in SEGMENTS.values() for name in names]
    segindex = [i for i, names in enumerate(SEGMENTS.values()) for _ in names]
    u = mda.Universe.empty(2 * len(resnames), n_residues=len(resnames), n_segments=len(SEGMENTS),
                           atom_resindex=np.repeat(np.arange(len(resnames)), 2),
                           residue_segindex=segindex, trajectory=True)
    u.add_TopologyAttr("names", ["CA", "CB"] * len(resnames))
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", [i + 1 for names in SEGMENTS.values() for i in range(len(names))])
    u.add_TopologyAttr("segids", list(SEGMENTS))
    u.add_TopologyAttr("chainIDs", [segid[1] for segid, names in SEGMENTS.items()
                                    for _ in names for _ in range(2)])
    return u


def test_regex_motifs_keep_segid_keys_and_chain_then_motif_order():
    groups = motif_groups(_universe(), {"motifs": {"pro": "P", "nxip": "N.IP"}}, {},
                          "protein and segid SA SB")
    assert list(groups) == ["motif_pro_SA_4", "motif_nxip_SA_1", "motif_pro_SB_4", "motif_nxip_SB_1"]
    assert groups["motif_nxip_SA_1"].residues.resids.tolist() == [1, 2, 3, 4]
    assert len(groups["motif_nxip_SA_1"]) == 8


def test_selection_gap_never_creates_a_motif():
    # Without SC resid 3 the selected residues read NHIP, but they are not contiguous.
    groups = motif_groups(_universe(), {"motifs": {"nxip": "N.IP"}}, {},
                          "protein and not (segid SC and resid 3)")
    assert list(groups) == ["motif_nxip_SA_1", "motif_nxip_SB_1"]


def test_reference_uses_the_shared_definition_scope():
    shared = {"patch": {"pattern": "N.IP", "chains": ["B"]}}
    groups = motif_groups(_universe(), {"motifs": {"nxip": {"motif": "patch"}}}, shared, "protein")
    assert list(groups) == ["motif_nxip_SB_1"]


def test_default_motifs_apply_when_unset():
    assert list(motif_groups(_universe(), {}, {}, "protein")) == ["motif_nxip_SA_1", "motif_nxip_SB_1"]


def _document(motifs):
    from mdscope.studio.model import new_document
    doc = new_document()
    doc["selections"]["motifs"] = {"patch": {"sequence": "NHIP"}}
    doc["analyses"]["local_environment"] = {"enabled": False, "motifs": motifs}
    return doc


def test_regex_and_reference_entries_resolve():
    resolved = resolve_config(_document({"a": "N.IP", "b": {"motif": "patch"}}), check_paths=False)
    assert resolved["analyses"]["local_environment"]["motifs"]["b"] == {"motif": "patch"}


@pytest.mark.parametrize("motifs", [
    "N.IP", {"x": "["}, {"x": ".*"}, {"x": {"motif": "unknown"}},
    {"x": {"motif": "patch", "chains": ["A"]}},
])
def test_invalid_local_environment_motifs_fail_at_yaml_boundary(motifs):
    with pytest.raises(ValueError, match="motif"):
        resolve_config(_document(motifs), check_paths=False)
