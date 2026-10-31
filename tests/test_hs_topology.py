"""Bond-verified HS attachment discovery, including shared protein chain IDs."""
import MDAnalysis as mda
import numpy as np
import pytest

from mdscope.analysis.mechanisms.hs_topology import discover_covalent_hs


def _universe():
    names = ["OG", "OG"]
    resnames = ["SER", "SER"]
    resids = [39, 309]
    chains = ["PROA", "PROB"]
    owners = [0, 1]
    bonds = []
    for chain, parent, resid0 in (("PROA", 0, 1), ("PROB", 1, 100)):
        previous = parent
        for j, name in enumerate(("BXYL", "AGLCNA", "AIDOA")):
            atom = len(names)
            names.append("C1")
            resnames.append(name)
            resids.append(resid0 + j)
            chains.append(chain)  # Same chainID as protein; cannot group by chainID.
            owners.append(len(resnames)-1)
            bonds.append((previous, atom))
            previous = atom
    u = mda.Universe.empty(
        len(names), n_residues=len(resnames),
        atom_resindex=owners, trajectory=True,
    )
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", resids)
    u.add_TopologyAttr("chainIDs", chains)
    u.add_TopologyAttr("bonds", bonds)
    return u


def test_bond_verified_hs_auto_mapping():
    u = _universe()
    hs, sites, audit = discover_covalent_hs(u, {
        "glycan_selection": "resname BXYL AGLCNA AIDOA",
        "attachment_selection": "resname SER and name OG",
        "expected_chains": 2,
        "min_glycan_residues": 3,
    })
    assert set(hs) == {"PROA_S39", "PROB_S309"}
    assert all(len(g.residues) == 3 for g in hs.values())
    assert sites["PROA_S39"].resids.tolist() == [39]
    assert audit["PROB_S309"]["protein_atom"] == "OG"


def test_hs_auto_rejects_incorrect_chain_count():
    with pytest.raises(ValueError, match="Expected 6"):
        discover_covalent_hs(_universe(), {
            "glycan_selection": "resname BXYL AGLCNA AIDOA",
            "attachment_selection": "resname SER and name OG",
            "expected_chains": 6,
        })


def test_hs_auto_rejects_missing_covalent_attachment():
    with pytest.raises(ValueError, match="selections matched no glycan or attachment atoms"):
        discover_covalent_hs(_universe(), {
            "glycan_selection": "resname BXYL AGLCNA AIDOA",
            "attachment_selection": "resname SER and name CB",
        })
