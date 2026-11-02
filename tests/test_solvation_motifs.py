"""Sequence-defined solvation targets must resolve topology identities safely."""
import numpy as np
import pytest
import MDAnalysis as mda

from mdscope.observables.solvation import _targets


def _protein():
    u=mda.Universe.empty(8,n_residues=4,
                         atom_resindex=[0,0,1,1,2,2,3,3],
                         trajectory=True)
    u.add_TopologyAttr("names",["CA","CB"]*4)
    u.add_TopologyAttr("resnames",["ASP","GLU","ASP","GLU"])
    u.add_TopologyAttr("resids",[1,2,3,4])
    u.add_TopologyAttr("chainIDs",["PROA"]*4+["PROB"]*4)
    return u


def test_motif_targets_are_chain_specific_whole_residues():
    u=_protein()
    result=_targets(u,{"targets":{"motif":{"acidic_patch":{
        "sequence":"DE","selection":"protein and name CA"}}}})
    assert len(result)==2
    assert {row[1:3] for row in result}=={("motif","acidic_patch")}
    assert all(len(row[3])==4 for row in result)
    assert len(set(row[0] for row in result))==2


def test_motif_missing_sequence_fails():
    u=_protein()
    with pytest.raises(ValueError,match="not found"):
        _targets(u,{"targets":{"motif":{"acidic_patch":{"sequence":"DDEDVVG"}}}})


def test_motif_invalid_regex_or_overlap_fails():
    u=_protein()
    with pytest.raises(ValueError,match="exact amino-acid sequence"):
        _targets(u,{"targets":{"motif":{"bad":{"sequence":"D.*E"}}}})
    with pytest.raises(ValueError,match="overlaps"):
        _targets(u,{"targets":{"motif":{"first":{"sequence":"DE"},
                                          "second":{"sequence":"DE"}}}})



def test_sequence_mapping_supports_explicit_chemical_site_subset():
    u=_protein()
    sites=_targets(u,{'targets':{'motif':{'acidic':{'sequence':'DE','atom_selection':'name CB'}}}})
    assert len(sites)==2
    assert all(len(group)==2 and set(group.names)=={'CB'} for _,_,_,group in sites)
    with pytest.raises(ValueError,match='no heavy atoms'):
        _targets(u,{'targets':{'motif':{'acidic':{'sequence':'DE','atom_selection':'name O*'}}}})
