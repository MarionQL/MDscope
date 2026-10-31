"""Regression tests for MDscope repository issues #1–#5."""
import numpy as np
from mdscope.core.lipid_geometry import classify_leaflets
from mdscope.observables.lipid_enrichment import shell_fold, enrichment_frame
from mdscope.observables.bilayer_normal import fit_local_normals
from mdscope.interactions.lipid_hydrogen_bonds import ordered_pair
from mdscope.observables.membrane_order import OrderSite,_profiles,_metric_name
from mdscope.pipeline.progress import track

def test_lipid_order_leaflets_and_position_profiles():
    sites=[
        OrderSite(atom_index=k,residue_index=k//4,resname="POPC",tail="sn1",
                  position=k%4+1,atom_name="C"+str(k))
        for k in range(16)
    ]
    leaf={0:"upper",1:"upper",2:"lower",3:"lower"}
    p=_profiles(sites,{"per_leaflet":True,"per_lipid_type":True},leaf)
    names={x.profile_id for x in p}
    assert "lipid_order__sn1" in names
    assert "lipid_order__sn1__upper" in names
    assert "lipid_order__sn1__lower" in names
    assert "lipid_order__POPC__sn1__upper" in names
    assert _metric_name("lipid_order__POPC__sn1__upper")=="lipid_order_POPC_sn1_upper_mean"

def test_circular_leaflet_assignment():
    # Bilayer spans the z boundary: each leaflet remains a separate cluster.
    z=np.array([96,97,95,98,7,8,6,9],float)
    upper=classify_leaflets(z,100.)
    assert len(np.unique(upper[:4]))==1
    assert len(np.unique(upper[4:]))==1
    assert upper[0]!=upper[4]

def test_enrichment_has_correct_bulk_denominator():
    species=np.array(["POPC","POPC","POPC","POPE",
                      "POPC","POPE","POPE","POPE"])
    upper=np.array([1,1,1,1,0,0,0,0],bool)
    positions=np.array([[1,0,20],[2,0,20],[20,0,20],[20,10,20],
                        [1,0,40],[2,0,40],[20,0,40],[20,10,40]],float)
    result=enrichment_frame(species,upper,positions,np.array([0,0,30]),
                            np.array([100,100,100,90,90,90]),0,5,2)
    u={(row["species"],row["leaflet"]):row for row in result}
    assert np.isclose(u["POPC","upper"]["enrichment_fold"],4/3)
    assert np.isclose(u["POPE","upper"]["enrichment_fold"],0)
    assert np.isclose(u["POPC","lower"]["enrichment_fold"],2)
    assert np.isnan(shell_fold(0,0,1,10))

def test_local_plane_angle():
    x,y=np.meshgrid(np.arange(4)*4+10,np.arange(4)*4+10)
    plane=np.column_stack([x.ravel(),y.ravel(),(.1*x+12).ravel()])
    upper=plane.copy(); upper[:,2]+=35
    xyz=np.concatenate([plane,upper])
    mask=np.array([0]*len(plane)+[1]*len(upper),bool)
    angles,neighbors=fit_local_normals(xyz,mask,np.array([100,100,100,90,90,90]),11,5,2)
    good=angles[np.isfinite(angles)]
    assert len(good)>18
    np.testing.assert_allclose(good,np.degrees(np.arctan(.1)),atol=.15)

def test_hbond_species_pair_is_undirected():
    assert ordered_pair("POPC","POPE")==ordered_pair("POPE","POPC")
    assert ordered_pair("POPC","POPC")=="POPC__POPC"

def test_progress_disabled():
    assert list(track(range(9),total=9,enabled=False))==list(range(9))
