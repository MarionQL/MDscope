"""End-to-end synthetic PDB/DCD exercising issue #3 and #4 analysis output contracts."""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from mdscope.analysis.config import Replicate
from mdscope.analysis.observables.lipid_enrichment import run_lipid_enrichment
from mdscope.analysis.observables.bilayer_normal import run_bilayer_normal
from mdscope.analysis.interactions.lipid_hydrogen_bonds import tabulate_events

def _tiny_lipid_universe(tmp_path):
    u=mda.Universe.empty(9,n_residues=9,atom_resindex=np.arange(9),
                         trajectory=True)
    u.add_TopologyAttr("names",["CA"]+["P"]*8)
    u.add_TopologyAttr("resnames",["ALA"]+["POPC"]*4+["POPE"]*4)
    u.add_TopologyAttr("resids",np.arange(1,10))
    u.add_TopologyAttr("masses",[12.]*9)
    xy=np.array([[10,10],[14,10],[10,14],[14,14]])
    points=np.concatenate([np.column_stack([xy,np.ones(4)*10]),
                           np.column_stack([xy,np.ones(4)*30])])
    coords=np.vstack([[12.,12.,20.],points])
    u.atoms.positions=coords
    u.dimensions=np.array([100,100,100,90,90,90],float)
    pdb=tmp_path/"membrane.pdb"
    with mda.Writer(str(pdb),n_atoms=9) as w:w.write(u)
    dcd=tmp_path/"membrane.dcd"
    with mda.Writer(str(dcd),n_atoms=9) as w:
        for k in range(5):
            u.atoms.positions=coords+np.array([.05*k,0,0])
            w.write(u)
    return pdb,dcd

def test_lipid_analyses_real_small_trajectory(tmp_path):
    pdb,dcd=_tiny_lipid_universe(tmp_path)
    selections={"POPC":"resname POPC and name P","POPE":"resname POPE and name P"}
    cfg={"selections":{"protein":"resname ALA and name CA"},
         "progress":{"enabled":False},
         "trajectory":{"start":None,"stop":None,"stride":1},
         "time":{"dt_ps":None},
         "analyses":{
             "lipid_enrichment":{"species_selections":selections,"axis":"z",
                "protein_selection":"resname ALA and name CA",
                "inner_radius_A":0.,"outer_radius_A":9.},
             "bilayer_normal":{"species_selections":selections,"axis":"z",
                "neighbor_radius_A":10.,"min_neighbors":3}}}
    rep=Replicate("test",1,pdb,dcd,"rep1")
    out=tmp_path/"test"/"rep1"
    a=run_lipid_enrichment(rep,cfg,out)
    b=run_bilayer_normal(rep,cfg,out)
    assert a["frames"]==5 and b["frames"]==5
    e=pd.read_parquet(out/"extensions/lipid_enrichment/timeseries.parquet")
    n=pd.read_parquet(out/"extensions/bilayer_normal/timeseries.parquet")
    assert len(e)==len(n)==5
    assert "lipid_enrichment_POPC_upper_fold" in e.columns
    assert np.all(np.isfinite(n["bilayer_normal_valid_fraction"]))
    assert (out/"extensions/bilayer_normal/timeseries.parquet").exists()

def test_lipid_hb_tabulation_excludes_same_residue():
    u=mda.Universe.empty(4,n_residues=2,atom_resindex=[0,0,1,1],trajectory=True)
    u.add_TopologyAttr("names",["O","H","O","H"])
    u.add_TopologyAttr("resnames",["POPC","POPE"])
    u.add_TopologyAttr("resids",[1,2])
    arr=np.array([[0,0,1,2,2.7,180],[0,0,1,0,2.0,170]],float)
    table=tabulate_events(arr,u.atoms,[0],[10.],exclude_same_residue=True)
    assert len(table)==1
    assert table.iloc[0]["species_pair"]=="POPC__POPE"
