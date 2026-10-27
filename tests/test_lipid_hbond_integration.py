"""Bond-aware lipid-hydrogen-bond integration with explicit donor/acceptor selections."""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from mdscope.analysis.config import Replicate
from mdscope.analysis.lipid_hydrogen_bonds import run_lipid_hydrogen_bonds

def test_run_lipid_hb_on_bonded_pair(tmp_path):
    u=mda.Universe.empty(4,n_residues=2,atom_resindex=[0,0,1,1],trajectory=True)
    u.add_TopologyAttr("names",["O","H","O","H"])
    u.add_TopologyAttr("resnames",["POPC","POPE"])
    u.add_TopologyAttr("resids",[1,2])
    u.add_TopologyAttr("masses",[16.,1.,16.,1.])
    u.add_TopologyAttr("charges",[-.8,.4,-.8,.4])
    u.add_TopologyAttr("bonds",[(0,1),(2,3)])
    u.dimensions=np.array([50.,50.,50.,90.,90.,90.])
    pos=np.array([[5.,10.,10.],[6.,10.,10.],[7.7,10.,10.],[8.7,10.,10.]])
    u.atoms.positions=pos
    pdb=tmp_path/"lipids.pdb"
    with mda.Writer(str(pdb),n_atoms=4,bonds="all") as writer:
        writer.write(u)
    dcd=tmp_path/"lipids.dcd"
    with mda.Writer(str(dcd),n_atoms=4) as writer:
        for k in range(4):
            u.atoms.positions=pos + np.array([.01*k,0,0])
            writer.write(u)
    source=mda.Universe(str(pdb),str(dcd))
    assert len(source.bonds)==2
    cfg={"trajectory":{"start":None,"stop":None,"stride":1},
         "progress":{"enabled":False},"time":{"dt_ps":None},
         "analyses":{"lipid_hydrogen_bonds":{
             "lipid_selection":"resname POPC POPE","hydrogens_sel":"name H",
             "acceptors_sel":"name O", "donors_sel":None,
             "d_h_cutoff_A":1.2,"d_a_cutoff_A":3.0,"angle_cutoff_deg":150.,
             "exclude_same_residue":True,"update_selections":False}}}
    rep=Replicate("test",1,pdb,dcd,"rep1")
    out=tmp_path/"output"
    result=run_lipid_hydrogen_bonds(rep,cfg,out)
    events=pd.read_parquet(out/"interactions/lipid_hydrogen_bonds/events.parquet")
    assert result["frames"]==4
    assert len(events)==4
    assert set(events["species_pair"])=={"POPC__POPE"}
    assert (out/"interactions/lipid_hydrogen_bonds/count_timeseries.parquet").exists()
