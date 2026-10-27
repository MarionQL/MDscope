"""Intermolecular lipid–lipid H bonds by lipid-species pair.

Uses one MDAnalysis HBA pass for all selected lipid atoms. A pair is counted
only when donor and acceptor belong to different residues; lipid/water and
lipid/protein hydrogen bonds are excluded by selection.
"""
from __future__ import annotations
from pathlib import Path
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
from .config import ensure_dir,atomic_parquet,time_ps
from .interaction_tables import atom_label,chain_label,count_timeseries

def ordered_pair(a,b):
    """Canonical unordered species pair, so donor/acceptor direction doesn't split totals."""
    return "__".join(sorted((str(a),str(b))))

def tabulate_events(event_array,atoms,frames,times,exclude_same_residue=True):
    cols=["time_ps","donor","hydrogen","acceptor","donor_species","acceptor_species",
          "species_pair","scope","interaction_class","distance_A","angle_deg"]
    if np.asarray(event_array).size==0:return pd.DataFrame(columns=cols)
    data=np.asarray(event_array)
    donor=data[:,1].astype(int)
    hydrogen=data[:,2].astype(int)
    acceptor=data[:,3].astype(int)
    kept=np.ones(len(data),dtype=bool)
    if exclude_same_residue:
        kept=atoms[donor].resindices!=atoms[acceptor].resindices
    data=data[kept]
    donor=donor[kept];hydrogen=hydrogen[kept];acceptor=acceptor[kept]
    frame_map=dict(zip(np.asarray(frames,int),np.asarray(times,float)))
    sp1=np.asarray(atoms[donor].resnames,str)
    sp2=np.asarray(atoms[acceptor].resnames,str)
    classes=[ordered_pair(a,b) for a,b in zip(sp1,sp2)]
    return pd.DataFrame({
        "time_ps":[frame_map[int(f)] for f in data[:,0]],
        "donor":[atom_label(a) for a in atoms[donor]],
        "hydrogen":[atom_label(a) for a in atoms[hydrogen]],
        "acceptor":[atom_label(a) for a in atoms[acceptor]],
        "donor_species":sp1,"acceptor_species":sp2,"species_pair":classes,
        "scope":["intraspecies" if a==b else "interspecies" for a,b in zip(sp1,sp2)],
        "interaction_class":classes,
        "distance_A":data[:,4].astype(float),
        "angle_deg":data[:,5].astype(float),
    },columns=cols)

def run_lipid_hydrogen_bonds(rep,cfg,repout):
    spec=cfg["analyses"]["lipid_hydrogen_bonds"]
    universe=mda.Universe(str(rep.topology),str(rep.trajectory))
    selection=str(spec.get("lipid_selection","resname POPC POPE POPI DOPC DPPC PSM CHL1 BGAL GM1a CER160 CER1"))
    atoms=universe.select_atoms(selection)
    if not len(atoms):raise ValueError("lipid_hydrogen_bonds.lipid_selection matched no atoms")
    hba=HydrogenBondAnalysis(
        universe,
        donors_sel=spec.get("donors_sel"),
        hydrogens_sel=spec.get("hydrogens_sel"),
        acceptors_sel=spec.get("acceptors_sel"),
        d_h_cutoff=float(spec.get("d_h_cutoff_A",1.2)),
        d_a_cutoff=float(spec.get("d_a_cutoff_A",3.0)),
        d_h_a_angle_cutoff=float(spec.get("angle_cutoff_deg",150.0)),
        update_selections=bool(spec.get("update_selections",False)),
    )
    # Constrain every donor H and acceptor to selected lipids; with bonded
    # topology donors follow the selected hydrogen. User selections override.
    if spec.get("hydrogens_sel") is None:
        hba.hydrogens_sel=hba.guess_hydrogens(selection)
    if spec.get("acceptors_sel") is None:
        hba.acceptors_sel=hba.guess_acceptors(selection)
    if not hba.hydrogens_sel or not hba.acceptors_sel:
        raise ValueError("Could not identify lipid donor hydrogens/acceptors; provide explicit selections")
    tr=cfg["trajectory"]
    try:
        hba.run(start=tr.get("start"),stop=tr.get("stop"),step=int(tr.get("stride",1)))
    except Exception as exc:
        raise RuntimeError("Lipid H-bond analysis requires explicit H atoms/bond connectivity; verify selections") from exc
    times=np.asarray(hba.times,float)
    if cfg.get("time",{}).get("dt_ps") is not None:
        times=cfg["time"].get("origin_ps",0.)+np.arange(len(times))*float(cfg["time"]["dt_ps"])
    df=tabulate_events(hba.results.hbonds,universe.atoms,hba.frames,times,
                      exclude_same_residue=bool(spec.get("exclude_same_residue",True)))
    out=ensure_dir(Path(repout)/"interactions"/"lipid_hydrogen_bonds")
    atomic_parquet(out/"events.parquet",df)
    counts=count_timeseries(df,times,"lipid_hydrogen_bond","interaction_class",
                             ["donor","hydrogen","acceptor","interaction_class"])
    # Preserve zero-count rows even when no bonds occur at any frame.
    atomic_parquet(out/"count_timeseries.parquet",counts)
    return {"frames":len(times),"events":len(df),"species_pairs":df["species_pair"].nunique()}
