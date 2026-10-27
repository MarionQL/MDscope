"""Per-frame protein-lipid/HS contacts and motif proximity.\n\n    Leaflet-normalized enrichment is owned by lipid_enrichment.py, not duplicated here.

All scalar observables are written into MDscope's normal summary time series,
so they inherit identical hierarchical bootstrap and plotting code.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import distance_array
from .biophysical_common import heavy,box_of,residue_pairs,sequence_motifs
from .config import atomic_parquet,ensure_dir,time_ps

def run_local_environment(rep,cfg,repout):
    import MDAnalysis as mda
    settings=cfg["analyses"].get("local_environment",{})
    u=mda.Universe(str(rep.topology),str(rep.trajectory))
    prot=heavy(u.select_atoms(settings.get("protein_selection") or cfg["selections"]["protein"]))
    lip=heavy(u.select_atoms(settings.get("lipid_selection","resname POPC POPE POPI POPS DOPC DPPC PSM CHL1 CER160 CER1 BGAL GM1a")))
    if not len(lip):raise ValueError("No lipids matched local_environment.lipid_selection")
    if not len(prot):raise ValueError("No protein atoms matched local_environment.protein_selection")
    head=heavy(u.select_atoms(settings.get("headgroup_selection","resname POPC POPE POPI POPS DOPC DPPC PSM and name P")))
    glycan=heavy(u.select_atoms(settings.get("glycan_selection","resname AGLCNA AIDOA BGLCA BGAL BXYL")))
    ions=heavy(u.select_atoms(settings.get("ion_selection","resname SOD NA CLA CL POT K CAL CA")))
    anchor=heavy(u.select_atoms(settings.get("anchor_selection","protein and name CA and (resid 151-170 or resid 349-368)")))
    if not len(anchor):anchor=prot
    cutoff=float(settings.get("contact_cutoff_A",5.0))
    motifs=sequence_motifs(prot.residues,settings.get("motifs",{"nxip":"N.IP","gxxxg":"G...G","hs_candidate":"SG"}))
    motif_groups={}
    for hit in motifs:
        key="motif_"+hit["motif"]+"_"+hit["chain"]+"_"+str(hit["start_resid"])
        motif_groups[key]=heavy(u.atoms[np.isin(u.atoms.resindices,hit["resindices"])])
    if not motifs:
        print("WARNING: no configured motifs matched sequence")
    grouped={}
    for name in np.unique(lip.resnames):
        grouped[str(name)]=lip[lip.resnames==name]
    # Use identical lipid species denominator across time, with the same selection.
    rows=[]; per_res=[]; per_motif=[]; meta=[]
    for key,grp in motif_groups.items():
        meta.append({"motif_id":key,"sequence":"".join([]),"resids":",".join(map(str,grp.residues.resids))})
    tr=cfg["trajectory"]
    for frame_i,ts in enumerate(u.trajectory[slice(tr.get("start"),tr.get("stop"),int(tr.get("stride",1)))]):
        t=time_ps(ts,frame_i,cfg); box=box_of(u)
        row={"time_ps":t}
        cp=residue_pairs(prot,lip,cutoff,box)
        contacted=np.unique(cp[:,0]) if len(cp) else np.empty(0,int)
        row["local_protein_lipid_contact_residues_count"]=len(contacted)
        if len(cp):
            a,b=np.unique(cp[:,0],return_counts=True)
            for resid,n in zip(a,b):
                per_res.append((t,int(resid),int(n)))
        for species,atoms in grouped.items():
            hits=residue_pairs(prot,atoms,cutoff,box)
            row["local_"+species+"_contact_residues_count"]=len(np.unique(hits[:,0])) if len(hits) else 0
        for name,group in motif_groups.items():
            if len(group):
                d=distance_array(group.positions,lip.positions,box=box)
                row[name+"_min_lipid_distance_A"]=float(d.min())
                row[name+"_lipid_contact_fraction"]=float(np.mean(d.min(axis=1)<=cutoff))
                # Geometric occlusion: not equivalent to solvent-accessible area.
                others=prot[~np.isin(prot.indices,group.indices)]
                if len(others):
                    d2=distance_array(group.positions,others.positions,box=box)
                    row[name+"_protein_occlusion_fraction"]=float(np.mean(np.any(d2<cutoff,axis=1)))
        if len(glycan):
            contacts=residue_pairs(glycan,lip,cutoff,box)
            row["local_hs_lipid_contact_pairs_count"]=len(contacts)
            if len(ions):
                row["local_hs_ion_contacts_count"]=len(residue_pairs(glycan,ions,3.5,box))
        rows.append(row)
    out=ensure_dir(repout/"extensions"/"local_environment")
    atomic_parquet(out/"timeseries.parquet",pd.DataFrame(rows))
    if per_res:atomic_parquet(out/"residue_contacts.parquet",pd.DataFrame(per_res,columns=["time_ps","resindex","n_lipid_residues"]))
    if meta:atomic_parquet(out/"motifs.parquet",pd.DataFrame(meta))
    return {"frames":len(rows),"metrics":len(pd.DataFrame(rows).columns)-1}
