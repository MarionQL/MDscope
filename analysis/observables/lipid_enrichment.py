"""Leaflet-normalized lipid enrichment/depletion around a protein selection.

Raw counts and replicate-wise scalar series feed existing MDscope statistics.
A lipid is counted once using a single representative anchor per residue.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from ..config import atomic_parquet,ensure_dir,time_ps
from ..lipid_geometry import anchor_groups,centers,lateral_distances,periodic_center,aligned_leaflets
from ..progress import progress

def shell_fold(observed,local_total,bulk_count,bulk_total):
    """Local mole fraction divided by bulk mole fraction in SAME leaflet."""
    if min(local_total,bulk_total)<=0 or bulk_count<=0:
        return np.nan
    return float((observed/local_total)/(bulk_count/bulk_total))

def enrichment_frame(species,upper,positions,protein_center,box,inner,outer,normal):
    """Return long-form per-leaflet species metrics and raw denominators."""
    distance=lateral_distances(positions,protein_center,box,normal)
    local=(distance>=inner)&(distance<outer)
    categories=np.unique(species)
    rows=[]
    for leaf,mask in (("upper",upper),("lower",~upper)):
        bulk_total=int(mask.sum())
        near=mask&local
        n_local=int(near.sum())
        for name in categories:
            smask=(species==name)
            count=int((near&smask).sum())
            bulk=int((mask&smask).sum())
            fold=shell_fold(count,n_local,bulk,bulk_total)
            rows.append({"leaflet":leaf,"species":name,
                         "count_local":count,"count_bulk":bulk,
                         "n_local":n_local,"n_bulk":bulk_total,
                         "local_fraction":count/n_local if n_local else np.nan,
                         "bulk_fraction":bulk/bulk_total if bulk_total else np.nan,
                         "enrichment_fold":fold,
                         "depleted":bool(fold<1.) if np.isfinite(fold) else False})
    return rows

def run_lipid_enrichment(rep,cfg,repout):
    spec=cfg["analyses"]["lipid_enrichment"]
    u=mda.Universe(str(rep.topology),str(rep.trajectory))
    tr=cfg["trajectory"]
    start=tr.get("start");stop=tr.get("stop");stride=int(tr.get("stride",1))
    u.trajectory[0 if start is None else int(start)]
    species,indices,groups,leaf=anchor_groups(u,spec)
    protein=u.select_atoms(str(spec.get("protein_selection") or cfg["selections"]["protein"]))
    if not len(protein):raise ValueError("lipid_enrichment.protein_selection is empty")
    outer=float(spec.get("outer_radius_A",15.0))
    inner=float(spec.get("inner_radius_A",0.0))
    if not 0<=inner<outer:raise ValueError("Require 0 <= inner_radius_A < outer_radius_A")
    normal={"x":0,"y":1,"z":2}[str(spec.get("axis","z")).lower()]
    rows=[];tables=[]
    frames=u.trajectory[slice(start,stop,stride)]
    with progress(frames,total=len(frames),desc=f"{rep.system}/{rep.label}: lipid enrichment",
                  enabled=cfg.get("progress",{}).get("enabled",True)) as iterator:
        for k,ts in enumerate(iterator):
            t=time_ps(ts,k,cfg)
            box=np.asarray(u.dimensions,float)
            anchors = centers(groups)
            current_leaf, changes = aligned_leaflets(anchors[:, normal], float(box[normal]), leaf)
            # A Cartesian center_of_geometry is incorrect for a protein dimer
            # straddling a periodic boundary. Use a circular lattice center.
            center = periodic_center(protein.positions, box)
            vals=enrichment_frame(species,current_leaf,anchors,center,box,inner,outer,normal)
            row={"time_ps":t, "lipid_enrichment_leaflet_changes": changes}
            for value in vals:
                name=f"lipid_enrichment_{value['species']}_{value['leaflet']}"
                row[name+"_fold"]=value["enrichment_fold"]
                row[name+"_local_fraction"]=value["local_fraction"]
                row[name+"_contact_count"]=value["count_local"]
                tables.append({"time_ps":t,**value})
            rows.append(row)
    out=ensure_dir(repout/"extensions"/"lipid_enrichment")
    atomic_parquet(out/"timeseries.parquet",pd.DataFrame(rows))
    atomic_parquet(out/"counts.parquet",pd.DataFrame(tables))
    atomic_parquet(out/"selection_audit.parquet",pd.DataFrame({
        "species":species,"resindex":indices,"leaflet":np.where(leaf,"upper","lower")}))
    return {"frames":len(rows),"n_lipids":len(groups),"n_species":len(np.unique(species))}
