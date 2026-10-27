"""Exploratory, predeclared biochemically motivated within-replicate coupling screen.

Does not return p-values. Rates correlations by reproducibility across replicates,
and refuses to pool frames across replicates as independent observations.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from .config import ensure_dir,atomic_parquet

PAIRS=[
 ("lipid_proximity","local_hs_lipid_contact_pairs_count","local_protein_lipid_contact_residues_count"),
 ("ion_screening","local_hs_ion_contacts_count","local_hs_lipid_contact_pairs_count"),
 ("lipid_contact_rigidity","local_protein_lipid_contact_residues_count","radius_of_gyration_A"),
 ("membrane_coupling","local_protein_lipid_contact_residues_count","membrane_thickness_A")
]

def block_means(x,block):
    n=len(x)//block
    return x[:n*block].reshape(n,block).mean(axis=1) if n else np.empty(0)

def correlation(x,y,min_blocks=8):
    a=np.asarray(x,float);b=np.asarray(y,float)
    good=np.isfinite(a)&np.isfinite(b)
    a,b=a[good],b[good]
    if len(a)<min_blocks or np.std(a)==0 or np.std(b)==0:return np.nan
    return float(np.corrcoef(a,b)[0,1])

def screen_couplings(cfg,root,system_repouts):
    spec=cfg.get("analyses",{}).get("coupling_screen",{})
    pairs=spec.get("pairs",[])
    pairs=([{"name":n,"x":x,"y":y} for n,x,y in PAIRS]+pairs) if spec.get("include_defaults",True) else pairs
    window=int(spec.get("block_frames",50))
    if window<1:raise ValueError("block_frames must be positive")
    rows=[]
    for system,reps in system_repouts.items():
        for repout in reps:
            source=repout/"summary"/"timeseries.parquet"
            if not source.exists():continue
            df=pd.read_parquet(source)
            for pair in pairs:
                name,x,y=str(pair["name"]),str(pair["x"]),str(pair["y"])
                if x not in df or y not in df:continue
                clean=df[[x,y]].replace([np.inf,-np.inf],np.nan).dropna()
                if len(clean)<window*8:continue
                xs=block_means(clean[x].to_numpy(float),window)
                ys=block_means(clean[y].to_numpy(float),window)
                r=correlation(xs,ys)
                rows.append({"system":system,"replicate":repout.name,"mechanism":name,
                    "x":x,"y":y,"pearson_block_r":r,"block_frames":window,
                    "n_blocks":len(xs),"interpretation":"exploratory_association_not_causation"})
    if not rows:return
    out=ensure_dir(root/"project_summary"/"mechanisms")
    frame=pd.DataFrame(rows)
    atomic_parquet(out/"replicate_correlations.parquet",frame)
    summaries=[]
    for keys,d in frame.groupby(["system","mechanism","x","y"],dropna=False):
        rs=d["pearson_block_r"].to_numpy(float);rs=rs[np.isfinite(rs)]
        if not len(rs):continue
        summaries.append(dict(zip(("system","mechanism","x","y"),keys),
            n_replicates=len(rs),median_r=float(np.median(rs)),
            min_r=float(rs.min()),max_r=float(rs.max()),
            fraction_same_sign=float(max(np.mean(rs>0),np.mean(rs<0))),
            note="associational; no causal or statistical significance claim"))
    if summaries:atomic_parquet(out/"mechanism_consistency.parquet",pd.DataFrame(summaries))
    return out
