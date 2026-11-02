"""Concise manuscript-facing plots constructed ONLY from the standard statistics tables.

Do not independently fit uncertainty or pool correlated trajectory frames.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..config import ensure_dir, presentation_enabled
from .figure_inventory import current_statistics

def plot_lipid_enrichment_matrix(cfg,root):
    if not current_statistics(cfg, root, 'time_series'):
        return
    source=Path(root)/"project_summary"/"statistics"/"time_series"/"replicate_statistics.parquet"
    if not source.exists() or not presentation_enabled(cfg,"plots","lipid_enrichment"):return
    frame=pd.read_parquet(source)
    if not {"system","metric","mean"}.issubset(frame):return
    selected=frame[frame["metric"].str.startswith("lipid_enrichment_")
                   &frame["metric"].str.endswith("_fold")].copy()
    if selected.empty:return
    selected["species"]=selected["metric"].str.replace("^lipid_enrichment_","",regex=True).str.replace("_(upper|lower)_fold$","",regex=True)
    selected["leaflet"]=selected["metric"].str.extract("_(upper|lower)_fold$")[0]
    # Replicate means are averaged equally; not a pooled per-frame average.
    selected["column"]=selected["system"].astype(str)+" | "+selected["leaflet"].astype(str)
    grouped=selected.groupby(["species","column"],dropna=False)["mean"].mean().reset_index()
    table=grouped.pivot(index="species",columns="column",values="mean")
    if table.empty:return
    values=table.to_numpy(float)
    signed=np.full_like(values,np.nan,dtype=float)
    good=np.isfinite(values)&(values>=0)
    # Plot 0 observed enrichment on a labeled finite clipping floor.
    signed[good]=np.log2(np.maximum(values[good],1e-3))
    fig,ax=plt.subplots(figsize=(max(7,len(table.columns)*1.5),max(3.5,len(table.index)*.58+2)))
    cmap=plt.get_cmap("RdBu_r").copy();cmap.set_bad("lightgrey")
    image=ax.imshow(signed,aspect="auto",vmin=-2.5,vmax=2.5,cmap=cmap)
    ax.set_xticks(range(len(table.columns)),table.columns,rotation=45,ha="right")
    ax.set_yticks(range(len(table.index)),table.index)
    ax.set_title("Lipid enrichment/depletion near protein (replicate means)")
    ax.set_xlabel("Condition and leaflet")
    ax.set_ylabel("Lipid species")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i,j]):
                ax.text(j,i,f"{values[i,j]:.2f}×",ha="center",va="center",fontsize=9)
    cb=fig.colorbar(image,ax=ax,pad=.025);cb.set_label("log₂ enrichment (display clipped)")
    fig.tight_layout()
    dest=ensure_dir(Path(root)/"project_summary"/"figures"/"membrane")
    fig.savefig(dest/"lipid_enrichment_heatmap.png",dpi=int(cfg.get("plots",{}).get("dpi",250)))
    plt.close(fig)

