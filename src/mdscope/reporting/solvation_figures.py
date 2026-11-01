"""Persisted solvation evidence for the compact PDF report.

Detailed, potentially huge tables remain in the companion review ZIP.
No trajectory or neighbor-search calculations occur here.
"""
from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd


SUMMARY_TABLES = (
    "target_summary", "shell_candidates", "shell_stability",
    "survival", "episodes", "bridges", "water_hbond_bridges",
    "hydration_sasa", "water_density_3d",
)


def solvation_inventory(system_repouts):
    rows = []
    for system, reps in system_repouts.items():
        for rep in reps:
            directory = Path(rep)/"extensions"/"solvation"
            if not directory.is_dir():
                continue
            for path in sorted(directory.glob("*.parquet")):
                # Inspect Parquet metadata, never materialize huge event files.
                import pyarrow.parquet as pq
                rows.append({
                    "system":system, "replicate":Path(rep).name,
                    "table":path.stem, "rows":pq.ParquetFile(path).metadata.num_rows,
                    "path":path})
    return pd.DataFrame(rows)


def solvation_rdf_figure(root, system_repouts, *, max_panels=6):
    path=Path(root)/"project_summary"/"statistics"/"solvation"/"system_rdf.parquet"
    if not path.exists():
        return None
    rdf=pd.read_parquet(path)
    if rdf.empty:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups=list(rdf.groupby(["solvent","target_class","site_class"],sort=True))
    groups=groups[:max_panels]
    n=len(groups)
    fig,axes=plt.subplots(n,1,figsize=(7.2,max(2.4,2.15*n)),squeeze=False)
    for ax,(key,group) in zip(axes.flat,groups):
        for system,curve in group.groupby("system"):
            curve=curve.sort_values("r_inner_A")
            x=(curve.r_inner_A.to_numpy(float)+curve.r_outer_A.to_numpy(float))/2
            y=curve.mean_g_bulk.to_numpy(float)
            sd=curve.sd_g_bulk.to_numpy(float)
            ax.plot(x,y,label=f"{system} (n={int(curve.n_replicates.min())})")
            ax.fill_between(x,y-sd,y+sd,alpha=.18)
        ax.axhline(1,color="gray",ls=":",lw=.8)
        ax.set_title(" / ".join(map(str,key)),fontsize=9)
        ax.set_ylabel("g(r)")
        ax.legend(fontsize=7,frameon=False)
    axes[-1,0].set_xlabel("Distance (Å)")
    fig.tight_layout()
    buf=BytesIO()
    fig.savefig(buf,format="png",dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf,len(groups),len(list(rdf.groupby(["solvent","target_class","site_class"])))


def solvation_timeseries_figure(system_repouts, *, max_panels=4):
    data=[]
    for system,reps in system_repouts.items():
        for rep in reps:
            path=Path(rep)/"extensions"/"solvation"/"species_timeseries.parquet"
            if path.exists():
                tab=pd.read_parquet(path)
                if {"time_ps","solvent","target_class","site_class",
                    "mean_n_solvent"}.issubset(tab):
                    data.append(tab.assign(system=system,replicate=Path(rep).name))
    if not data:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df=pd.concat(data,ignore_index=True)
    groups=list(df.groupby(["solvent","target_class","site_class"],sort=True))
    selected=groups[:max_panels]
    fig,axes=plt.subplots(len(selected),1,figsize=(7.2,max(2.6,2.25*len(selected))),
                          squeeze=False)
    for ax,(key,group) in zip(axes.flat,selected):
        for (system,rep),curve in group.groupby(["system","replicate"]):
            curve=curve.sort_values("time_ps")
            ax.plot(curve.time_ps/1000,curve.mean_n_solvent,
                    label=f"{system}/{rep}",lw=.9)
        ax.set_title(" / ".join(map(str,key)),fontsize=9)
        ax.set_ylabel("Mean solvent count")
        ax.legend(fontsize=6,frameon=False,ncol=2)
    axes[-1,0].set_xlabel("Time (ns)")
    fig.tight_layout()
    buf=BytesIO()
    fig.savefig(buf,format="png",dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf,len(selected),len(groups)


def solvation_stat_rows(root):
    base=Path(root)/"project_summary"/"statistics"/"time_series"
    p=base/"system_statistics.parquet"
    if not p.exists():
        return pd.DataFrame()
    df=pd.read_parquet(p)
    if "group" not in df:
        return pd.DataFrame()
    return df[df.group.fillna("").astype(str).str.startswith("solvation/")]


def solvation_comparison_rows(root):
    p=Path(root)/"project_summary"/"statistics"/"time_series"/"system_comparisons.parquet"
    if not p.exists():
        return pd.DataFrame()
    df=pd.read_parquet(p)
    if "group" not in df:
        return pd.DataFrame()
    return df[df.group.fillna("").astype(str).str.startswith("solvation/")]
