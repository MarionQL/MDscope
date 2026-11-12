"""Optional solvation figures from persisted numerical products."""
from __future__ import annotations

from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd


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
            if "mean_g_aqueous_reference" in curve:
                aqueous = curve.mean_g_aqueous_reference.to_numpy(float)
                if np.isfinite(aqueous).any():
                    ax.plot(x, aqueous, ls="--", label=f"{system}: aqueous-region reference")
        ax.axhline(1,color="gray",ls=":",lw=.8)
        ax.set_title(" / ".join(map(str,key)),fontsize=9)
        ax.set_ylabel("g(r); solid: full-cell reference")
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


def solvation_survival_figure(system_repouts, *, max_curves=8):
    """Show descriptive censored contact-survival curves, if recorded."""
    import pyarrow.parquet as pq
    curves=[]
    for system,reps in system_repouts.items():
        for rep in reps:
            path=Path(rep)/"extensions"/"solvation"/"survival.parquet"
            if not path.exists():
                continue
            cols=pq.ParquetFile(path).schema_arrow.names
            duration_column = "duration_ps" if "duration_ps" in cols else "time_ps"
            if not {duration_column,"survival_probability"}.issubset(cols):
                continue
            df=pd.read_parquet(path).rename(columns={duration_column: "duration_ps"})
            keys=[k for k in ("solvent","target_id") if k in df]
            groups=df.groupby(keys,sort=False) if keys else [((),df)]
            for key,group in groups:
                if len(curves)>=max_curves:
                    break
                curves.append((f"{system}/{rep.name}/{key}",group))
    if not curves:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.2,3.7))
    for label,df in curves:
        df=df.sort_values("duration_ps")
        ax.step(df.duration_ps,df.survival_probability,where="post",label=label,lw=1,marker=".")
    ax.set_xlim(left=0.)
    if all(curve.duration_ps.eq(0).all() for _, curve in curves):
        ax.set_xlim(0., 1.)
        ax.set_xticks([0.])
        ax.text(.5, .5, "Only zero-span sampled episodes;\nframe spacing limits lifetime resolution.",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
    ax.set_xlabel("Observed contact duration (ps)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0,1.02)
    ax.legend(fontsize=6,frameon=False)
    fig.tight_layout()
    buf=BytesIO()
    fig.savefig(buf,format="png",dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf,len(curves)



def plot_solvation(cfg, root, system_repouts):
    """Render a few solvation diagnostics when explicitly requested."""
    from ..config import ensure_dir, plot_enabled
    from .figure_inventory import current_statistics

    if (not cfg.get("plots", {}).get("families", {}).get("solvation", False)
            or not plot_enabled(cfg, "solvation")):
        return
    figures = {
        "shell_rdf": solvation_rdf_figure(root, system_repouts)
        if current_statistics(cfg, root, "solvation_rdfs") else None,
        "water_and_ions": solvation_timeseries_figure(system_repouts),
        "contact_survival": solvation_survival_figure(system_repouts),
    }
    for name, figure in figures.items():
        if figure is not None:
            directory = ensure_dir(Path(root) / "project_summary/figures/solvation")
            (directory / f"{name}.png").write_bytes(figure[0].getvalue())
