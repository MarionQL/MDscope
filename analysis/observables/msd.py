"""Existing membrane mean-squared displacement, isolated from APL/thickness.

Streamed from the same trajectory pass and checkpoint state as before.
Diffusion is a descriptive single-origin MSD slope fit on each replicate.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.exceptions import NoDataError
from ..config import atomic_parquet

def _residue_centers(atoms):
    """Vectorized center per selected residue, preserving MDAnalysis residue order.

    MDAnalysis implements compound centers in NumPy rather than making one
    Python AtomGroup/COM call per lipid per frame. Preserve the original
    center-of-geometry fallback for missing or invalid masses.
    """
    try:
        result = np.asarray(atoms.center_of_mass(compound="residues"), dtype=float)
    except (AttributeError, ValueError, TypeError, NoDataError):
        result = np.asarray(atoms.center_of_geometry(compound="residues"), dtype=float)
    invalid = ~np.all(np.isfinite(result), axis=1)
    if np.any(invalid):
        geometric = np.asarray(atoms.center_of_geometry(compound="residues"), dtype=float)
        result[invalid] = geometric[invalid]
    return result

def initialize_msd(universe,state,spec):
    msd = spec.get("msd", {}) or {}
    state.msd = bool(msd.get("enabled", False))
    state.msd_include_3d = bool(msd.get("include_3d", False))
    state.msd_remove_drift = bool(msd.get("remove_drift", False))
    if state.msd:
        if universe.dimensions is None or not np.all(np.asarray(
            universe.dimensions[:3],
            float
        ) > 0):
            raise ValueError("membrane MSD requires valid periodic box dimensions")
        atoms = universe.select_atoms(str(msd.get("selection") or spec["headgroup_selection"]))
        if not len(atoms):
            raise ValueError("membrane.msd.selection matched zero atoms")
        state.msd_atoms = atoms
        state.msd_resindices = np.asarray(atoms.residues.ix, dtype=int)
        state.msd_upper_mask = np.isin(state.msd_resindices, state.upper_resix)
        state.msd_lower_mask = np.isin(state.msd_resindices, state.lower_resix)
        if not np.all(state.msd_upper_mask | state.msd_lower_mask):
            raise ValueError(
                "membrane.msd.selection includes residues outside the leaflet classification"
            )
        state.msd_last_positions = _residue_centers(state.msd_atoms)
        state.msd_displacements = np.zeros_like(state.msd_last_positions)

def msd_metrics(universe,state):
    out={}
    if state.msd:
        current = _residue_centers(state.msd_atoms)
        # Accumulate minimum-image step vectors so lipids crossing
        # periodic boundaries do not create artificial MSD jumps.
        step = np.asarray(minimize_vectors(
            current - state.msd_last_positions,
            universe.dimensions
        ), dtype=float)
        if state.msd_remove_drift:
            # Subtract common translation in the lipid reference frame.
            step -= step.mean(axis=0, keepdims=True)
        state.msd_displacements += step
        state.msd_last_positions = current
        lateral = [i for i in range(3) if i != state.ai]
        lateral_sq = np.sum(state.msd_displacements[:, lateral] ** 2, axis=1)
        out["membrane_msd_lateral_A2"] = float(np.mean(lateral_sq))
        out["membrane_msd_lateral_upper_A2"] = float(np.mean(lateral_sq[state.msd_upper_mask]))
        out["membrane_msd_lateral_lower_A2"] = float(np.mean(lateral_sq[state.msd_lower_mask]))
        if state.msd_include_3d:
            total_sq = np.sum(state.msd_displacements ** 2, axis=1)
            out["membrane_msd_3d_A2"] = float(np.mean(total_sq))
            out["membrane_msd_3d_upper_A2"] = float(np.mean(total_sq[state.msd_upper_mask]))
            out["membrane_msd_3d_lower_A2"] = float(np.mean(total_sq[state.msd_lower_mask]))
    return out

def fit_lateral_diffusion(times_ps,msd_A2,start_fraction=.2,end_fraction=.8):
    """D[nm²/ns] = slope[Å²/ps] * (1000 ps/ns)/(100 Å²/nm²)/4."""
    t=np.asarray(times_ps,float);y=np.asarray(msd_A2,float)
    good=np.isfinite(t)&np.isfinite(y)
    t,y=t[good],y[good]
    if len(t)<6 or t[-1]<=t[0]:
        return {"diffusion_nm2_ns":np.nan,"r_squared":np.nan,"n_fit":0}
    begin=t[0]+(t[-1]-t[0])*float(start_fraction)
    end=t[0]+(t[-1]-t[0])*float(end_fraction)
    mask=(t>=begin)&(t<=end)
    if mask.sum()<4:
        return {"diffusion_nm2_ns":np.nan,"r_squared":np.nan,"n_fit":int(mask.sum())}
    slope,intercept=np.polyfit(t[mask],y[mask],1)
    pred=slope*t[mask]+intercept
    ss=float(np.sum((y[mask]-pred)**2))
    total=float(np.sum((y[mask]-y[mask].mean())**2))
    return {"diffusion_nm2_ns":float(slope*2.5),
            "r_squared":float(1-ss/total) if total>0 else np.nan,
            "n_fit":int(mask.sum()),"slope_A2_ps":float(slope),
            "fit_start_ps":float(begin),"fit_end_ps":float(end)}

def summarize_diffusion(repout:Path,cfg):
    """Only consume saved MDscope membrane MSD; never reopen trajectories."""
    path=Path(repout)/"summary"/"timeseries.parquet"
    if not path.exists():return
    data=pd.read_parquet(path)
    if "time_ps" not in data:return
    settings=cfg.get("analyses",{}).get("membrane",{}).get("msd",{})
    begin=float(settings.get("slope_start_fraction",.2))
    end=float(settings.get("slope_end_fraction",.8))
    if not 0<=begin<end<=1:
        raise ValueError("membrane.msd fit fractions require 0 <= start < end <= 1")
    rows=[]
    for col in data.columns:
        if col.startswith("membrane_msd_lateral") and col.endswith("_A2"):
            stats=fit_lateral_diffusion(data["time_ps"],data[col],begin,end)
            rows.append({"metric":col,"replicate":Path(repout).name,"system":Path(repout).parent.name,
                         "method":"single_origin_msd_slope","fit_start_fraction":begin,
                         "fit_end_fraction":end,**stats})
    if rows:atomic_parquet(Path(repout)/"summary"/"diffusion_coefficients.parquet",pd.DataFrame(rows))

def summarize_project_diffusion(root:Path,system_repouts:dict,cfg):
    """Collect independent replica slope estimates and create publication-ready D tables."""
    chunks=[]
    for system,reps in system_repouts.items():
        for repout in reps:
            f=Path(repout)/"summary"/"diffusion_coefficients.parquet"
            if f.exists():chunks.append(pd.read_parquet(f))
    if not chunks:return
    all_data=pd.concat(chunks,ignore_index=True)
    out=Path(root)/"project_summary"/"diffusion"
    out.mkdir(parents=True,exist_ok=True)
    atomic_parquet(out/"replicate_coefficients.parquet",all_data)
    grouped=(all_data.groupby(["system","metric"],sort=True)["diffusion_nm2_ns"]
             .agg(n_replicates="count",mean="mean",sd="std",median="median",minimum="min",maximum="max")
             .reset_index())
    grouped["interpretation"]="replicate_level_descriptive_no_significance"
    atomic_parquet(out/"system_coefficients.parquet",grouped)
    # Coefficients are replicate-derived summaries, unlike framewise observables:
    # show independent replicate points; never violin-plot frames as replicas.
    if cfg.get("plots",{}).get("enabled",True):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(max(6,1.4*len(grouped)),4.8))
        labels=[];ticks=[]
        for i,((system,metric),data) in enumerate(all_data.groupby(["system","metric"],sort=True)):
            vals=data["diffusion_nm2_ns"].dropna().to_numpy(float)
            if not len(vals):continue
            ax.scatter(np.full(len(vals),i),vals,s=45,alpha=.85)
            ax.plot([i-.22,i+.22],[vals.mean(),vals.mean()],lw=2.5,color="black")
            labels.append(system+"\n"+metric.replace("membrane_msd_lateral","lateral"))
            ticks.append(i)
        ax.set_xticks(ticks,labels,rotation=25,ha="right")
        ax.set_ylabel("Lateral diffusion coefficient (nm²/ns)")
        ax.set_title("Diffusion from existing single-origin MSD: replicate estimates")
        fig.tight_layout()
        fig.savefig(out/"diffusion_coefficients.png",dpi=int(cfg.get("plots",{}).get("dpi",250)))
        plt.close(fig)
