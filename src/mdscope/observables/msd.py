"""Existing membrane mean-squared displacement, isolated from APL/thickness.

Streamed from the same trajectory pass and checkpoint state as before.
Diffusion is a descriptive single-origin MSD slope fit on each replicate.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from ..core.geometry import minimum_image_displacements
from MDAnalysis.exceptions import NoDataError
from ..config import atomic_parquet
from .multi_origin_msd import diffusion_fit_sensitivity

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
        step = np.asarray(minimum_image_displacements(
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

def summarize_project_diffusion(root:Path,system_repouts:dict):
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


def summarize_project_multi_origin(root, system_repouts, cfg=None):
    """Aggregate independent per-replicate multi-origin fits without frame pseudoreplication."""
    settings = (cfg or {}).get("analyses", {}).get("membrane", {}).get("msd", {})
    diagnose = settings.get("multi_origin_fit_diagnostics", False)
    diagnostics = []
    tables = []
    for system, reps in system_repouts.items():
        for repout in reps:
            path = Path(repout) / "summary" / "multi_origin_diffusion.parquet"
            if path.exists():
                fits = pd.read_parquet(path)
                tables.append(fits)
            else:
                fits = None
            if diagnose:
                curve_path = path.with_name("multi_origin_msd.parquet")
                curves = pd.read_parquet(curve_path) if curve_path.is_file() else None
                metrics = [] if fits is None else fits["metric"].tolist()
                for metric in metrics or [None]:
                    identity = {"system": system, "replicate": Path(repout).name, "metric": metric}
                    curve = None if curves is None or metric is None else curves.loc[curves.metric == metric]
                    if curve is None or curve.empty:
                        diagnostics.append(dict(identity, fit_status="missing_curve_or_fit", screen="unavailable"))
                        continue
                    curve = curve.sort_values("lag_ps")
                    for row in diffusion_fit_sensitivity(
                        curve.lag_ps, curve.msd_A2, curve.n_time_origins,
                        start=float(settings.get("slope_start_fraction", .2)),
                        end=float(settings.get("slope_end_fraction", .8)),
                    ):
                        diagnostics.append(identity | row)
    target = Path(root) / "project_summary" / "diffusion"
    diagnostic_path = target / "multi_origin_fit_diagnostics.parquet"
    if diagnostics:
        target.mkdir(parents=True, exist_ok=True)
        atomic_parquet(diagnostic_path, pd.DataFrame(diagnostics))
    else:
        diagnostic_path.unlink(missing_ok=True)
    if not tables:
        return
    all_rows = pd.concat(tables, ignore_index=True)
    target.mkdir(parents=True, exist_ok=True)
    atomic_parquet(target / "multi_origin_replicates.parquet", all_rows)
    grouped = (all_rows.groupby(["system", "metric"], sort=True)["diffusion_nm2_ns"]
               .agg(n_replicates="count", mean="mean", sd="std", median="median")
               .reset_index())
    grouped["interpretation"] = "replicate_level_descriptive_no_significance"
    atomic_parquet(target / "multi_origin_systems.parquet", grouped)

