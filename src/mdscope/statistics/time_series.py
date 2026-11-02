"""
Shared autocorrelation, hierarchical bootstrap, permutation, FDR, and
comparison statistics.
"""

from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import statistics_family_config, atomic_parquet, ensure_dir
from ..interactions import interaction_key_columns, interaction_tables, sample_dt
from ..observables.membrane_order import temporal_membrane_order_profiles
from ..observables.density_1d import temporal_density_profiles
from ..observables.residue_membrane_distance import temporal_residue_membrane_profiles
from ..observables.structure_profiles import temporal_structure_profiles
from ..observables.rmsf import temporal_rmsf_profiles


from .common import finite, _temporal, _block_bootstrap, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta

def collect_time_series(repout: Path) -> dict[tuple[str | None, str], pd.DataFrame]:
    out = {}
    for path, group_col in(
        (Path(repout) / "summary" / "timeseries.parquet", None),
        (Path(repout) / "summary" / "group_timeseries.parquet", "group")
    ):
        if not path.exists():
            continue
        data = pd.read_parquet(path)
        groups = [(None, data)] if group_col is None else data.groupby(group_col, sort=False)
        for group_id, subset in groups:
            for metric in subset.columns:
                if metric not in {"time_ps", "group"} and len(finite(subset[metric])):
                    out[(
                        None if group_id is None else str(group_id),
                        metric
                    )] = subset[["time_ps", metric]].rename(columns={metric: "value"})
    # Reuse the SAME temporal block bootstrap / hierarchical statistics as
    # every other MDscope time series. Never resample individual lipid molecules
    # or water contacts as if they were independent simulation replicates.
    solvation = Path(repout) / "extensions" / "solvation" / "species_timeseries.parquet"
    if solvation.exists():
        data = pd.read_parquet(solvation)
        required = {"solvent", "target_class", "site_class", "time_ps",
                    "mean_n_solvent"}
        if not required.issubset(data.columns):
            raise ValueError(
                "solvation species_timeseries.parquet is incomplete: "
                + ", ".join(sorted(required - set(data.columns))))
        for (solvent, kind, label), group in data.groupby(
                ["solvent", "target_class", "site_class"], sort=False):
            group_id = f"solvation/{solvent}/{kind}/{label}"
            if group.duplicated("time_ps").any():
                raise ValueError(f"duplicate solvation time points: {group_id}")
            group = group.sort_values("time_ps")
            for metric in ("mean_n_solvent", "hydrated_fraction"):
                if metric in group and len(finite(group[metric])):
                    out[(group_id, metric)] = group[["time_ps", metric]].rename(
                        columns={metric: "value"})
    return out


def analyze_time_series(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]) -> None:
    cfg = statistics_family_config(cfg, "time_series")
    if cfg is None:
        return
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    draws, level = int(stats.get("bootstrap_draws", 2000)), float(stats.get("ci", 0.95))
    rng = np.random.default_rng(int(stats.get("seed", 2026)))
    resample = stats.get("hierarchy", {}).get("resample_replicates", True)
    outdir = ensure_dir(root / "project_summary" / "statistics" / "time_series")
    rep_rows, acf_rows, cache = [], [], {}
    for system, reps in system_repouts.items():
        for repout in reps:
            local_rows, local_acf = [], []
            for (group, metric), table in collect_time_series(repout).items():
                values = pd.to_numeric(table["value"], errors="coerce").to_numpy(float)
                times = pd.to_numeric(table["time_ps"], errors="coerce").to_numpy(float)
                valid = np.isfinite(values) & np.isfinite(times)
                intervals = np.diff(times)
                temporal_valid = bool(valid.all() and len(intervals) > 0
                    and np.all(intervals > 0)
                    and np.allclose(intervals, intervals[0], rtol=1e-5, atol=1e-8))
                values = values[valid]
                times = times[valid]
                if not len(values):
                    continue
                meta = _temporal(values, times, cfg)
                mean_boot = _block_bootstrap(
                    values,
                    meta["block_length_frames"],
                    draws,
                    rng,
                    "mean"
                ) if temporal_valid else np.full(draws, np.nan)
                med_boot = _block_bootstrap(
                    values,
                    meta["block_length_frames"],
                    draws,
                    rng,
                    "median"
                ) if temporal_valid else np.full(draws, np.nan)
                if not temporal_valid:
                    meta.update(acf=np.array([]), tau_int_frames=np.nan,
                                tau_int_ps=np.nan, effective_n=np.nan)
                mlo, mhi = _ci(mean_boot, level)
                dlo, dhi = _ci(med_boot, level)
                row = {
                    "system": system,
                    "replicate": repout.name,
                    "temporal_sampling_status": "regular_complete" if temporal_valid else "unavailable_gapped_or_irregular",
                    "group": group,
                    "metric": metric,
                    **_describe(values),
                       **{
                           k: v for k,
                           v in meta.items() if k != "acf"
                       }, "bootstrap_draws": draws, "ci_level": level,
                    "mean_ci_low": mlo,
                    "mean_ci_high": mhi,
                    "median_ci_low": dlo,
                    "median_ci_high": dhi,
                }
                rep_rows.append(row)
                local_rows.append(row)
                cache[(system, repout.name, group, metric)] = {
                    "mean": mean_boot,
                    "median": med_boot,
                    "temporal_valid": temporal_valid
                }
                for lag, value in enumerate(meta["acf"]):
                    a = {
                        "system": system,
                        "replicate": repout.name,
                        "group": group,
                        "metric": metric,
                        "lag_frames": lag,
                         "lag_ps": lag * meta["sample_dt_ps"], "autocorrelation": float(value)}
                    acf_rows.append(a)
                    local_acf.append(a)
            if local_rows:
                atomic_parquet(
                    repout / "summary" / "time_series_statistics.parquet",
                    pd.DataFrame(local_rows)
                )
                atomic_parquet(
                    repout / "summary" / "autocorrelation.parquet",
                    pd.DataFrame(local_acf)
                )
    rep = pd.DataFrame(rep_rows)
    if rep.empty:
        return
    atomic_parquet(outdir / "replicate_statistics.parquet", rep)
    atomic_parquet(outdir / "autocorrelation.parquet", pd.DataFrame(acf_rows))
    sys_rows, sys_boot = [], {}
    for(
        system,
        group,
        metric
    ), sub in rep.groupby(["system", "group", "metric"], dropna=False, sort=False):
        gid = None if pd.isna(group) else str(group)
        rows = sub.to_dict("records")
        keys = [(system, r["replicate"], gid, metric) for r in rows]
        temporal_valid = all(cache[k]["temporal_valid"] for k in keys)
        mean_boot = (_hierarchical([cache[k]["mean"] for k in keys], draws, rng, "mean", resample)
                     if temporal_valid else np.full(draws, np.nan))
        med_boot = _hierarchical(
            [cache[k]["median"] for k in keys],
            draws,
            rng,
            "median",
            resample
        ) if temporal_valid else np.full(draws, np.nan)
        sys_boot[(system, gid, metric)] = {"mean": mean_boot, "median": med_boot, "temporal_valid": temporal_valid}
        mlo, mhi = _ci(mean_boot, level)
        dlo, dhi = _ci(med_boot, level)
        sys_rows.append({
            "system": system,
            "group": gid,
            "metric": metric,
            "n_replicates": len(sub),
            "temporal_sampling_status": "regular_complete" if temporal_valid else "unavailable_gapped_or_irregular",
            "mean": float(sub["mean"].mean()),
                         "median": float(sub["median"].median()), "mean_ci_low": mlo,
                         "mean_ci_high": mhi,
                         "median_ci_low": dlo, "median_ci_high": dhi, **_range_meta(
                             rows
                         ), "bootstrap_draws": draws, "ci_level": level})
    sys = pd.DataFrame(sys_rows)
    atomic_parquet(outdir / "system_statistics.parquet", sys)
    comps = []
    for group, metric in sys[["group", "metric"]].drop_duplicates().itertuples(
        index=False,
        name=None
    ):
        gid = None if pd.isna(group) else str(group)
        available = [s for s in system_repouts if (s, gid, metric) in sys_boot]
        for s1, s2 in combinations(available, 2):
            sub = rep[
                (rep.metric == metric)
                & (rep.group.fillna("") == (gid or ""))
            ]
            amap = {
                str(row.replicate): float(row.mean)
                for row in sub[sub.system == s1].itertuples()
            }
            bmap = {
                str(row.replicate): float(row.mean)
                for row in sub[sub.system == s2].itertuples()
            }
            amed = {
                str(row.replicate): float(row.median)
                for row in sub[sub.system == s1].itertuples()
            }
            bmed = {
                str(row.replicate): float(row.median)
                for row in sub[sub.system == s2].itertuples()
            }
            paired = (
                stats.get("hypothesis_tests", {}).get("design", "independent")
                == "paired"
            )
            regular = (sys_boot[(s1, gid, metric)]["temporal_valid"]
                       and sys_boot[(s2, gid, metric)]["temporal_valid"])
            if not regular:
                # A paired nan-reducer would silently omit invalid replicas.
                # Neither design supports temporal intervals after that loss.
                dmean = np.full(draws, np.nan)
                dmed = np.full(draws, np.nan)
            elif paired:
                mean_draws_a = {
                    rep_id: cache[(s1, rep_id, gid, metric)]["mean"]
                    for rep_id in amap
                }
                mean_draws_b = {
                    rep_id: cache[(s2, rep_id, gid, metric)]["mean"]
                    for rep_id in bmap
                }
                median_draws_a = {
                    rep_id: cache[(s1, rep_id, gid, metric)]["median"]
                    for rep_id in amed
                }
                median_draws_b = {
                    rep_id: cache[(s2, rep_id, gid, metric)]["median"]
                    for rep_id in bmed
                }
                dmean = _hierarchical_difference(
                    mean_draws_a,
                    mean_draws_b,
                    draws,
                    rng,
                    "mean",
                    resample,
                )
                dmed = _hierarchical_difference(
                    median_draws_a,
                    median_draws_b,
                    draws,
                    rng,
                    "median",
                    resample,
                )
            else:
                dmean = (
                    sys_boot[(s2, gid, metric)]["mean"]
                    - sys_boot[(s1, gid, metric)]["mean"]
                )
                dmed = (
                    sys_boot[(s2, gid, metric)]["median"]
                    - sys_boot[(s1, gid, metric)]["median"]
                )
            mlo, mhi = _ci(dmean, level)
            dlo, dhi = _ci(dmed, level)
            pgt, plt = _direction(dmean)
            dpgt, dplt = _direction(dmed)
            if regular:
                mean_difference, median_difference = float(np.mean(dmean)), float(np.median(dmed))
            elif paired:
                labels = sorted(set(amap) & set(bmap))
                mean_difference = float(np.mean([bmap[k]-amap[k] for k in labels]))
                median_difference = float(np.median([bmed[k]-amed[k] for k in labels]))
            else:
                mean_difference = float(np.mean(list(bmap.values())) - np.mean(list(amap.values())))
                median_difference = float(np.median(list(bmed.values())) - np.median(list(amed.values())))
            comps.append({"system_1": s1, "system_2": s2, "group": gid, "metric": metric,
                          "mean_difference_system2_minus_system1": mean_difference,
                          "mean_difference_ci_low": mlo, "mean_difference_ci_high": mhi,
                          "mean_probability_difference_gt_0": pgt,
                          "mean_probability_difference_lt_0": plt,
                          "mean_p_value": _permutation(amap, bmap, cfg, rng),
                          "median_difference_system2_minus_system1": median_difference,
                          "median_difference_ci_low": dlo, "median_difference_ci_high": dhi,
                          "median_probability_difference_gt_0": dpgt,
                          "median_probability_difference_lt_0": dplt,
                          "median_p_value": _permutation(amed, bmed, cfg, rng),
                          "bootstrap_draws": draws, "ci_level": level})
    comp = pd.DataFrame(comps)
    if not comp.empty:
        family_columns = ["system_1", "system_2"]
        if stats.get("hypothesis_families"):
            membership = {metric: name for name, metrics in stats["hypothesis_families"].items() for metric in metrics}
            comp["hypothesis_family"] = comp.metric.map(membership).fillna("exploratory")
            family_columns.append("hypothesis_family")
        if stats.get("hypothesis_tests", {}).get("multiple_testing", "fdr_bh") == "fdr_bh":
            comp = _fdr(comp, "mean_p_value", "mean_q_value", family_columns)
            comp = _fdr(comp, "median_p_value", "median_q_value", family_columns)
        atomic_parquet(outdir / "system_comparisons.parquet", comp)



