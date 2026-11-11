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
from ..interactions.interaction_tables import event_lifetimes
from ..interactions import interaction_key_columns, interaction_tables, sample_dt
from ..observables.membrane_order import temporal_membrane_order_profiles
from ..observables.density_1d import temporal_density_profiles
from ..observables.residue_membrane_distance import temporal_residue_membrane_profiles
from ..observables.structure_profiles import temporal_structure_profiles
from ..observables.rmsf import temporal_rmsf_profiles


from .common import finite, _temporal, _block_indices, _block_bootstrap, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta, _regular_times, _point_difference, _adaptive_intervals

LOCAL_GEOMETRY = {
    "sdc4_local_upper_height_A", "sdc4_local_upper_curvature_Ainv", "sdc4_local_thickness_A",
}


def _hs_contact_observables(series, group_id):
    """Derive occupancy and conditional counts even from older motif tables."""
    for name in series.columns:
        if not (name.startswith("sdc4_motif_") and name.endswith("_glycan_contacts")):
            continue
        count = pd.to_numeric(series[name], errors="coerce").where(lambda x: np.isfinite(x))
        prefix = name.removesuffix("glycan_contacts")
        yield (group_id, prefix + "hs_contacts_during_contact"), pd.DataFrame({
            "time_ps": series.time_ps, "value": count.where(count.gt(0)),
        })
        if prefix + "hs_any_contact_state" not in series:
            yield (group_id, prefix + "hs_any_contact_state"), pd.DataFrame({
                "time_ps": series.time_ps, "value": count.gt(0).astype(float).where(count.notna()),
            })


def _has_conditional_support(metric):
    return metric in LOCAL_GEOMETRY or (
        metric.startswith("sdc4_surface_site_") and metric.endswith("_leaflet_height_A")
    ) or metric.endswith("_hs_contacts_during_contact") or (
        metric.startswith("sdc4_order_")
        and metric.endswith(("near_minus_SCH", "far_minus_SCH", "near_minus_far_SCH"))
    )


def _conditional_geometry_bootstrap(values, times, cfg, draws, rng):
    """Preserve full-grid observation gaps when resampling local geometry."""
    present = np.isfinite(values)
    missing = np.full(draws, np.nan)
    if not _regular_times(times) or present.sum() < 2:
        return missing, missing, None, False
    center = values[present].mean()
    block = max(
        _temporal(present.astype(float), times, cfg)["block_length_frames"],
        _temporal(np.where(present, values - center, 0), times, cfg)["block_length_frames"],
    )
    # An interval needs support in more than one block. Withhold it if any
    # sampled full-length trace has no defined conditional estimate.
    if len(np.unique(np.flatnonzero(present) // block)) < 2:
        return missing, missing, block, False
    means, medians = np.empty(draws), np.empty(draws)
    for start in range(0, draws, 32):
        stop = min(start + 32, draws)
        sample = values[_block_indices(len(values), block, stop - start, rng)]
        if not np.isfinite(sample).any(axis=1).all():
            return missing, missing, block, False
        means[start:stop] = np.nanmean(sample, axis=1)
        medians[start:stop] = np.nanmedian(sample, axis=1)
    return means, medians, block, True


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
            out.update(_hs_contact_observables(subset, group_id))
            for metric in subset.columns:
                if (group_id is None and _has_conditional_support(metric)
                        and not metric.endswith("_hs_contacts_during_contact")):
                    available = pd.to_numeric(subset[metric], errors="coerce").map(np.isfinite)
                    out[(None, metric + "_available_fraction")] = pd.DataFrame({
                        "time_ps": subset.time_ps, "value": available.astype(float),
                    })
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
    bridges = Path(repout)/"extensions/solvation/bridge_pair_timeseries.parquet"
    if bridges.exists():
        import json
        data = pd.read_parquet(bridges)
        required = {"time_ps", "solvent", "target_a", "target_b", "bridge_entity_count", "bridge_present_fraction"}
        if not required <= set(data):
            raise ValueError(f"Incomplete bridge-pair time series: {bridges}")
        for (solvent, left, right), group in data.groupby(["solvent", "target_a", "target_b"], sort=False):
            if group.duplicated("time_ps").any():
                raise ValueError("Duplicate bridge-pair observation times")
            identity = "solvation_bridges/"+json.dumps([solvent, left, right], separators=(",", ":"))
            for metric in ("bridge_present_fraction", "bridge_entity_count",
                           "observed_minus_spatial_bridge_fraction"):
                if metric not in group or not group[metric].notna().any():
                    continue
                out[(identity, metric)] = group.sort_values("time_ps")[["time_ps", metric]].rename(columns={metric: "value"})
    matrix = Path(repout)/'extensions/lipid_contacts/matrix.parquet'
    if matrix.exists():
        data = pd.read_parquet(matrix)
        for (leaf, a, b), group in data.groupby(['leaflet','species_a','species_b'], sort=False):
            for metric in ('neighbors_per_molecule','enrichment'):
                out[(f'lipid_contacts/{leaf}/{a}/{b}', metric)] = group[['time_ps',metric]].rename(columns={metric:'value'})
    return out


def analyze_time_series(
    cfg: dict, root: Path, system_repouts: dict[str, list[Path]], *, output_name="time_series"
) -> None:
    cfg = statistics_family_config(cfg, "time_series")
    if cfg is None:
        return
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    scalar_window = stats.get("scalar_window_ps")
    draws, level = int(stats.get("bootstrap_draws", 2000)), float(stats.get("ci", 0.95))
    system_draws = draws if stats.get("hierarchy", {}).get("enabled", True) else 0
    seed = int(stats.get("seed", 2026))
    rng = np.random.default_rng(seed)
    adaptive = stats.get("adaptive_hierarchy") if system_draws else None
    adaptive_rng = np.random.default_rng(np.random.SeedSequence([seed, 1])) if adaptive else None
    resample = stats.get("hierarchy", {}).get("resample_replicates", True)
    adaptive_provenance = ({
        **{f"bootstrap_mc_{key}": value for key, value in adaptive.items()},
        "bootstrap_mc_seed": seed,
        "bootstrap_mc_first_stage_draws": draws,
    } if adaptive else {})
    outdir = ensure_dir(root / "project_summary" / "statistics" / output_name)
    rep_rows, acf_rows, cache, hs_episodes = [], [], {}, []
    for system, reps in system_repouts.items():
        for repout in reps:
            local_rows, local_acf = [], []
            for (group, metric), table in collect_time_series(repout).items():
                window_covered = True
                if scalar_window is not None:
                    lower, upper = scalar_window
                    window_covered = table.time_ps.min() <= lower and table.time_ps.max() >= upper
                    table = table.loc[table.time_ps.between(lower, upper)]
                values = pd.to_numeric(table["value"], errors="coerce").to_numpy(float)
                times = pd.to_numeric(table["time_ps"], errors="coerce").to_numpy(float)
                valid = np.isfinite(values) & np.isfinite(times)
                if (group is None and metric.startswith("sdc4_motif_") and "_hs_" in metric
                        and metric.endswith("_contact_state") and window_covered
                        and valid.all() and _regular_times(times) and np.isin(values, [0, 1]).all()):
                    occupied = pd.DataFrame({"time_ps": times[values == 1], "observable": metric})
                    gap_factor = cfg.get("analyses", {}).get("lifetimes", {}).get("gap_factor", 1.5)
                    episodes = event_lifetimes(occupied, ["observable"], times, gap_factor)
                    if not episodes.empty:
                        hs_episodes.append(episodes.assign(system=system, replicate=repout.name))
                temporal_valid = bool(window_covered and valid.all() and _regular_times(times))
                conditional = window_covered and group is None and _has_conditional_support(metric) and not temporal_valid
                if conditional:
                    conditional_boot = _conditional_geometry_bootstrap(
                        np.where(np.isfinite(values), values, np.nan), times, cfg, draws, rng
                    )
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
                ) if temporal_valid and draws else np.full(draws, np.nan)
                med_boot = _block_bootstrap(
                    values,
                    meta["block_length_frames"],
                    draws,
                    rng,
                    "median"
                ) if temporal_valid and draws else np.full(draws, np.nan)
                sampling_status = "regular_complete" if temporal_valid else "unavailable_gapped_or_irregular"
                if not window_covered:
                    sampling_status = "unavailable_window_coverage"
                if not temporal_valid:
                    meta.update(acf=np.array([]), tau_int_frames=np.nan,
                                tau_int_ps=np.nan, effective_n=np.nan)
                if conditional:
                    mean_boot, med_boot, block, temporal_valid = conditional_boot
                    if block is not None:
                        meta.update(block_length_frames=block,
                                    blocks_per_pseudostudy=int(np.ceil(len(valid) / block)))
                    sampling_status = (
                        "regular_conditional_on_availability" if temporal_valid
                        else "unavailable_conditional_block_support"
                    )
                mlo, mhi = _ci(mean_boot, level)
                dlo, dhi = _ci(med_boot, level)
                row = {
                    "system": system,
                    "replicate": repout.name,
                    "temporal_sampling_status": sampling_status,
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
            if local_rows and output_name == "time_series":
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
    if hs_episodes:
        atomic_parquet(outdir / "hs_contact_episodes.parquet", pd.concat(hs_episodes, ignore_index=True))
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
        require_complete = scalar_window is not None or _has_conditional_support(metric)
        complete_replica_coverage = not require_complete or len(keys) == len(system_repouts[system])
        temporal_valid = complete_replica_coverage and all(cache[k]["temporal_valid"] for k in keys)
        if adaptive:
            def sample(n):
                return np.column_stack([
                    _hierarchical([cache[k][stat] for k in keys], n, adaptive_rng, stat, resample)
                    for stat in ("mean", "median")
                ]) if temporal_valid else np.empty((0, 2))

            (mean_boot, med_boot), precision = _adaptive_intervals(sample, level, adaptive)
        else:
            mean_boot, med_boot = (
                _hierarchical([cache[k][stat] for k in keys], system_draws, rng, stat, resample)
                if temporal_valid else np.full(system_draws, np.nan)
                for stat in ("mean", "median")
            )
            precision = {}
        sys_boot[(system, gid, metric)] = {
            "temporal_valid": temporal_valid,
            **({"keys": keys} if adaptive else {"mean": mean_boot, "median": med_boot}),
        }
        mlo, mhi = _ci(mean_boot, level)
        dlo, dhi = _ci(med_boot, level)
        sampling = sorted(set(sub.temporal_sampling_status))
        if not complete_replica_coverage:
            sampling.append("incomplete_replica_coverage")
        sys_rows.append({
            "system": system,
            "group": gid,
            "metric": metric,
            "n_replicates": len(sub),
            "temporal_sampling_status": "|".join(sampling),
            "mean": float(sub["mean"].mean()),
                         "median": float(sub["median"].median()), "mean_ci_low": mlo,
                         "mean_ci_high": mhi,
                         "median_ci_low": dlo, "median_ci_high": dhi, **_range_meta(
                             rows
                         ), "bootstrap_draws": len(mean_boot), "ci_level": level,
            **adaptive_provenance, **precision})
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
            left = sys_boot[(s1, gid, metric)]
            right = sys_boot[(s2, gid, metric)]
            if adaptive:
                def sample(n):
                    if not regular or (paired and not amap.keys() & bmap.keys()):
                        return np.empty((0, 2))
                    columns = []
                    for stat in ("mean", "median"):
                        a = {k[1]: cache[k][stat] for k in left["keys"]}
                        b = {k[1]: cache[k][stat] for k in right["keys"]}
                        if paired:
                            values = _hierarchical_difference(a, b, n, adaptive_rng, stat, resample)
                        else:
                            values = (_hierarchical(b.values(), n, adaptive_rng, stat, resample)
                                      - _hierarchical(a.values(), n, adaptive_rng, stat, resample))
                        columns.append(values)
                    return np.column_stack(columns)

                (dmean, dmed), precision = _adaptive_intervals(sample, level, adaptive)
            elif not regular:
                # A paired nan-reducer would silently omit invalid replicas.
                dmean = np.full(system_draws, np.nan)
                dmed = np.full(system_draws, np.nan)
            elif paired:
                dmean, dmed = (
                    _hierarchical_difference(
                        {rep_id: cache[(s1, rep_id, gid, metric)][stat] for rep_id in amap},
                        {rep_id: cache[(s2, rep_id, gid, metric)][stat] for rep_id in bmap},
                        system_draws, rng, stat, resample)
                    for stat in ("mean", "median")
                )
            else:
                dmean = right["mean"] - left["mean"]
                dmed = right["median"] - left["median"]
            mlo, mhi = _ci(dmean, level)
            dlo, dhi = _ci(dmed, level)
            pgt, plt = _direction(dmean)
            dpgt, dplt = _direction(dmed)
            if regular and system_draws and len(dmean):
                mean_difference, median_difference = float(np.mean(dmean)), float(np.median(dmed))
            else:
                mean_difference = _point_difference(amap,bmap,paired)
                median_difference = _point_difference(amed,bmed,paired,"median")
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
                          "ci_status": "unavailable_temporal_support" if not regular else None,
                          "bootstrap_draws": len(dmean) if adaptive else system_draws,
                          "ci_level": level, **adaptive_provenance, **precision})
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



