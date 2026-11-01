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

from ..config import atomic_parquet, ensure_dir
from ..interactions import interaction_key_columns, interaction_tables, sample_dt
from ..observables.membrane_order import temporal_membrane_order_profiles
from ..observables.density_1d import temporal_density_profiles
from ..observables.residue_membrane_distance import temporal_residue_membrane_profiles
from ..observables.structure_profiles import temporal_structure_profiles
from ..observables.rmsf import temporal_rmsf_profiles


from .common import _temporal, _block_weights, _hierarchical, _hierarchical_difference, _ci, _direction, _permutation, _fdr, _range_meta

def _interaction_bootstrap(
    times,
    presence_times,
    lifetimes: pd.DataFrame,
    cfg: dict,
    draws: int,
    rng: np.random.Generator,
    cache: dict,
    rep_key
) -> tuple[dict, dict]:
    # Bootstrap the binary presence trajectory. Event starts/lifetimes
    # are weighted at their original frames so independently sampled
    # blocks cannot create fake joined events.
    times = np.asarray(times, float)
    n = len(times)
    dt = sample_dt(times)
    presence = np.zeros(n, np.float32)
    if len(presence_times):
        idx = np.searchsorted(times, np.asarray(presence_times, float))
        valid = (idx < n) & np.isclose(times[np.minimum(idx, n - 1)], presence_times)
        presence[idx[valid]] = 1
    meta = _temporal(presence, times, cfg)
    key = (rep_key, n, meta["block_length_frames"])
    # Reuse the same block weights for all interactions from one
    # replicate/block length; this is a major speed/memory win when
    # thousands of pairs are present.
    if key not in cache:
        if len(cache) >= 8:
            cache.pop(next(iter(cache)))
        cache[key] = _block_weights(n, key[2], draws, rng)
    weights = cache[key]
    duration_ns = n * dt / 1000 if dt else np.nan
    occupancy = weights.dot(presence) / max(1, n)
    starts = np.zeros(n, np.float32)
    complete_frames, durations = [], []
    if lifetimes is not None and not lifetimes.empty:
        for row in lifetimes.itertuples():
            if bool(getattr(row, "left_censored", False)):
                continue
            i = int(np.searchsorted(times, float(row.start_ps)))
            if i >= n or not np.isclose(times[i], float(row.start_ps)):
                continue
            starts[i] += 1
            if bool(getattr(row, "completed", False)):
                complete_frames.append(i)
                durations.append(float(row.lifetime_ps))
    event_rate = (
        weights.dot(starts) / duration_ns if np.isfinite(duration_ns)
        and duration_ns > 0 else np.full(draws, np.nan)
    )
    if complete_frames:
        ew = weights[:, np.asarray(complete_frames, int)].astype(float)
        durations = np.asarray(durations, float)
        counts = ew.sum(axis=1)
        mean_life = np.divide(
            ew.dot(durations),
            counts,
            out=np.full(draws, np.nan),
            where=counts > 0
        )
        order = np.argsort(durations)
        sorted_w, sorted_d = ew[:, order], durations[order]
        totals = sorted_w.sum(axis=1)
        med_life = sorted_d[np.argmax(np.cumsum(sorted_w, axis=1) >= totals[:, None] / 2, axis=1)]
        med_life[totals == 0] = np.nan
    else:
        durations = np.array([], float)
        mean_life = med_life = np.full(draws, np.nan)
    point = {"occupancy": float(presence.mean()) if n else np.nan,
             "event_rate_per_ns": float(starts.sum() / duration_ns) if np.isfinite(duration_ns)
             and duration_ns > 0 else np.nan,
             "mean_lifetime_ps": float(np.mean(durations)) if len(durations) else np.nan,
             "median_lifetime_ps": float(np.median(durations)) if len(durations) else np.nan,
             "total_lifetime_ps": float(
                 presence.sum() * dt
             ) if dt else np.nan, "n_event_starts": int(
                 starts.sum()
             ), "n_completed_events": len(durations)}
    boot = {
        "occupancy": occupancy,
        "event_rate_per_ns": event_rate,
        "mean_lifetime_ps": mean_life,
        "median_lifetime_ps": med_life,
            "total_lifetime_ps": occupancy * n * dt if dt else np.full(draws, np.nan)}
    return boot, {**point, **{k: v for k, v in meta.items() if k != "acf"}}


def analyze_interactions(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]) -> None:
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    draws, level = int(stats.get("bootstrap_draws", 2000)), float(stats.get("ci", 0.95))
    rng = np.random.default_rng(int(stats.get("seed", 2026)) + 2)
    resample = stats.get("hierarchy", {}).get("resample_replicates", True)
    outdir = ensure_dir(root / "project_summary" / "statistics" / "interactions")
    sources, kinds = {}, set()
    excluded_kinds = set(stats.get('interaction_exclude_kinds', []))
    for system, reps in system_repouts.items():
        for repout in reps:
            tables = interaction_tables(repout)
            path = repout / "summary" / "timeseries.parquet"
            times = (
                pd.read_parquet(
                    path
                )["time_ps"].to_numpy(float) if path.exists() else np.array([], float)
            )
            sources[(system, repout.name)] = (times, tables)
            kinds.update(tables)
    rep_rows, sys_rows, comps = [], [], []
    weight_cache = {}
    for kind in sorted(kinds - excluded_kinds):
        interactions, prepared = set(), {}
        for key, (times, tables) in sources.items():
            tab = tables.get(kind, {})
            presence = tab.get("presence", pd.DataFrame()).copy()
            life = tab.get("lifetimes", pd.DataFrame()).copy()
            base = presence if not presence.empty else life
            keys = interaction_key_columns(kind, base)
            missing = [column for column in keys if column not in base.columns]
            if missing:
                raise ValueError(
                    f"Interaction table {kind!r} is missing identity columns: "
                    + ", ".join(missing)
                )
            for table in (presence, life):
                if not table.empty:
                    table["interaction"] = (
                        table[keys].astype(str).agg(" | ".join, axis=1)
                        if keys
                        else "interaction"
                    )
                    interactions.update(table["interaction"])
            presence_groups = (
                {identity: group["time_ps"].to_numpy(float)
                 for identity, group in presence.groupby("interaction", sort=False)}
                if not presence.empty else {}
            )
            life_groups = (
                {identity: group for identity, group in life.groupby("interaction", sort=False)}
                if not life.empty else {}
            )
            prepared[key] = (times, presence_groups, life_groups)
        for interaction in sorted(interactions):
            repdata = {}
            for system, reps in system_repouts.items():
                for repout in reps:
                    times, presence_groups, life_groups = prepared[(system, repout.name)]
                    if not len(times):
                        continue
                    pt = presence_groups.get(interaction)
                    if pt is None:
                        pt = np.array([], float)
                    ls = life_groups.get(interaction)
                    if ls is None:
                        ls = pd.DataFrame()
                    boot, meta = _interaction_bootstrap(
                        times,
                        pt,
                        ls,
                        cfg,
                        draws,
                        rng,
                        weight_cache,
                        (system, repout.name)
                    )
                    repdata[(system, repout.name)] = (boot, meta)
                    for metric, values in boot.items():
                        lo, hi = _ci(values, level)
                        rep_rows.append({
                            "system": system,
                            "replicate": repout.name,
                            "interaction_type": kind,
                            "interaction": interaction,
                            "metric": metric,
                                         "value": meta[metric], "ci_low": lo, "ci_high": hi,
                                         **{k: meta[k] for k in(
                                             "tau_int_frames",
                                             "tau_int_ps",
                                             "effective_n",
                                             "block_length_frames",
                                             "block_length_ps",
                                             "blocks_per_pseudostudy",
                                             "n_event_starts",
                                             "n_completed_events"
                                         )},
                                         "bootstrap_draws": draws, "ci_level": level})
            sys_boot = {}
            for system, reps in system_repouts.items():
                items = [(
                    r.name,
                    repdata[(system, r.name)]
                ) for r in reps if (system, r.name) in repdata]
                if not items:
                    continue
                sys_boot[system] = {}
                for metric in(
                    "occupancy",
                    "event_rate_per_ns",
                    "mean_lifetime_ps",
                    "median_lifetime_ps",
                    "total_lifetime_ps"
                ):
                    valid = [(rep, data) for rep, data in items if np.isfinite(data[1][metric])]
                    if not valid:
                        continue
                    combine = "median" if metric == "median_lifetime_ps" else "mean"
                    boot = _hierarchical(
                        [x[1][0][metric] for x in valid],
                        draws,
                        rng,
                        combine,
                        resample
                    )
                    sys_boot[system][metric] = boot
                    lo, hi = _ci(boot, level)
                    points = np.asarray([x[1][1][metric] for x in valid], float)
                    sys_rows.append({
                        "system": system,
                        "interaction_type": kind,
                        "interaction": interaction,
                        "metric": metric,
                                     "n_replicates_total": len(items),
                                     "n_replicates_observed": len(valid),
                                     "mean": float(np.mean(points)),
                                     "median": float(np.median(points)),
                                     "mean_ci_low": lo, "mean_ci_high": hi, **_range_meta(
                                         [x[1][1] for x in valid]
                                     ), "bootstrap_draws": draws, "ci_level": level})
            available = [s for s in system_repouts if s in sys_boot]
            for s1, s2 in combinations(available, 2):
                for metric in sorted(set(sys_boot[s1]) & set(sys_boot[s2])):
                    paired = (
                        stats.get("hypothesis_tests", {}).get(
                            "design",
                            "independent",
                        )
                        == "paired"
                    )
                    if paired:
                        rep_draws_1 = {
                            rep: data[0][metric]
                            for (system, rep), data in repdata.items()
                            if system == s1 and np.isfinite(data[1][metric])
                        }
                        rep_draws_2 = {
                            rep: data[0][metric]
                            for (system, rep), data in repdata.items()
                            if system == s2 and np.isfinite(data[1][metric])
                        }
                        reducer = (
                            "median"
                            if metric == "median_lifetime_ps"
                            else "mean"
                        )
                        delta = _hierarchical_difference(
                            rep_draws_1,
                            rep_draws_2,
                            draws,
                            rng,
                            reducer,
                            resample,
                        )
                    else:
                        delta = sys_boot[s2][metric] - sys_boot[s1][metric]
                    if not len(delta):
                        continue
                    lo, hi = _ci(delta, level)
                    pgt, plt = _direction(delta)
                    amap = {
                        rep: data[1][metric] for (sys, rep),
                        data in repdata.items() if sys == s1 and np.isfinite(data[1][metric])
                    }
                    bmap = {
                        rep: data[1][metric] for (sys, rep),
                        data in repdata.items() if sys == s2 and np.isfinite(data[1][metric])
                    }
                    comps.append({
                        "system_1": s1,
                        "system_2": s2,
                        "interaction_type": kind,
                        "interaction": interaction,
                        "metric": metric,
                                  "difference_system2_minus_system1": float(np.nanmean(delta)),
                                  "difference_ci_low": lo, "difference_ci_high": hi,
                                  "probability_difference_gt_0": pgt,
                                  "probability_difference_lt_0": plt,
                                  "p_value": _permutation(amap, bmap, cfg, rng),
                                  "bootstrap_draws": draws, "ci_level": level})
    for name, rows in(
        ("replicate_interaction_metrics.parquet", rep_rows),
        ("system_interaction_metrics.parquet", sys_rows)
    ):
        if rows:
            atomic_parquet(outdir / name, pd.DataFrame(rows))
    comp = pd.DataFrame(comps)
    if not comp.empty:
        if stats.get(
            "hypothesis_tests",
            {}
        ).get(
            "multiple_testing",
            "fdr_bh"
        ) == "fdr_bh": comp = _fdr(
            comp,
            "p_value",
            "q_value",
            ["system_1", "system_2", "interaction_type", "metric"]
        )
        atomic_parquet(outdir / "system_interaction_comparisons.parquet", comp)


