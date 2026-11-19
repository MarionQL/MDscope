"""
Shared autocorrelation, hierarchical bootstrap, permutation, FDR, and
comparison statistics.
"""

from __future__ import annotations

from time import perf_counter
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import statistics_family_config
from ..statistics_selection import interaction_selection, selected_interaction_metrics, metric_estimator_enabled
from ..core.persistence import atomic_json, atomic_parquet, ensure_dir
from ..interactions.interaction_tables import interaction_key_columns, interaction_tables, sample_dt


from .common import _temporal, _block_weights, _hierarchical, _hierarchical_difference, _ci, _direction, _permutation, _fdr, _range_meta, _point_difference, selected_rng as _selected_rng


def _interaction_labels(table: pd.DataFrame, keys: list[str]) -> np.ndarray:
    """Format a label once per unique identity, preserving persisted key semantics."""
    if not keys:
        return np.full(len(table), "interaction", dtype=object)
    # Normalize absent keys before factorization: pandas may materialize
    # Python None, NaN, or pd.NA for the same saved nullable field.
    columns = [
        np.frompyfunc(
            lambda value: "nan" if pd.isna(value) else str(value), 1, 1
        )(table[key].to_numpy(dtype=object))
        for key in keys
    ]
    codes, unique = pd.factorize(pd.MultiIndex.from_arrays(columns), sort=False)
    labels = np.asarray([" | ".join(row) for row in unique], dtype=object)
    return labels[codes]


def _presence_frame_indices(times, presence_times):
    indices = np.searchsorted(times, np.asarray(presence_times, float))
    valid = (indices < len(times)) & np.isclose(
        times[np.minimum(indices, len(times) - 1)], presence_times
    )
    return indices, valid


def _interaction_bootstrap(
    times,
    presence_times,
    lifetimes: pd.DataFrame,
    cfg: dict,
    draws: int,
    rng: np.random.Generator,
    cache: dict,
    rep_key,
    zero_meta_cache: dict | None = None,
    selected_seed: int | None = None,
) -> tuple[dict, dict]:
    # Bootstrap the binary presence trajectory. Event starts/lifetimes
    # are weighted at their original frames so independently sampled
    # blocks cannot create fake joined events.
    times = np.asarray(times, float)
    n = len(times)
    dt = sample_dt(times)
    presence = np.zeros(n, np.float32)
    if n and len(presence_times):
        idx, valid = _presence_frame_indices(times, presence_times)
        presence[idx[valid]] = 1
    # Interactions absent from a replicate produce the same zero-occupancy
    # temporal metadata for every identity on that replicate's time grid.
    # Cache only that deterministic metadata: bootstrap weight generation,
    # draws, and hierarchical RNG order remain untouched.
    absent = not len(presence_times) and (lifetimes is None or lifetimes.empty)
    if absent and zero_meta_cache is not None and rep_key in zero_meta_cache:
        meta = zero_meta_cache[rep_key]
    else:
        meta = _temporal(presence, times, cfg)
        if absent and zero_meta_cache is not None:
            zero_meta_cache[rep_key] = meta
    key = (rep_key, n, meta["block_length_frames"], draws)
    # Reuse the same block weights for all interactions from one
    # replicate/block length; this is a major speed/memory win when
    # thousands of pairs are present.
    if key not in cache:
        if len(cache) >= 8:
            cache.pop(next(iter(cache)))
        weight_rng = (rng if selected_seed is None else
                      _selected_rng(selected_seed, "block", *rep_key, n, key[2]))
        cache[key] = _block_weights(n, key[2], draws, weight_rng) if draws else np.empty((0, n), np.uint16)
    weights = cache[key]
    duration_ns = n * dt / 1000 if dt else np.nan
    occupancy = weights.dot(presence) / max(1, n)
    starts = np.zeros(n, np.float32)
    complete_frames, durations = [], []
    if lifetimes is not None and not lifetimes.empty:
        event_times = lifetimes["start_ps"].to_numpy(float)
        indices, valid = _presence_frame_indices(times, event_times) if n else (
            np.empty(len(event_times), dtype=int), np.zeros(len(event_times), dtype=bool)
        )
        if "left_censored" in lifetimes:
            valid &= ~lifetimes["left_censored"].to_numpy(dtype=bool)
        np.add.at(starts, indices[valid], 1)
        completed = (lifetimes["completed"].to_numpy(dtype=bool)
                     if "completed" in lifetimes else np.zeros(len(event_times), dtype=bool))
        keep = valid & completed
        complete_frames = indices[keep].tolist()
        durations = lifetimes.loc[keep, "lifetime_ps"].to_numpy(float).tolist()
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
    cfg = statistics_family_config(cfg, "interactions")
    if cfg is None:
        return
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    selectors = interaction_selection(cfg)
    if selectors is not None and not selectors:
        return
    min_occupancy = float(stats.get("interaction_min_occupancy", 0.0))
    if not 0 <= min_occupancy <= 1:
        raise ValueError("statistics.interaction_min_occupancy must be between 0 and 1")
    draws, level = int(stats.get("bootstrap_draws", 2000)), float(stats.get("ci", 0.95))
    system_draws = draws if stats.get("hierarchy", {}).get("enabled", True) else 0
    seed = int(stats.get("seed", 2026)) + 2
    rng = np.random.default_rng(seed)
    resample = stats.get("hierarchy", {}).get("resample_replicates", True)
    # No output directory is created until saved identities are validated.
    outdir = root / "project_summary" / "statistics" / "interactions"
    sources, kinds = {}, set()
    detail_enabled = os.environ.get("MDSCOPE_PROFILE", "0").strip().lower() not in {
        "", "0", "false", "off", "no",
    }
    detail = {"kinds": []}
    loading_started = perf_counter()
    excluded_kinds = set(stats.get('interaction_exclude_kinds', []))
    included = (set(stats["interaction_include_kinds"])
                if stats.get("interaction_include_kinds") is not None else None)
    if selectors is not None:
        wanted_kinds = {kind for kind, _, _ in selectors} - excluded_kinds
        included = wanted_kinds if included is None else included & wanted_kinds
    verbose = cfg.get("progress", {}).get("enabled", True)
    for system, reps in system_repouts.items():
        for repout in reps:
            if verbose:
                print(f"Interaction statistics: loading {system}/{repout.name}", flush=True)
            tables = interaction_tables(
                repout,
                exclude_kinds=excluded_kinds,
                table_names=("presence", "lifetimes"),
                include_kinds=included,
            )
            path = repout / "summary" / "timeseries.parquet"
            times = (
                pd.read_parquet(path, columns=["time_ps"])["time_ps"].to_numpy(float)
                if path.exists() else np.array([], float)
            )
            sources[(system, repout.name)] = (times, tables)
            kinds.update(tables)
    if selectors is not None:
        requested = {(kind, identity) for kind, identity, _ in selectors}
        existing = {(kind, None) for kind in kinds - excluded_kinds}
        # Kind-wide selections need only directory discovery, not another
        # 10k-identity string materialization. Exact labels are checked
        # against the saved canonical key columns before any output write.
        explicit = {(kind, identity) for kind, identity in requested
                    if identity is not None}
        for kind in {kind for kind, _ in explicit} & (kinds - excluded_kinds):
            wanted = {identity for group, identity in explicit if group == kind}
            for _, tables in sources.values():
                tab = tables.get(kind, {})
                for name in ("presence", "lifetimes"):
                    data = tab.get(name)
                    if data is None or data.empty:
                        continue
                    keys = interaction_key_columns(kind, data)
                    missing_columns = [key for key in keys if key not in data.columns]
                    if missing_columns:
                        raise ValueError(
                            f"Interaction {kind} missing identity columns: {missing_columns}")
                    labels = _interaction_labels(data, keys)
                    existing.update((kind, identity) for identity in wanted
                                    if identity in set(labels))
        missing = requested - existing
        if missing:
            raise ValueError(
                "Unknown saved interaction identities or excluded kinds: "
                + ", ".join(map(str, sorted(missing, key=str)))
            )
    ensure_dir(outdir)
    if detail_enabled:
        detail["input_loading_seconds"] = perf_counter() - loading_started
    rep_rows, sys_rows, comps = [], [], []
    weight_cache = {}
    zero_meta_cache = {}
    for kind in sorted(kinds - excluded_kinds):
        if verbose:
            print(f"Interaction statistics {kind}: grouping input records", flush=True)
        grouping_started = perf_counter()
        interactions, retained, prepared = set(), set(), {}
        requested_ids = (None if selectors is None or any(
            group == kind and identity is None for group, identity, _ in selectors
        ) else {identity for group, identity, _ in selectors if group == kind})
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
                    table["interaction"] = _interaction_labels(table, keys)
                    if requested_ids is not None:
                        table.drop(table.index[~table["interaction"].isin(requested_ids)],
                                   inplace=True)
                    interactions.update(table["interaction"])
            if min_occupancy and len(times) and not presence.empty:
                indices, valid = _presence_frame_indices(times, presence["time_ps"])
                occupied_frames = pd.Series(
                    indices[valid], index=presence.loc[valid, "interaction"].to_numpy()
                ).groupby(level=0, sort=False).nunique()
                retained.update(occupied_frames.index[occupied_frames / len(times) >= min_occupancy])
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
        if min_occupancy:
            if verbose:
                print(
                    f"Interaction statistics {kind}: retaining {len(retained)}/{len(interactions)} "
                    f"identities with occupancy >= {min_occupancy:.1%} in at least one replicate",
                    flush=True,
                )
            interactions = retained
        if detail_enabled:
            grouped_seconds = perf_counter() - grouping_started
        kind_started = perf_counter()
        log_every = int(stats.get("interaction_progress_every", 1000))
        for completed, interaction in enumerate(sorted(interactions), 1):
            selected_metrics = selected_interaction_metrics(selectors, kind, interaction)
            if not selected_metrics:
                continue
            metric_infer = {
                metric: metric_estimator_enabled(
                    cfg, "interactions", (kind, interaction, metric), "ci"
                ) for metric in selected_metrics
            }
            metric_tests = {
                metric: metric_estimator_enabled(
                    cfg, "interactions", (kind, interaction, metric), "hypothesis_tests"
                ) for metric in selected_metrics
            }
            infer_any = any(metric_infer.values())
            local_draws = draws if infer_any else 0
            local_system_draws = system_draws if infer_any else 0
            if verbose and (completed == 1 or completed % log_every == 0):
                print(f"Interaction statistics {kind}: processing identity "
                      f"{completed}/{len(interactions)} | elapsed {perf_counter()-kind_started:.1f} s",
                      flush=True)
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
                        local_draws,
                        rng,
                        weight_cache,
                        (system, repout.name),
                        zero_meta_cache,
                        seed if selectors is not None else None,
                    )
                    repdata[(system, repout.name)] = (boot, meta)
                    for metric, values in boot.items():
                        if metric not in selected_metrics:
                            continue
                        lo, hi = _ci(values, level) if metric_infer[metric] else (np.nan, np.nan)
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
                                         "bootstrap_draws": local_draws if metric_infer[metric] else 0, "ci_level": level})
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
                    if not valid or metric not in selected_metrics:
                        continue
                    infer_metric = metric_infer[metric]
                    combine = "median" if metric == "median_lifetime_ps" else "mean"
                    boot = _hierarchical(
                        [x[1][0][metric] for x in valid],
                        local_system_draws if infer_metric else 0,
                        (rng if selectors is None else
                         _selected_rng(seed, "hierarchy", kind, interaction, system, metric)),
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
                                     ), "bootstrap_draws": local_system_draws if infer_metric else 0, "ci_level": level})
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
                            local_system_draws if metric_infer[metric] else 0,
                            (rng if selectors is None else
                             _selected_rng(seed, "paired", kind, interaction, s1, s2, metric)),
                            reducer,
                            resample,
                        )
                    else:
                        delta = sys_boot[s2][metric] - sys_boot[s1][metric]
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
                                  "difference_system2_minus_system1": (float(np.nanmean(delta)) if local_system_draws and metric_infer[metric] and len(delta)
                                      else _point_difference(amap, bmap, paired,
                                            "median" if metric == "median_lifetime_ps" else "mean")),
                                  "difference_ci_low": lo, "difference_ci_high": hi,
                                  "probability_difference_gt_0": pgt,
                                  "probability_difference_lt_0": plt,
                                  "p_value": (_permutation(
                                      amap, bmap, cfg,
                                      rng if selectors is None else
                                      _selected_rng(seed, "permutation", kind, interaction,
                                                    s1, s2, metric))
                                      if metric_tests[metric] else np.nan),

                                  "bootstrap_draws": local_system_draws if metric_infer[metric] else 0, "ci_level": level})
        if detail_enabled:
            detail["kinds"].append({
                "kind": kind, "unique_identities": len(interactions),
                "grouping_seconds": grouped_seconds,
                "per_identity_seconds": perf_counter() - kind_started,
            })
    writing_started = perf_counter()
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

    if detail_enabled:
        detail["persistence_seconds"] = perf_counter() - writing_started
        detail["bootstrap_draws"] = draws
        detail["total_unique_identities_across_kinds"] = sum(
            entry["unique_identities"] for entry in detail["kinds"]
        )
        atomic_json(outdir / "interaction_phase_profile.json", detail)



