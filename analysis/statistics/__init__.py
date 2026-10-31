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
from ..membrane_order import temporal_membrane_order_profiles
from ..timeseries import temporal_density_profiles
from ..profiles import temporal_residue_membrane_profiles, temporal_structure_profiles
from ..timeseries import temporal_rmsf_profiles


def finite(values) -> np.ndarray:
    x = np.asarray(pd.to_numeric(values, errors="coerce"), float)
    return x[np.isfinite(x)]


def autocorrelation(values, max_lag=None) -> np.ndarray:
    # FFT autocorrelation is O(N log N), much faster than correlating
    # every lag directly for long MD trajectories.
    x = finite(values)
    if len(x) < 2:
        return np.array([1.0])
    x -= x.mean()
    variance = np.dot(x, x)
    if variance == 0:
        return np.ones(1)
    lag = min(len(x) - 1, int(max_lag) if max_lag is not None else max(1, min(len(x) // 2, 5000)))
    size = 1 << (2 * len(x) - 1).bit_length()
    fft = np.fft.rfft(x, size)
    acf = np.fft.irfft(
        fft * np.conjugate(fft),
        size
    )[:lag + 1] / np.arange(len(x), len(x) - lag - 1, -1)
    return acf / acf[0]


def _temporal(values, times, cfg: dict) -> dict:
    acf = autocorrelation(values, cfg["statistics"].get("block", {}).get("max_acf_lag"))
    positive = acf[1:]
    stop = np.flatnonzero(positive <= 0)
    tau = max(0.5, float(0.5 + np.sum(positive[:stop[0]] if len(stop) else positive)))
    n = len(finite(values))
    block_cfg = cfg["statistics"].get("block", {}) or {}
    requested = block_cfg.get("length", "auto")
    hierarchy = cfg["statistics"].get("hierarchy", {})
    if not hierarchy.get("block_within_replicate", True):
        block = 1
    elif requested != "auto":
        block = max(1, int(requested))
    else:
        multiplier = float(block_cfg.get("tau_multiplier", 2.0))
        block = max(1, int(np.ceil(multiplier * tau)))
    block = min(max(1, n), block) if n else 1
    dt = sample_dt(times)
    return {
        "acf": acf, "tau_int_frames": tau, "tau_int_ps": tau * dt if dt else np.nan,
        "effective_n": n / (2 * tau) if n else 0.0,
        "block_length_frames": block, "block_length_ps": block * dt if dt else np.nan,
        "blocks_per_pseudostudy": int(np.ceil(n / block)) if n else 0, "sample_dt_ps": dt,
    }


def _block_indices(n: int, block: int, draws: int, rng: np.random.Generator) -> np.ndarray:
    # A pseudo-study is rebuilt from contiguous blocks so bootstrap
    # draws preserve short-range temporal correlation.
    if n < 1:
        return np.empty((draws, 0), np.int32)
    block = min(max(1, int(block)), n)
    nblocks = int(np.ceil(n / block))
    starts = rng.integers(
        0,
        n - block + 1,
        size=(draws, nblocks),
        dtype=np.int32
    ) if n > block else np.zeros((draws, nblocks), np.int32)
    return (starts[..., None] + np.arange(block, dtype=np.int32)).reshape(draws, -1)[:, :n]


def _block_bootstrap(
    values,
    block: int,
    draws: int,
    rng: np.random.Generator,
    reducer="mean"
) -> np.ndarray:
    values = np.asarray(values, float)
    shape = (draws,) + values.shape[1:]
    out = np.empty(shape, float)
    reduce = np.nanmedian if reducer == "median" else np.nanmean
    # Generate draws in small batches so profile bootstraps do not
    # materialize drawsxframesxfeatures all at once.
    batch = 32
    for start in range(0, draws, batch):
        size = min(batch, draws - start)
        out[start:start + size] = reduce(
            values[_block_indices(len(values), block, size, rng)],
            axis=1
        )
    return out


def _block_weights(n: int, block: int, draws: int, rng: np.random.Generator) -> np.ndarray:
    # Interaction statistics need event locations as well as occupancy.
    # Frame multiplicities reproduce a block bootstrap without copying
    # binary trajectories.
    idx = _block_indices(n, block, draws, rng)
    dtype = np.uint16 if int(np.ceil(n / max(1, block))) < 65535 else np.uint32
    return np.stack([np.bincount(row, minlength=n) for row in idx]).astype(dtype, copy=False)


def _hierarchical(
    rep_draws,
    draws: int,
    rng: np.random.Generator,
    reducer="mean",
    resample=True
) -> np.ndarray:
    # Second bootstrap level: choose replicates, then choose one
    # temporal-bootstrap realization from each selected replicate.
    arrays = [np.asarray(x, float) for x in rep_draws if len(x)]
    if not arrays:
        return np.array([])
    stack = np.stack(arrays)
    nrep = len(stack)
    reps = rng.integers(
        0,
        nrep,
        size=(draws, nrep)
    ) if resample else np.broadcast_to(np.arange(nrep), (draws, nrep))
    draw_idx = rng.integers(0, stack.shape[1], size=(draws, nrep))
    selected = stack[reps, draw_idx]
    return (np.nanmedian if reducer == "median" else np.nanmean)(selected, axis=1)


def _hierarchical_difference(
    rep_draws_a: dict[str, np.ndarray],
    rep_draws_b: dict[str, np.ndarray],
    draws: int,
    rng: np.random.Generator,
    reducer: str = "mean",
    resample: bool = True,
) -> np.ndarray:
    """Bootstrap a paired system difference while preserving replicate pairs."""
    labels = sorted(set(rep_draws_a) & set(rep_draws_b))
    if not labels:
        return np.array([])
    stack_a = np.stack([np.asarray(rep_draws_a[label], float) for label in labels])
    stack_b = np.stack([np.asarray(rep_draws_b[label], float) for label in labels])
    if stack_a.shape[2:] != stack_b.shape[2:]:
        raise ValueError("Paired bootstrap arrays have incompatible shapes")
    nrep = len(labels)
    selected_reps = (
        rng.integers(0, nrep, size=(draws, nrep))
        if resample
        else np.broadcast_to(np.arange(nrep), (draws, nrep))
    )
    draw_a = rng.integers(0, stack_a.shape[1], size=(draws, nrep))
    draw_b = rng.integers(0, stack_b.shape[1], size=(draws, nrep))
    differences = (
        stack_b[selected_reps, draw_b]
        - stack_a[selected_reps, draw_a]
    )
    reduce = np.nanmedian if reducer == "median" else np.nanmean
    return reduce(differences, axis=1)


def _ci(draws, level: float) -> tuple[float, float]:
    x = finite(draws)
    if not len(x):
        return np.nan, np.nan
    alpha = (1 - level) / 2
    return tuple(map(float, np.quantile(x, [alpha, 1 - alpha])))


def _describe(values) -> dict:
    x = finite(values)
    if not len(x):
        return {}
    q1, med, q3 = np.quantile(x, [0.25, 0.5, 0.75])
    return {
        "n": len(x),
        "mean": x.mean(),
        "sd": x.std(ddof=1) if len(x) > 1 else np.nan,
        "median": med,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "min": x.min(),
        "max": x.max(),
    }


def _direction(delta) -> tuple[float, float]:
    x = finite(delta)
    return (float(np.mean(x > 0)), float(np.mean(x < 0))) if len(x) else (np.nan, np.nan)


def _permutation(
    a: dict[str, float],
    b: dict[str, float],
    cfg: dict,
    rng: np.random.Generator,
    reducer="mean"
) -> float:
    # Replicates-not trajectory frames-are the independent units for
    # hypothesis testing. Exact enumeration is used when feasible.
    tests = cfg["statistics"].get("hypothesis_tests", {}) or {}
    if not tests.get("enabled", True):
        return np.nan
    func = np.mean if reducer == "mean" else np.median
    exact_limit = int(tests.get("exact_max_assignments", 100000))
    nperm = int(tests.get("permutation_draws", 10000))
    if tests.get("design", "independent") == "paired":
        labels = sorted(set(a) & set(b))
        d = finite([b[k] - a[k] for k in labels])
        if not len(d):
            return np.nan
        observed, total = abs(float(func(d))), 2 ** len(d)
        if total <= exact_limit:
            values = [func(d * np.where(
                ((mask >> np.arange(len(d))) & 1),
                1.0,
                -1.0
            )) for mask in range(total)]
            return float(np.mean(np.abs(values) >= observed - 1e-15))
        values = func(rng.choice([-1.0, 1.0], size=(nperm, len(d))) * d, axis=1)
        return float((1 + np.sum(np.abs(values) >= observed)) / (nperm + 1))
    av, bv = finite(list(a.values())), finite(list(b.values()))
    if not len(av) or not len(bv):
        return np.nan
    pooled, na = np.concatenate([av, bv]), len(av)
    observed = abs(float(func(bv) - func(av)))
    total = math.comb(len(pooled), na)
    if total <= exact_limit:
        values = []
        for chosen in combinations(range(len(pooled)), na):
            mask = np.zeros(len(pooled), bool)
            mask[list(chosen)] = True
            values.append(func(pooled[~mask]) - func(pooled[mask]))
        return float(np.mean(np.abs(values) >= observed - 1e-15))
    extreme = 0
    for _ in range(nperm):
        order = rng.permutation(len(pooled))
        extreme += abs(float(func(pooled[order[na:]]) - func(pooled[order[:na]]))) >= observed
    return float((extreme + 1) / (nperm + 1))


def _fdr(table: pd.DataFrame, p_col: str, q_col: str, families: list[str]) -> pd.DataFrame:
    if table.empty or p_col not in table:
        return table
    table = table.copy()
    table[q_col] = np.nan
    for _, idx in table.groupby(families, dropna=False, sort=False).groups.items():
        p = table.loc[idx, p_col].to_numpy(float)
        valid = np.flatnonzero(np.isfinite(p))
        if not len(valid):
            continue
        order = np.argsort(p[valid])
        ranked = p[valid][order]
        q = np.minimum.accumulate((ranked * len(ranked) / np.arange(
            1,
            len(ranked) + 1
        ))[::-1])[::-1]
        restored = np.empty_like(q)
        restored[order] = np.clip(q, 0, 1)
        table.loc[np.asarray(idx)[valid], q_col] = restored
    return table


def _range_meta(rows) -> dict:
    out = {}
    for key in(
        "tau_int_frames",
        "tau_int_ps",
        "effective_n",
        "block_length_frames",
        "block_length_ps",
        "blocks_per_pseudostudy"
    ):
        x = finite([row[key] for row in rows if key in row])
        if len(x):
            out.update({
                f"{key}_median": float(np.median(x)),
                f"{key}_min": float(np.min(x)),
                f"{key}_max": float(np.max(x)),
            })
    return out


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
    return out


def analyze_time_series(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]) -> None:
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
                )
                med_boot = _block_bootstrap(
                    values,
                    meta["block_length_frames"],
                    draws,
                    rng,
                    "median"
                )
                mlo, mhi = _ci(mean_boot, level)
                dlo, dhi = _ci(med_boot, level)
                row = {
                    "system": system,
                    "replicate": repout.name,
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
                    "median": med_boot
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
        mean_boot = _hierarchical([cache[k]["mean"] for k in keys], draws, rng, "mean", resample)
        med_boot = _hierarchical(
            [cache[k]["median"] for k in keys],
            draws,
            rng,
            "median",
            resample
        )
        sys_boot[(system, gid, metric)] = {"mean": mean_boot, "median": med_boot}
        mlo, mhi = _ci(mean_boot, level)
        dlo, dhi = _ci(med_boot, level)
        sys_rows.append({
            "system": system,
            "group": gid,
            "metric": metric,
            "n_replicates": len(sub),
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
            if paired:
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
            comps.append({"system_1": s1, "system_2": s2, "group": gid, "metric": metric,
                          "mean_difference_system2_minus_system1": float(np.mean(dmean)),
                          "mean_difference_ci_low": mlo, "mean_difference_ci_high": mhi,
                          "mean_probability_difference_gt_0": pgt,
                          "mean_probability_difference_lt_0": plt,
                          "mean_p_value": _permutation(amap, bmap, cfg, rng),
                          "median_difference_system2_minus_system1": float(np.median(dmed)),
                          "median_difference_ci_low": dlo, "median_difference_ci_high": dhi,
                          "median_probability_difference_gt_0": dpgt,
                          "median_probability_difference_lt_0": dplt,
                          "median_p_value": _permutation(amed, bmed, cfg, rng),
                          "bootstrap_draws": draws, "ci_level": level})
    comp = pd.DataFrame(comps)
    if not comp.empty:
        if stats.get("hypothesis_tests", {}).get("multiple_testing", "fdr_bh") == "fdr_bh":
            comp = _fdr(comp, "mean_p_value", "mean_q_value", ["system_1", "system_2"])
            comp = _fdr(comp, "median_p_value", "median_q_value", ["system_1", "system_2"])
        atomic_parquet(outdir / "system_comparisons.parquet", comp)


def _rmsf_profile(
    coords: np.ndarray,
    residues,
    residue_indices=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return residue RMSF in first-occurrence/topology order."""
    residue_labels = np.asarray(residues, str)
    grouping = (
        np.asarray(residue_indices, int)
        if residue_indices is not None
        else residue_labels
    )
    residx, unique_groups = pd.factorize(grouping, sort=False)
    labels = np.asarray(
        [residue_labels[np.flatnonzero(residx == i)[0]] for i in range(len(unique_groups))],
        dtype=str,
    )
    if len(coords):
        var = np.var(coords, axis=0, ddof=0)
        atom = np.sqrt(np.sum(var, axis=1))
    else:
        atom = np.full(coords.shape[1], np.nan)
    counts = np.bincount(residx, minlength=len(labels))
    residue_rmsf = np.divide(
        np.bincount(residx, weights=atom, minlength=len(labels)),
        counts,
        out=np.full(len(labels), np.nan, dtype=float),
        where=counts > 0,
    )
    return residue_rmsf, labels, residx


def _prepare_profile(spec: dict, cfg: dict, draws: int, rng: np.random.Generator) -> dict:
    times = np.asarray(spec["times"], float)
    if spec.get("reducer") != "rmsf":
        values = np.asarray(spec["values"], float)
        reducer = str(spec.get("reducer", "mean"))
        signal = np.sqrt(np.nanmean((values - np.nanmean(values, axis=0)) ** 2, axis=1))
        meta = _temporal(np.nan_to_num(signal), times, cfg)
        return {
            "positions": np.asarray(spec["positions"], float),
            "labels": np.asarray(spec["labels"], str),
                "point": (np.nanmedian if reducer == "median" else np.nanmean)(values, axis=0),
                "draws": _block_bootstrap(
                    values,
                    meta["block_length_frames"],
                    draws,
                    rng,
                    reducer
                ), "temporal": meta}
    coords = np.asarray(spec["coords"], np.float32)
    residues = np.asarray(spec["atom_residues"], str)
    signal = np.sqrt(np.mean(np.sum((coords - np.mean(coords, axis=0)) ** 2, axis=2), axis=1))
    meta = _temporal(signal, times, cfg)
    point, labels, residx = _rmsf_profile(
        coords,
        residues,
        spec.get("atom_residue_indices"),
    )
    weights = _block_weights(
        len(coords),
        meta["block_length_frames"],
        draws,
        rng
    ).astype(np.float32)
    mean = (
        np.einsum("dn,nac->dac", weights, coords, optimize=True)
        / len(coords)
    )
    second = (
        np.einsum("dn,nac->dac", weights, coords * coords, optimize=True)
        / len(coords)
    )
    atom = np.sqrt(
        np.sum(
            np.maximum(second - mean * mean, 0.0),
            axis=2,
        )
    )
    boot = (
        np.stack(
            [np.bincount(residx, weights=row, minlength=len(labels)) for row in atom]
        ) / np.maximum(np.bincount(residx, minlength=len(labels)), 1)
    )
    return {
        "positions": np.arange(len(labels), dtype=float),
        "labels": labels,
        "point": point,
        "draws": boot,
        "temporal": meta
    }


def collect_temporal_profiles(repout: Path, cfg: dict) -> dict[str, dict]:
    out = temporal_rmsf_profiles(repout, cfg)
    out.update(temporal_density_profiles(repout, cfg))
    out.update(temporal_structure_profiles(repout, cfg))
    out.update(temporal_residue_membrane_profiles(repout, cfg))
    out.update(temporal_membrane_order_profiles(repout, cfg))
    return out


def _profile_label_map(item: dict) -> dict[str, int]:
    labels = [str(label) for label in item["labels"]]
    if len(labels) != len(set(labels)):
        raise ValueError(
            "Profile labels must be unique within each profile for "
            "cross-replicate alignment"
        )
    return {label: i for i, label in enumerate(labels)}


def _common_profile_axis(items: list[dict]) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Align profile columns by label instead of silently truncating them."""
    if not items:
        return np.array([], float), np.array([], str), []
    maps = [_profile_label_map(item) for item in items]
    reference_labels = [str(label) for label in items[0]["labels"]]
    common = [
        label
        for label in reference_labels
        if all(label in mapping for mapping in maps[1:])
    ]
    if not common:
        return np.array([], float), np.array([], str), []
    reference_map = maps[0]
    positions = np.asarray(
        [items[0]["positions"][reference_map[label]] for label in common],
        float,
    )
    indices = [
        np.asarray([mapping[label] for label in common], dtype=int)
        for mapping in maps
    ]
    return positions, np.asarray(common, dtype=str), indices


def analyze_profiles(
    cfg: dict,
    root: Path,
    system_repouts: dict[str, list[Path]],
) -> None:
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    draws = int(stats.get("bootstrap_draws", 2000))
    level = float(stats.get("ci", 0.95))
    rng = np.random.default_rng(int(stats.get("seed", 2026)) + 1)
    resample = stats.get("hierarchy", {}).get("resample_replicates", True)
    outdir = ensure_dir(root / "project_summary" / "statistics" / "profiles")
    sources = {
        (system, repout.name): collect_temporal_profiles(repout, cfg)
        for system, reps in system_repouts.items()
        for repout in reps
    }
    profiles = (
        sorted(set().union(*(set(value) for value in sources.values())))
        if sources
        else []
    )
    rep_points = []
    rep_summary = []
    sys_points = []
    sys_summary = []
    comps = []

    for profile in profiles:
        prepared = {}
        for system, reps in system_repouts.items():
            for repout in reps:
                source = sources[(system, repout.name)]
                if profile not in source:
                    continue
                item = _prepare_profile(source[profile], cfg, draws, rng)
                prepared[(system, repout.name)] = item
                rep_points.extend(
                    {
                        "system": system,
                        "replicate": repout.name,
                        "profile": profile,
                        "position": float(position),
                        "position_label": str(label),
                        "value": float(value),
                    }
                    for position, label, value in zip(
                        item["positions"],
                        item["labels"],
                        item["point"],
                    )
                )
                lo, hi = _ci(np.nanmean(item["draws"], axis=1), level)
                rep_summary.append({
                    "system": system,
                    "replicate": repout.name,
                    "profile": profile,
                    **_describe(item["point"]),
                    "mean_profile_ci_low": lo,
                    "mean_profile_ci_high": hi,
                    **{
                        key: value
                        for key, value in item["temporal"].items()
                        if key != "acf"
                    },
                    "bootstrap_draws": draws,
                    "ci_level": level,
                })

        system_profiles = {}
        for system, reps in system_repouts.items():
            rep_items = [
                (repout.name, prepared[(system, repout.name)])
                for repout in reps
                if (system, repout.name) in prepared
            ]
            if not rep_items:
                continue
            items = [item for _, item in rep_items]
            positions, labels, column_indices = _common_profile_axis(items)
            if not len(labels):
                continue
            aligned_draws = [
                item["draws"][:, indices]
                for item, indices in zip(items, column_indices)
            ]
            aligned_points = [
                item["point"][indices]
                for item, indices in zip(items, column_indices)
            ]
            boot = _hierarchical(
                aligned_draws,
                draws,
                rng,
                "mean",
                resample,
            )
            point = np.nanmean(np.stack(aligned_points), axis=0)
            alpha = (1 - level) / 2
            low, high = np.nanquantile(boot, [alpha, 1 - alpha], axis=0)
            for index, (position, label) in enumerate(zip(positions, labels)):
                sys_points.append({
                    "system": system,
                    "profile": profile,
                    "position": float(position),
                    "position_label": str(label),
                    "n_replicates": len(items),
                    "mean": float(point[index]),
                    "median": float(
                        np.nanmedian([values[index] for values in aligned_points])
                    ),
                    "mean_ci_low": float(low[index]),
                    "mean_ci_high": float(high[index]),
                    "bootstrap_draws": draws,
                    "ci_level": level,
                })
            glo, ghi = _ci(np.nanmean(boot, axis=1), level)
            temporal_rows = [item["temporal"] for item in items]
            sys_summary.append({
                "system": system,
                "profile": profile,
                "n_replicates": len(items),
                "mean_profile_value": float(np.nanmean(point)),
                "mean_profile_ci_low": glo,
                "mean_profile_ci_high": ghi,
                **_range_meta(temporal_rows),
                "bootstrap_draws": draws,
                "ci_level": level,
            })
            system_profiles[system] = {
                "draws": boot,
                "positions": positions,
                "labels": labels,
                "rep_draws": {
                    rep_name: item["draws"][:, indices]
                    for (rep_name, item), indices in zip(
                        rep_items,
                        column_indices,
                    )
                },
            }

        available = [system for system in system_repouts if system in system_profiles]
        for system_1, system_2 in combinations(available, 2):
            sys_1 = system_profiles[system_1]
            sys_2 = system_profiles[system_2]
            map_1 = {str(label): i for i, label in enumerate(sys_1["labels"])}
            map_2 = {str(label): i for i, label in enumerate(sys_2["labels"])}
            common_labels = [
                str(label)
                for label in sys_1["labels"]
                if str(label) in map_2
            ]
            if not common_labels:
                continue
            idx_1 = np.asarray([map_1[label] for label in common_labels], int)
            idx_2 = np.asarray([map_2[label] for label in common_labels], int)
            paired = (
                stats.get("hypothesis_tests", {}).get("design", "independent")
                == "paired"
            )
            if paired:
                rep_draws_1 = {
                    rep: values[:, idx_1]
                    for rep, values in sys_1["rep_draws"].items()
                }
                rep_draws_2 = {
                    rep: values[:, idx_2]
                    for rep, values in sys_2["rep_draws"].items()
                }
                delta = _hierarchical_difference(
                    rep_draws_1,
                    rep_draws_2,
                    draws,
                    rng,
                    "mean",
                    resample,
                )
            else:
                delta = (
                    sys_2["draws"][:, idx_2]
                    - sys_1["draws"][:, idx_1]
                )
            if not len(delta):
                continue
            alpha = (1 - level) / 2
            low, high = np.nanquantile(delta, [alpha, 1 - alpha], axis=0)

            for index, label in enumerate(common_labels):
                position = float(sys_1["positions"][idx_1[index]])
                amap = {}
                bmap = {}
                for (system, replicate), item in prepared.items():
                    label_map = _profile_label_map(item)
                    if label not in label_map:
                        continue
                    value = float(item["point"][label_map[label]])
                    if not np.isfinite(value):
                        continue
                    if system == system_1:
                        amap[replicate] = value
                    elif system == system_2:
                        bmap[replicate] = value
                pgt, plt = _direction(delta[:, index])
                comps.append({
                    "system_1": system_1,
                    "system_2": system_2,
                    "profile": profile,
                    "position": position,
                    "position_label": label,
                    "difference_system2_minus_system1": float(
                        np.nanmean(delta[:, index])
                    ),
                    "difference_ci_low": float(low[index]),
                    "difference_ci_high": float(high[index]),
                    "probability_difference_gt_0": pgt,
                    "probability_difference_lt_0": plt,
                    "p_value": _permutation(amap, bmap, cfg, rng),
                    "bootstrap_draws": draws,
                    "ci_level": level,
                })

    outputs = [
        ("replicate_profiles.parquet", rep_points),
        ("replicate_profile_statistics.parquet", rep_summary),
        ("system_profiles.parquet", sys_points),
        ("system_profile_statistics.parquet", sys_summary),
    ]
    for name, rows in outputs:
        if rows:
            atomic_parquet(outdir / name, pd.DataFrame(rows))

    comp = pd.DataFrame(comps)
    if not comp.empty:
        if (
            stats.get("hypothesis_tests", {}).get(
                "multiple_testing",
                "fdr_bh",
            )
            == "fdr_bh"
        ):
            comp = _fdr(
                comp,
                "p_value",
                "q_value",
                ["system_1", "system_2", "profile"],
            )
        atomic_parquet(outdir / "system_comparisons.parquet", comp)

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
    for kind in sorted(kinds):
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
            prepared[key] = (times, presence, life)
        for interaction in sorted(interactions):
            repdata = {}
            for system, reps in system_repouts.items():
                for repout in reps:
                    times, presence, life = prepared[(system, repout.name)]
                    if not len(times):
                        continue
                    pt = presence.loc[
                        presence.interaction == interaction,
                        "time_ps"
                    ].to_numpy(float) if not presence.empty else np.array([], float)
                    ls = (
                        life[life.interaction == interaction] if not life.empty else pd.DataFrame()
                    )
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


def analyze_project(cfg: dict, root: Path, system_repouts: dict[str, list[Path]]) -> None:
    analyze_time_series(cfg, root, system_repouts)
    analyze_profiles(cfg, root, system_repouts)
    analyze_interactions(cfg, root, system_repouts)
