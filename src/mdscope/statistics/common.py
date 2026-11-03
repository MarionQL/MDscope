"""
Shared autocorrelation, hierarchical bootstrap, permutation, FDR, and
comparison statistics.
"""

from __future__ import annotations

import math
import os
from concurrent.futures import ThreadPoolExecutor
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


def _bootstrap_workers() -> int:
    """Respect the CPU allocation rather than the physical node size."""
    affinity = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    requested = os.environ.get("MDSCOPE_STATS_WORKERS") or os.environ.get("SLURM_CPUS_PER_TASK")
    return max(1, min(affinity, int(requested))) if requested else max(1, affinity if "SLURM_JOB_ID" in os.environ else 1)


def _reduce_bootstrap_batch(values, indices, reducer):
    reduce = np.nanmedian if reducer == "median" else np.nanmean
    return reduce(values[indices], axis=1)


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
    batch = 32
    workers = min(_bootstrap_workers(), max(1, math.ceil(draws / batch)))
    # Generate indices on the calling thread, in exactly the same order as
    # the serial implementation. Workers only perform deterministic reductions.
    if workers == 1:
        for start in range(0, draws, batch):
            size = min(batch, draws - start)
            out[start:start + size] = _reduce_bootstrap_batch(
                values, _block_indices(len(values), block, size, rng), reducer
            )
        return out

    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = []
        for start in range(0, draws, batch):
            size = min(batch, draws - start)
            indices = _block_indices(len(values), block, size, rng)
            pending.append((start, size, pool.submit(
                _reduce_bootstrap_batch, values, indices, reducer
            )))
            # Bound outstanding batches to avoid holding the entire bootstrap
            # index tensor in memory.
            if len(pending) >= workers:
                first, count, future = pending.pop(0)
                out[first:first + count] = future.result()
        for first, count, future in pending:
            out[first:first + count] = future.result()
    return out


def _block_weights(n: int, block: int, draws: int, rng: np.random.Generator) -> np.ndarray:
    # Interaction statistics need event locations as well as occupancy.
    # Frame multiplicities reproduce a block bootstrap without copying
    # binary trajectories.
    idx = _block_indices(n, block, draws, rng)
    dtype = np.uint16 if int(np.ceil(n / max(1, block))) < 65535 else np.uint32
    weights = np.empty((draws, n), dtype=dtype)
    for i, row in enumerate(idx):
        weights[i] = np.bincount(row, minlength=n)
    return weights


def _hierarchical(
    rep_draws,
    draws: int,
    rng: np.random.Generator,
    reducer="mean",
    resample=True
) -> np.ndarray:
    if draws == 0:
        first = next(iter(rep_draws), np.array([]))
        return np.empty((0,) + np.asarray(first).shape[1:])
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
    if draws == 0:
        return np.empty((0,) + np.asarray(rep_draws_a[labels[0]]).shape[1:])
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


def _point_difference(a, b, paired=False, reducer="mean"):
    """Finite-observation descriptive contrast when temporal draws are unavailable."""
    reduce = np.median if reducer == "median" else np.mean
    if paired:
        differences = finite([b[k] - a[k] for k in sorted(set(a) & set(b))])
        return float(reduce(differences)) if len(differences) else np.nan
    left, right = finite(list(a.values())), finite(list(b.values()))
    return float(reduce(right) - reduce(left)) if len(left) and len(right) else np.nan


def _regular_times(times):
    times = np.asarray(times, float)
    intervals = np.diff(times)
    return bool(times.ndim == 1 and np.all(np.isfinite(times)) and len(intervals)
                and np.all(intervals > 0)
                and np.allclose(intervals, intervals[0], rtol=1e-5, atol=1e-8))


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


