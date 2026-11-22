"""Shared, replicate-balanced state model with an optional SDC4 preset.

One global feature scaler and one k-means model are fit across *all* systems
and replicates. A single state label has the same meaning in every system.
Only saved scalar time-series are read: no second trajectory traversal.

Sampling diagnostics report replicate-to-replicate variation and first/second
half JS divergence; frame-wise observations are never treated as independent
replicates for system-level confidence intervals.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from ..core.persistence import atomic_json, atomic_parquet, ensure_dir


def _lloyd(data, k, seed, restarts=8, max_iter=200):
    """Deterministic k-means++ with multiple restarts; no sklearn dependency."""
    if len(data) < k or len(np.unique(data, axis=0)) < k:
        raise ValueError(f"Need at least {k} distinct conformations for {k} states")
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(restarts):
        centers = [data[int(rng.integers(len(data)))]]
        for _ in range(1, k):
            d2 = np.min(np.stack([
                np.sum((data - center)**2, axis=1) for center in centers
            ]), axis=0)
            total = d2.sum()
            index = int(rng.choice(len(data), p=d2/total)) if total > 0 else 0
            centers.append(data[index])
        centers = np.asarray(centers, float)
        for _ in range(max_iter):
            distances = np.sum((data[:, None] - centers[None])**2, axis=2)
            labels = np.argmin(distances, axis=1)
            updated = centers.copy()
            for j in range(k):
                subset = data[labels == j]
                if len(subset):
                    updated[j] = subset.mean(axis=0)
                else:
                    nearest = distances.min(axis=1)
                    updated[j] = data[int(np.argmax(nearest))]
            if np.allclose(centers, updated, atol=1e-8, rtol=0):
                centers = updated
                break
            centers = updated
        inertia = float(np.sum(np.min(np.sum(
            (data[:, None] - centers[None])**2, axis=2
        ), axis=1)))
        if best is None or inertia < best[0]:
            best = inertia, centers
    return best[1], best[0]


def _assign(data, centers):
    return np.argmin(np.sum((data[:, None] - centers[None])**2, axis=2), axis=1)


def _moving_block_ci(states, k, draws, block, rng):
    """Temporal block bootstrap within a replicate; descriptive uncertainty."""
    n = len(states)
    blocks = max(1, int(np.ceil(n / block)))
    values = np.empty((draws, k))
    for i in range(draws):
        starts = rng.integers(0, n, size=blocks)
        indices = (starts[:, None] + np.arange(block)) % n
        sample = states[indices.ravel()[:n]]
        values[i] = np.bincount(sample, minlength=k) / n
    return np.quantile(values, [0.025, 0.975], axis=0)


def _write_matched_populations(cfg, out, assignments, system_repouts, n_states):
    from ..statistics.replicate_contrasts import replicate_mean_contrasts
    from ..statistics.windows import matched_grid

    grid = matched_grid(cfg, system_repouts)
    if grid is None:
        return
    full = (assignments.groupby(["system", "replicate", "state"]).size()
            .unstack("state", fill_value=0).reindex(columns=range(n_states), fill_value=0))
    full = full.div(full.sum(axis=1), axis=0).rename_axis("state", axis=1)
    full = full.stack().rename("value").reset_index()
    full_contrasts = replicate_mean_contrasts(full, ["state"], system_repouts, cfg["statistics"])
    if not full_contrasts.empty:
        atomic_parquet(out / "full_population_comparisons.parquet", full_contrasts)
    rows = []
    for (system, replicate), part in assignments.groupby(["system", "replicate"], sort=True):
        subset = part.loc[part.time_ps.between(grid[0], grid[-1])].sort_values("time_ps")
        complete = np.array_equal(subset.time_ps.to_numpy(float), grid)
        for state in range(n_states):
            count = int(subset.state.eq(state).sum())
            rows.append(dict(system=system, replicate=replicate, state=state,
                             frames=count, expected_frames=len(grid), observed_frames=len(subset),
                             window_start_ps=grid[0], window_end_ps=grid[-1],
                             value=count / len(grid) if complete else np.nan))
    values = pd.DataFrame(rows)
    if values.empty:
        return
    contrasts = replicate_mean_contrasts(values, ["state"], system_repouts, cfg["statistics"])
    atomic_parquet(out / "matched_replicate_populations.parquet", values)
    if not contrasts.empty:
        atomic_parquet(out / "matched_population_comparisons.parquet", contrasts)


def _fit_shared_states(cfg, root, system_repouts):
    spec = cfg.get("analyses", {}).get("shared_states", {})
    if not spec.get("enabled", False):
        return None
    from ..contracts.state_feature_spec import SDC4_FEATURES, validate_state_features
    from .state_features import load_state_features
    validate_state_features(spec)
    sdc4 = spec.get("preset", "sdc4") == "sdc4"
    if sdc4 and not cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False):
        raise ValueError("shared_states preset sdc4 requires analyses.sdc4_mechanisms.enabled")
    features = SDC4_FEATURES if spec.get("features") is None else spec["features"]
    k = int(spec.get("n_states", 3))
    seed = int(spec.get("seed", 2026))
    draws = int(spec.get("bootstrap_draws", 300))
    train_n = int(spec.get("training_frames_per_replicate", 2000))
    if k < 2 or draws < 20 or train_n < k:
        raise ValueError("shared_states requires n_states>=2, bootstrap_draws>=20 and training_frames_per_replicate>=n_states")
    rng = np.random.default_rng(seed)
    replicates = {}
    repout_paths = {}
    observation_indices = {}
    observation_counts = {}
    regular_sampling = {}
    training = []
    for system, paths in system_repouts.items():
        for repout in paths:
            path = Path(repout) / "extensions/sdc4_mechanisms/timeseries.parquet"
            if sdc4:
                if not path.is_file():
                    raise FileNotFoundError(path)
                data = pd.read_parquet(path)
            else:
                data = load_state_features(repout, features, spec["feature_sources"])
                path = Path(repout)
            missing = [name for name in features if name not in data]
            if missing:
                raise ValueError(f"Shared-state features absent in {path}: {missing}")
            if not np.isfinite(data["time_ps"].to_numpy(float)).all():
                raise ValueError(f"Nonfinite frame times in {path}")
            if data["time_ps"].duplicated().any():
                raise ValueError(f"Duplicate frame times in {path}")
            data = data.sort_values("time_ps").reset_index(drop=True)
            X = data[features].to_numpy(float)
            good = np.all(np.isfinite(X), axis=1)
            if not np.any(good):
                raise ValueError(f"No complete feature vectors in {path}")
            if good.mean() < float(spec.get("min_complete_fraction", 0.90)):
                raise ValueError(f"Too many missing state-model feature vectors in {path}: {good.mean():.1%} complete")
            clean = data.loc[good].reset_index(drop=True)
            key = (system, Path(repout).name)
            replicates[key] = clean
            observation_indices[key] = np.flatnonzero(good)
            observation_counts[key] = len(data)
            intervals = np.diff(data["time_ps"].to_numpy(float))
            regular_sampling[key] = len(intervals) > 0 and np.allclose(
                intervals, intervals[0], rtol=1e-5, atol=1e-8)
            repout_paths[key] = Path(repout)
            x = X[good]
            # Evenly spaced, equal-count training points per replicate.
            indices = np.linspace(0, len(x)-1, min(len(x), train_n)).astype(int)
            training.append(x[indices])
    # Replicate-balanced scaling: the same number of training frames from
    # each replicate; do not let long trajectories dominate state definitions.
    if not training:
        raise ValueError("Shared states require at least one saved replicate")
    n_train = min(map(len, training))
    training = np.concatenate([
        sample[np.linspace(0, len(sample)-1, n_train).astype(int)]
        for sample in training
    ])
    location = np.median(training, axis=0)
    scale = np.std(training, axis=0)
    if np.any(scale < 1e-8):
        constant = [features[i] for i in np.flatnonzero(scale < 1e-8)]
        raise ValueError(f"State features are constant across training data: {constant}")
    Z = (training - location) / scale
    centers, inertia = _lloyd(Z, k, seed)
    # Sort labels by first feature's centroid for deterministic interpretable
    # state indices (e.g., lower to higher ectodomain height).
    centers = centers[np.argsort(centers[:, 0])]
    # PCA is diagnostic/visualization only; clustering uses all features.
    _, _, vt = np.linalg.svd(Z - Z.mean(axis=0), full_matrices=False)
    pcs = vt[:min(2, len(features))]
    out = ensure_dir(Path(root) / "project_summary" / "shared_states")
    centroids = pd.DataFrame(centers * scale + location, columns=features)
    centroids.insert(0, "state", np.arange(k))
    atomic_parquet(out / "state_centroids.parquet", centroids)
    assignments, populations, diagnostics = [], [], []
    rep_fraction = defaultdict(list)
    outcomes = None
    if sdc4:
        from ..mechanisms.state_outcomes import SDC4StateOutcomes
        outcomes = SDC4StateOutcomes()
    for (system, replicate), df in replicates.items():
        x = df[features].to_numpy(float)
        labels = _assign((x - location) / scale, centers)
        times = df["time_ps"].to_numpy(float)
        projections = ((x - location) / scale) @ pcs.T
        table = pd.DataFrame({
            "system": system, "replicate": replicate,
            "time_ps": times, "state": labels,
            "pc1": projections[:, 0],
            "pc2": projections[:, 1] if projections.shape[1] > 1 else np.nan,
        })
        assignments.append(table)
        n = len(labels)
        counts = np.bincount(labels, minlength=k)
        fractions = counts / n
        block = max(2, int(spec.get("block_frames") or np.sqrt(n)))
        complete = observation_counts[(system, replicate)] == n
        temporal_valid = complete and regular_sampling[(system, replicate)]
        if temporal_valid:
            lo, hi = _moving_block_ci(labels, k, draws, min(block, n), rng)
        else:
            # Dropping feature frames or unequal time intervals makes adjacent
            # retained rows unsuitable for frame-based temporal resampling.
            lo = hi = np.full(k, np.nan)
        rep_fraction[system].append(fractions)
        for j in range(k):
            populations.append({
                "system": system, "replicate": replicate, "state": j,
                "frames": int(counts[j]), "total_frames": n,
                "fraction": float(fractions[j]),
                "within_replicate_block_ci_low": float(lo[j]),
                "within_replicate_block_ci_high": float(hi[j]),
            })
        half = n // 2
        if half > 0 and n-half > 0:
            a = np.bincount(labels[:half], minlength=k) / half
            b = np.bincount(labels[half:], minlength=k) / (n-half)
            js = float(jensenshannon(a, b, base=2.)**2)
        else:
            js = np.nan
        adjacent = np.diff(observation_indices[(system, replicate)]) == 1
        changed = np.diff(labels) != 0
        transitions = int(np.count_nonzero(changed & adjacent))
        diagnostics.append({
            "system": system, "replicate": replicate, "frames": n,
            "input_frames": observation_counts[(system, replicate)],
            "missing_feature_frames": observation_counts[(system, replicate)] - n,
            "complete_feature_fraction": n / observation_counts[(system, replicate)],
            "transitions_across_missing_gaps": int(np.count_nonzero(changed & ~adjacent)),
            "temporal_resampling_status": "available" if temporal_valid else "unavailable_missing_features_or_irregular_times",
            "first_second_half_js_divergence": js,
            "observed_transitions": transitions,
            "block_frames": block,
            "duration_ps": float(times[-1] - times[0]) if n > 1 else 0.,
        })
        if not sdc4:
            continue
        outcomes.add(df, system=system, replicate=replicate,
                     repout=repout_paths[(system, replicate)],
                     labels=labels, times=times, n_states=k, temporal_valid=temporal_valid)
    if outcomes is not None:
        outcomes.write(out, rng, draws)
    frame_assignments = pd.concat(assignments, ignore_index=True)
    atomic_parquet(out / "frame_assignments.parquet", frame_assignments)
    _write_matched_populations(cfg, out, frame_assignments, system_repouts, k)
    atomic_parquet(out / "replicate_populations.parquet", pd.DataFrame(populations))
    atomic_parquet(out / "sampling_diagnostics.parquet", pd.DataFrame(diagnostics))
    systems = []
    for system, vectors in rep_fraction.items():
        matrix = np.asarray(vectors)
        nreps = len(matrix)
        for j in range(k):
            vals = matrix[:, j]
            if nreps >= 2:
                sample = rng.integers(0, nreps, size=(draws, nreps))
                means = vals[sample].mean(axis=1)
                ci = np.quantile(means, [.025, .975])
            else:
                ci = [np.nan, np.nan]
            systems.append({
                "system": system, "state": j, "n_replicates": nreps,
                "mean_replicate_fraction": float(vals.mean()),
                "replicate_sd": float(vals.std(ddof=1)) if nreps > 1 else np.nan,
                "replicate_bootstrap_ci_low": float(ci[0]),
                "replicate_bootstrap_ci_high": float(ci[1]),
            })
    atomic_parquet(out / "system_populations.parquet", pd.DataFrame(systems))
    atomic_json(out / "model.json", {
        "features": features, "n_states": k, "seed": seed,
        "preset": "sdc4" if sdc4 else None,
        "feature_sources": spec.get("feature_sources", {}), "alignment": "exact",
        "scaling": "standardize",
        "replicate_time_windows": [
            {"system": system, "replicate": replica, "first_time_ps": float(data.time_ps.iloc[0]),
             "last_time_ps": float(data.time_ps.iloc[-1]), "input_frames": observation_counts[(system, replica)],
             "complete_frames": len(data), "replica_source": str(repout_paths[(system, replica)])}
            for (system, replica), data in replicates.items()],
        "scaler_median": location.tolist(), "scaler_sd": scale.tolist(),
        "centroids_standardized": centers.tolist(),
        "pca_directions": pcs.tolist(), "training_inertia": inertia,
        "training_frames_per_replicate": n_train,
        "n_replicates": len(replicates),
        "missing_frame_policy": "Populations condition on complete feature frames; transitions exclude missing gaps. Temporal block intervals and shift nulls are omitted for incomplete or irregularly sampled traces.",
        "interpretation": "Geometric k-means groups, not kinetic/metastable or causal states. One global model. System CIs resample independent replicates; "
                          "within-replicate block CIs and half-trajectory JS are descriptive "
                          "and cannot establish equilibrium or causal effects.",
    })
    return out


def fit_shared_states(cfg, root, system_repouts):
    """Publish a complete model generation; preserve old results on failure."""
    import tempfile
    if not cfg.get("analyses", {}).get("shared_states", {}).get("enabled", False):
        return None
    root = ensure_dir(root)
    with tempfile.TemporaryDirectory(prefix=".shared-states-", dir=root) as stage:
        generated = _fit_shared_states(cfg, stage, system_repouts)
        return publish_shared_states(generated, root)


def publish_shared_states(generated, root):
    """Publish a prepared generation; keep prior scientific outputs archived."""
    from ..core.persistence import archive_output
    if generated is None:
        return None
    destination = ensure_dir(Path(root) / "project_summary") / "shared_states"
    previous = None
    if destination.exists():
        previous = archive_output(destination, destination.parent / "archive", "shared_states")
    try:
        generated.rename(destination)
    except OSError:
        if previous is not None:
            previous.rename(destination)
        raise
    return destination
