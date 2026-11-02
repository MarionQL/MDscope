"""Shared, replicate-balanced conformational state model for SDC4 systems.

One global feature scaler and one k-means model are fit across *all* systems
and replicates. A single state label has the same meaning in every system.
Only saved mechanism time-series are read: no second trajectory traversal.

Sampling diagnostics report replicate-to-replicate variation and first/second
half JS divergence; frame-wise observations are never treated as independent
replicates for system-level confidence intervals.
"""
from __future__ import annotations

from pathlib import Path
from collections import defaultdict
import json

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon

from ..config import atomic_json, atomic_parquet, ensure_dir


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


def fit_shared_states(cfg, root, system_repouts):
    spec = cfg.get("analyses", {}).get("shared_states", {})
    if not spec.get("enabled", False):
        return None
    if not cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False):
        raise ValueError("shared_states requires analyses.sdc4_mechanisms.enabled")
    features = spec.get("features") or [
        "sdc4_ectodomain_A_height_A", "sdc4_ectodomain_B_height_A",
        "sdc4_ectodomain_A_tilt_deg", "sdc4_ectodomain_B_tilt_deg",
        "sdc4_tm_crossing_angle_deg", "sdc4_local_thickness_A",
    ]
    if not isinstance(features, list) or not features or len(set(features)) != len(features):
        raise ValueError("shared_states.features must be a nonempty list of distinct columns")
    k = int(spec.get("n_states", 3))
    seed = int(spec.get("seed", 2026))
    draws = int(spec.get("bootstrap_draws", 300))
    train_n = int(spec.get("training_frames_per_replicate", 2000))
    if k < 2 or draws < 20 or train_n < k:
        raise ValueError("shared_states requires n_states>=2, bootstrap_draws>=20 and training_frames_per_replicate>=n_states")
    rng = np.random.default_rng(seed)
    replicates = {}
    repout_paths = {}
    training = []
    for system, paths in system_repouts.items():
        for repout in paths:
            path = Path(repout) / "extensions/sdc4_mechanisms/timeseries.parquet"
            if not path.is_file():
                raise FileNotFoundError(path)
            data = pd.read_parquet(path)
            missing = [name for name in features if name not in data]
            if missing:
                raise ValueError(f"Shared-state features absent in {path}: {missing}")
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
            repout_paths[key] = Path(repout)
            x = X[good]
            # Evenly spaced, equal-count training points per replicate.
            indices = np.linspace(0, len(x)-1, min(len(x), train_n)).astype(int)
            training.append(x[indices])
    # Replicate-balanced scaling: the same number of training frames from
    # each replicate; do not let long trajectories dominate state definitions.
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
    assignments, populations, diagnostics, conditioned = [], [], [], []
    rep_fraction = defaultdict(list)
    lipid_state_rows = []
    nxip_conditioned = []
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
        lo, hi = _moving_block_ci(labels, k, draws, min(block, n), rng)
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
        transitions = int(np.count_nonzero(np.diff(labels))) if n > 1 else 0
        diagnostics.append({
            "system": system, "replicate": replicate, "frames": n,
            "first_second_half_js_divergence": js,
            "observed_transitions": transitions,
            "block_frames": block,
            "duration_ps": float(times[-1] - times[0]) if n > 1 else 0.,
        })
        # Conditional NxIP outcomes: define states exclusively from the
        # configured geometric features; do not include NxIP as a predictor.
        # This table is descriptive within each independent replicate.
        for chain in ("A", "B"):
            prefix = f"sdc4_motif_NxIP_PRO{chain}_"
            endpoints = [
                col for col in df.columns if col.startswith(prefix)
                and (col.endswith("_steric_occlusion_fraction")
                     or col.endswith("_sasa_nm2")
                     or col.endswith("_glycan_contacts")
                     or col.endswith("_protein_neighbor_residues"))
            ]
            for endpoint in endpoints:
                values = df[endpoint].to_numpy(float)
                for state in range(k):
                    mask = (labels == state) & np.isfinite(values)
                    other = (labels != state) & np.isfinite(values)
                    nxip_conditioned.append({
                        "system": system, "replicate": replicate,
                        "chain": chain, "state": state, "endpoint": endpoint,
                        "frames_in_state": int(mask.sum()),
                        "frames_outside_state": int(other.sum()),
                        "mean_in_state": float(values[mask].mean()) if mask.any() else np.nan,
                        "mean_outside_state": float(values[other].mean()) if other.any() else np.nan,
                        "within_replicate_difference": (
                            float(values[mask].mean() - values[other].mean())
                            if mask.any() and other.any() else np.nan
                        ),
                    })
        # A frame-aligned downstream outcome; lipid contacts do not enter
        # geometric state training or assignments.
        lipid_path = repout_paths[(system, replicate)] / "extensions/sdc4_mechanisms/lipid_region_specificity.parquet"
        if lipid_path.is_file():
            lipid_data = pd.read_parquet(lipid_path)
            if not lipid_data.empty:
                merged = lipid_data.merge(
                    pd.DataFrame({"time_ps": times, "state": labels}),
                    on="time_ps", how="inner", validate="many_to_one",
                )
                for keys, subgroup in merged.groupby(
                    ["region", "species", "leaflet", "moiety", "state"],
                    sort=True,
                ):
                    valid = subgroup[subgroup["available_molecules"] > 0]
                    lipid_state_rows.append({
                        "system": system, "replicate": replicate,
                        "region": keys[0], "species": keys[1],
                        "leaflet": keys[2], "moiety": keys[3],
                        "state": int(keys[4]),
                        "frames_in_state": int(len(subgroup)),
                        "frames_with_species": int(len(valid)),
                        "mean_contact_fraction": (
                            float(valid["contact_fraction"].mean())
                            if len(valid) else np.nan
                        ),
                    })
        hs_names = sorted({
            name[len("sdc4_hs_"):-len("_hydration_water_per_heavy_atom")]
            for name in df.columns if name.startswith("sdc4_hs_")
            and name.endswith("_hydration_water_per_heavy_atom")
        })
        for hs in hs_names:
            prefix = f"sdc4_hs_{hs}_"
            if prefix + "membrane_contact_state" not in df:
                continue
            for state in range(k):
                for contact in (0, 1):
                    mask = (labels == state) & (
                        df[prefix + "membrane_contact_state"].to_numpy(int) == contact
                    )
                    if not np.any(mask):
                        continue
                    row = {
                        "system": system, "replicate": replicate, "hs_chain": hs,
                        "state": state, "membrane_contact_state": contact,
                        "frames": int(mask.sum()),
                        "hydration_water_per_heavy_atom_mean": float(
                            df.loc[mask, prefix + "hydration_water_per_heavy_atom"].mean()
                        ),
                    }
                    for col in df.columns:
                        if col.startswith(prefix + "ion_") and col.endswith("_count"):
                            row[col[len(prefix):] + "_mean"] = float(df.loc[mask, col].mean())
                    conditioned.append(row)
    nxip_df = pd.DataFrame(nxip_conditioned, columns=[
        "system", "replicate", "chain", "state", "endpoint",
        "frames_in_state", "frames_outside_state", "mean_in_state",
        "mean_outside_state", "within_replicate_difference",
    ])
    atomic_parquet(out / "state_conditioned_nxip.parquet", nxip_df)
    # A condition contrast is computed only from independent replicate-level
    # contrasts. With three trajectories per condition these intervals are
    # exploratory; no frame-based p-values or significance stars are emitted.
    comparisons = []
    if not nxip_df.empty:
        systems_for_comparison = list(dict.fromkeys(nxip_df["system"]))
        if len(systems_for_comparison) == 2:
            for (chain, state, endpoint), group in nxip_df.groupby(
                ["chain", "state", "endpoint"], sort=True
            ):
                arrays = [
                    group.loc[group["system"] == system,
                              "within_replicate_difference"].dropna().to_numpy(float)
                    for system in systems_for_comparison
                ]
                if min(map(len, arrays)) < 2:
                    continue
                difference = float(np.mean(arrays[1]) - np.mean(arrays[0]))
                indices0 = rng.integers(0, len(arrays[0]), (draws, len(arrays[0])))
                indices1 = rng.integers(0, len(arrays[1]), (draws, len(arrays[1])))
                boot = arrays[1][indices1].mean(axis=1) - arrays[0][indices0].mean(axis=1)
                lo, hi = np.quantile(boot, [.025, .975])
                comparisons.append({
                    "chain": chain, "state": int(state), "endpoint": endpoint,
                    "system_0": systems_for_comparison[0],
                    "system_1": systems_for_comparison[1],
                    "n_replicates_0": len(arrays[0]),
                    "n_replicates_1": len(arrays[1]),
                    "difference_system1_minus_system0": difference,
                    "replicate_bootstrap_ci_low": float(lo),
                    "replicate_bootstrap_ci_high": float(hi),
                    "inference": "exploratory; no multiplicity-adjusted test",
                })
    atomic_parquet(out / "state_conditioned_nxip_comparisons.parquet",
                   pd.DataFrame(comparisons, columns=[
                       "chain", "state", "endpoint", "system_0", "system_1",
                       "n_replicates_0", "n_replicates_1",
                       "difference_system1_minus_system0",
                       "replicate_bootstrap_ci_low", "replicate_bootstrap_ci_high",
                       "inference",
                   ]))
    atomic_parquet(
        out / "state_conditioned_lipid_specificity.parquet",
        pd.DataFrame(lipid_state_rows, columns=[
            "system", "replicate", "region", "species", "leaflet", "moiety",
            "state", "frames_in_state", "frames_with_species",
            "mean_contact_fraction",
        ]),
    )
    atomic_parquet(out / "frame_assignments.parquet", pd.concat(assignments, ignore_index=True))
    atomic_parquet(out / "replicate_populations.parquet", pd.DataFrame(populations))
    atomic_parquet(out / "sampling_diagnostics.parquet", pd.DataFrame(diagnostics))
    atomic_parquet(out / "state_conditioned_hs_solvation.parquet", pd.DataFrame(conditioned))
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
        "scaler_median": location.tolist(), "scaler_sd": scale.tolist(),
        "centroids_standardized": centers.tolist(),
        "pca_directions": pcs.tolist(), "training_inertia": inertia,
        "training_frames_per_replicate": n_train,
        "n_replicates": len(replicates),
        "interpretation": "One global model. System CIs resample independent replicates; "
                          "within-replicate block CIs and half-trajectory JS are descriptive "
                          "and cannot establish equilibrium or causal effects.",
    })
    return out

