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
from ..observables.density_1d import align_density_profiles, temporal_density_profiles
from ..observables.residue_membrane_distance import temporal_residue_membrane_profiles
from ..observables.structure_profiles import temporal_structure_profiles
from ..observables.rmsf import temporal_rmsf_profiles


from .common import _temporal, _block_bootstrap, _block_weights, _hierarchical, _hierarchical_difference, _ci, _describe, _direction, _permutation, _fdr, _range_meta, _regular_times, _point_difference
from .residue_regions import scan_residue_regions


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
    regular = _regular_times(times)
    if spec.get("reducer") != "rmsf":
        values = np.asarray(spec["values"], float)
        if values.ndim != 2 or len(values) != len(times):
            raise ValueError("Profile values must have one row per time sample")
        reducer = str(spec.get("reducer", "mean"))
        complete = np.all(np.isfinite(values), axis=0)
        if not complete.all() and np.isinf(values).any():
            values = np.where(np.isfinite(values), values, np.nan)
        if complete.all():
            # Retain the established reduction order for complete profiles.
            signal = np.sqrt(np.nanmean((values - np.nanmean(values, axis=0)) ** 2, axis=1))
        else:
            signal_values = values[:, complete]
            signal = (np.sqrt(np.mean((signal_values - np.mean(signal_values, axis=0)) ** 2, axis=1))
                      if complete.any() else np.zeros(len(times)))
        meta = _temporal(signal, times, cfg)
        if not regular or not complete.any():
            meta.update(tau_int_frames=np.nan, tau_int_ps=np.nan, effective_n=np.nan)
        boot = (_block_bootstrap(values, meta["block_length_frames"], draws, rng, reducer)
                if regular and complete.any() and draws else np.full((draws, values.shape[1]), np.nan))
        boot[:, ~complete] = np.nan
        return {
            "positions": np.asarray(spec["positions"], float),
            "labels": np.asarray(spec["labels"], str),
                "point": (np.nanmedian if reducer == "median" else np.nanmean)(values, axis=0),
                "draws": boot, "temporal": meta,
                "temporal_valid": complete & regular}
    coords = np.asarray(spec["coords"], np.float32)
    residues = np.asarray(spec["atom_residues"], str)
    signal = np.sqrt(np.mean(np.sum((coords - np.mean(coords, axis=0)) ** 2, axis=2), axis=1))
    meta = _temporal(signal, times, cfg)
    point, labels, residx = _rmsf_profile(
        coords,
        residues,
        spec.get("atom_residue_indices"),
    )
    if not regular or not np.all(np.isfinite(coords)) or not draws:
        if not regular or not np.all(np.isfinite(coords)):
            meta.update(tau_int_frames=np.nan, tau_int_ps=np.nan, effective_n=np.nan)
        return {"positions": np.arange(len(labels), dtype=float), "labels": labels,
                "point": point, "draws": np.full((draws, len(labels)), np.nan),
                "temporal": meta, "temporal_valid": np.full(len(labels), regular and np.all(np.isfinite(coords)))}
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
        "temporal": meta, "temporal_valid": np.ones(len(labels), bool)
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
    cfg = statistics_family_config(cfg, "profiles")
    if cfg is None:
        return
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    stats = cfg["statistics"]
    draws = int(stats.get("bootstrap_draws", 2000))
    system_draws = draws if stats.get("hierarchy", {}).get("enabled", True) else 0
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
        profile_sources = {
            key: source[profile] for key, source in sources.items() if profile in source
        }
        if profile.startswith("density_1d__"):
            profile_sources = align_density_profiles(profile_sources, profile)
        prepared = {}
        for system, reps in system_repouts.items():
            for repout in reps:
                key = (system, repout.name)
                if key not in profile_sources:
                    continue
                item = _prepare_profile(profile_sources[key], cfg, draws, rng)
                prepared[key] = item
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
                lo, hi = _ci(np.mean(item["draws"], axis=1), level)
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
                system_draws,
                rng,
                "mean",
                resample,
            )
            point = np.nanmean(np.stack(aligned_points), axis=0)
            complete = np.all(np.stack([item["temporal_valid"][indices]
                                       for item, indices in zip(items, column_indices)]), axis=0)
            complete &= len(items) == len(reps)
            boot[:, ~complete] = np.nan
            alpha = (1 - level) / 2
            low, high = (np.nanquantile(boot, [alpha, 1 - alpha], axis=0)
                         if system_draws else np.full((2, len(labels)), np.nan))
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
                    "bootstrap_draws": system_draws,
                    "ci_level": level,
                })
            glo, ghi = _ci(np.mean(boot, axis=1), level)
            temporal_rows = [item["temporal"] for item in items]
            sys_summary.append({
                "system": system,
                "profile": profile,
                "n_replicates": len(items),
                "mean_profile_value": float(np.nanmean(point)),
                "mean_profile_ci_low": glo,
                "mean_profile_ci_high": ghi,
                **_range_meta(temporal_rows),
                "bootstrap_draws": system_draws,
                "ci_level": level,
            })
            system_profiles[system] = {
                "draws": boot,
                "positions": positions,
                "labels": labels,
                "temporal_valid": complete,
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
                    system_draws,
                    rng,
                    "mean",
                    resample,
                )
            else:
                delta = (
                    sys_2["draws"][:, idx_2]
                    - sys_1["draws"][:, idx_1]
                )
            complete = sys_1["temporal_valid"][idx_1] & sys_2["temporal_valid"][idx_2]
            if not len(delta):
                delta = np.empty((0, len(common_labels)))
            delta[:, ~complete] = np.nan
            alpha = (1 - level) / 2
            low, high = (np.nanquantile(delta, [alpha, 1 - alpha], axis=0)
                         if len(delta) else np.full((2, len(common_labels)), np.nan))

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
                        np.nanmean(delta[:, index]) if len(delta) and complete[index]
                        else _point_difference(amap, bmap, paired)
                    ),
                    "difference_ci_low": float(low[index]),
                    "difference_ci_high": float(high[index]),
                    "probability_difference_gt_0": pgt,
                    "probability_difference_lt_0": plt,
                    "p_value": _permutation(amap, bmap, cfg, rng),
                    "bootstrap_draws": system_draws,
                    "ci_level": level,
                })

    if rep_points:
        replicate_profiles = pd.DataFrame(rep_points)
        atomic_parquet(outdir / "replicate_profiles.parquet", replicate_profiles)
        replicas = {
            system: [repout.name for repout in repouts]
            for system, repouts in system_repouts.items()
        }
        atomic_parquet(outdir / "residue_regions.parquet",
                       scan_residue_regions(replicate_profiles, replicas))

    outputs = [
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


