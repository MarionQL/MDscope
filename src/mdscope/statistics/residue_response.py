"""Exploratory, sequence-constrained cross-metric residue response screen.

Uses only saved replica means and completed profile contrasts. This is a
descriptive grouping of already estimated effects, NOT selective inference or
a confidence interval for a data-discovered region.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from ..statistics_contracts import RESIDUE_RESPONSE_PROFILE

_RESIDUE_LABEL = re.compile(r"^(?P<chain>[^:]+):[A-Za-z]+(?P<resid>-?\d+)$")
_KEYS = ["system_1", "system_2", "profile", "position_label"]


def residue_response_screen(replicates: pd.DataFrame, contrasts: pd.DataFrame,
                            expected_replicates: dict[str, list[str]],
                            *, min_length: int = 3,
                            minimum_cosine: float = 0.7,
                            profiles=None, chains=None):
    """Return (effect matrix, contiguous multi-metric regions).

    Per-profile scale = pooled within-(system,residue) replica standard
    deviation. This preserves the direction of physical-unit differences
    without allowing a metric's choice of units to dominate segmentation.
    Zero-variance and incomplete comparisons are marked unavailable.
    """
    columns = [*_KEYS, "chain", "resid", "difference_system2_minus_system1",
               "difference_ci_low", "difference_ci_high", "pooled_replica_sd",
               "standardized_effect", "ci_excludes_zero", "region_id"]
    region_columns = ["region_id", "system_1", "system_2", "chain",
                      "start_residue", "end_residue", "n_residues",
                      "n_metrics", "n_ci_supported_cells", "residue_labels",
                      "metrics", "exploratory"]
    if replicates.empty or contrasts.empty:
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=region_columns)
    required_rep = {"system", "replicate", "profile", "position_label", "value"}
    required_cmp = set(_KEYS) | {"difference_system2_minus_system1",
                                "difference_ci_low", "difference_ci_high"}
    if not required_rep.issubset(replicates) or not required_cmp.issubset(contrasts):
        raise ValueError("Residue response requires saved replica profiles and bootstrap contrasts")
    if profiles is not None:
        unsupported = [name for name in profiles if not RESIDUE_RESPONSE_PROFILE.match(name)]
        if unsupported:
            raise ValueError(f"Unsupported residue-response profiles: {unsupported}")
        replicates = replicates.loc[replicates.profile.isin(profiles)]
        contrasts = contrasts.loc[contrasts.profile.isin(profiles)]
    r = replicates.loc[replicates.profile.astype(str).str.match(RESIDUE_RESPONSE_PROFILE)].copy()
    c = contrasts.loc[contrasts.profile.astype(str).str.match(RESIDUE_RESPONSE_PROFILE)].copy()
    if r.empty or c.empty:
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=region_columns)
    for frame, identity in ((r, ["system", "replicate", "profile", "position_label"]),
                            (c, _KEYS)):
        if frame.duplicated(identity).any():
            raise ValueError("Ambiguous residue profile identity; repeated labels require topology-index mapping")
    r["value"] = pd.to_numeric(r["value"], errors="coerce")
    c = c.copy()
    # Require every configured independent replica at each protein position.
    support = r.groupby(["system", "profile", "position_label"]).agg(
        n=("value", "size"), finite=("value", lambda x: np.isfinite(x).sum()),
        sd=("value", "std"),
    )
    # Estimate metric scale from within-residue replica variability.
    valid_support = []
    scales = {}
    for (system, profile, label), item in support.iterrows():
        n_expected = len(expected_replicates.get(system, []))
        valid = n_expected >= 2 and item["n"] == n_expected and item["finite"] == n_expected
        valid_support.append((system, profile, label, valid))
    support_ok = pd.DataFrame(valid_support, columns=["system", "profile", "position_label", "ok"])
    for profile, subset in support.reset_index().groupby("profile"):
        finite = subset["sd"].to_numpy(float)
        positive = finite[np.isfinite(finite) & (finite > 0)]
        scales[profile] = float(np.median(positive)) if len(positive) else np.nan
    ok = support_ok.pivot_table(index=["profile", "position_label"], columns="system",
                                values="ok", aggfunc="first")
    lookup = ok.to_dict("index")
    records = []
    for item in c.itertuples(index=False):
        key = (str(item.profile), str(item.position_label))
        support_row = lookup.get(key, {})
        if not all(pd.notna(support_row.get(system)) and
                   bool(support_row[system])
                   for system in (str(item.system_1), str(item.system_2))):
            continue
        found = _RESIDUE_LABEL.fullmatch(str(item.position_label))
        if not found or found["chain"] == "?" or (chains is not None and found["chain"] not in chains):
            continue
        scale = scales.get(item.profile, np.nan)
        difference = float(item.difference_system2_minus_system1)
        low, high = float(item.difference_ci_low), float(item.difference_ci_high)
        standardized = difference / scale if np.isfinite(scale) and scale > 0 else np.nan
        records.append({
            **{name: str(getattr(item, name)) for name in _KEYS},
            "chain": found["chain"], "resid": int(found["resid"]),
            "difference_system2_minus_system1": difference,
            "difference_ci_low": low, "difference_ci_high": high,
            "pooled_replica_sd": scale,
            "standardized_effect": standardized,
            "ci_excludes_zero": bool(np.isfinite(low) and np.isfinite(high) and
                                    (low > 0 or high < 0)),
            "region_id": None,
        })
    result = pd.DataFrame(records, columns=columns)
    if result.empty:
        return result, pd.DataFrame(columns=region_columns)
    if result.duplicated(["system_1", "system_2", "chain", "resid", "profile"]).any():
        raise ValueError("Residue labels are not unique within a chain; check topology identity")
    regions = []
    for (left, right, chain), group in result.groupby(["system_1", "system_2", "chain"], sort=True):
        matrix = group.pivot(index="resid", columns="profile", values="standardized_effect").sort_index()
        # Require two observed measurements on each side of each boundary;
        # do not impute absent features as zero.
        residues = matrix.index.to_numpy(int)
        data = matrix.to_numpy(float)
        start = 0
        def emit(lo, hi):
            if hi - lo < min_length:
                return
            block = matrix.iloc[lo:hi]
            available_metrics = block.columns[block.notna().all(axis=0)]
            if len(available_metrics) < 2:
                return
            identifier = f"{left}__{right}__{chain}__{int(residues[lo])}-{int(residues[hi-1])}"
            mask = ((result.system_1 == left) & (result.system_2 == right) &
                    (result.chain == chain) & result.resid.isin(residues[lo:hi]))
            result.loc[mask, "region_id"] = identifier
            region_rows = result.loc[mask]
            regions.append({
                "region_id": identifier, "system_1": left, "system_2": right,
                "chain": chain, "start_residue": int(residues[lo]),
                "end_residue": int(residues[hi - 1]),
                "n_residues": int(hi - lo), "n_metrics": int(len(available_metrics)),
                "n_ci_supported_cells": int(region_rows.ci_excludes_zero.sum()),
                "residue_labels": ",".join(map(str, residues[lo:hi])),
                "metrics": ",".join(map(str, available_metrics)),
                "exploratory": True,
            })
        for i in range(1, len(residues)):
            shared = np.isfinite(data[i - 1]) & np.isfinite(data[i])
            if shared.sum() < 2 or residues[i] != residues[i-1] + 1:
                split = True
            else:
                a, b = data[i - 1, shared], data[i, shared]
                denom = float(np.linalg.norm(a) * np.linalg.norm(b))
                similarity = float(np.dot(a, b) / denom) if denom > 0 else -1.
                split = similarity < minimum_cosine
            if split:
                emit(start, i)
                start = i
        emit(start, len(residues))
    return result, pd.DataFrame(regions, columns=region_columns)


def write_residue_response(replicates, contrasts, expected_replicates, outdir, cfg):
    """Publish the configured derived screen without touching its profile inputs.

    The completion manifest is written last; readers reject interrupted updates.
    Superseded products are archived, including disabling and empty selections.
    """
    from pathlib import Path

    from ..config import residue_response_settings
    from ..core.persistence import archive_output, atomic_json, atomic_parquet
    settings = residue_response_settings(cfg)
    status = ("disabled" if not settings["enabled"] else
              "empty_selection" if settings["profiles"] == [] or settings["chains"] == []
              else "complete")
    matrix = regions = None
    if status == "complete":
        matrix, regions = residue_response_screen(
            replicates, contrasts, expected_replicates,
            min_length=settings["min_length"], minimum_cosine=settings["minimum_cosine"],
            profiles=settings["profiles"], chains=settings["chains"],
        )
        if matrix.empty:
            status = "no_eligible_data"
    outdir = Path(outdir)
    record = {
        "status": "running", "settings": settings, "exploratory": True,
        "independent_units": "simulation_replicates",
        "replicates": expected_replicates,
        "source_files": ["replicate_profiles.parquet", "system_comparisons.parquet"],
        "input_rows": {"replicate_profiles": len(replicates), "contrasts": len(contrasts)},
    }
    manifest = outdir / "residue_response_manifest.json"
    atomic_json(manifest, record)
    for name in ("residue_response_matrix.parquet", "residue_response_regions.parquet"):
        path = outdir / name
        if path.exists():
            archive_output(path, outdir / ".archives" / "residue_response", name)
    if status == "complete":
        atomic_parquet(outdir / "residue_response_matrix.parquet", matrix)
        atomic_parquet(outdir / "residue_response_regions.parquet", regions)
    record["status"] = status
    atomic_json(manifest, record)
    return record
