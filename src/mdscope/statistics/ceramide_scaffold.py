"""Descriptive, replica-balanced glycolipid/ceramide scaffold contrasts.

This is deliberately distinct from site-specific glycan chemistry and does
not bootstrap frames or treat molecules as independent condition replicates.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.persistence import atomic_parquet

SOURCE = "extensions/sdc4_mechanisms/ceramide_scaffold_pooled.parquet"
KEYS = ["region", "leaflet", "moiety"]
MEASURES = ["contact_fraction", "matched_local_contact_fraction"]


def summarize_scaffold_replicates(system_repouts):
    records = []
    for system, reps in system_repouts.items():
        for rep in reps:
            rep = Path(rep)
            path = rep / SOURCE
            if not path.is_file():
                raise FileNotFoundError(f"Required ceramide scaffold counts are missing: {path}")
            data = pd.read_parquet(path)
            required = set(["time_ps", "cohort", *KEYS, *MEASURES])
            if not required.issubset(data):
                raise ValueError(f"{path}: missing fields {sorted(required - set(data))}")
            if data.duplicated(["time_ps", "cohort", *KEYS]).any():
                raise ValueError(f"{path}: duplicate ceramide scaffold time/cohort/region")
            for identity, subset in data.groupby(["cohort", *KEYS], sort=True):
                for metric in MEASURES:
                    values = pd.to_numeric(subset[metric], errors="coerce").to_numpy(float)
                    finite = values[np.isfinite(values)]
                    records.append({
                        "system": system, "replicate": rep.name,
                        **dict(zip(["cohort", *KEYS], identity)),
                        "metric": metric, "n_frames": len(values),
                        "n_valid_frames": len(finite),
                        "mean": float(finite.mean()) if len(finite) else np.nan,
                    })
    return pd.DataFrame(records, columns=[
        "system", "replicate", "cohort", *KEYS, "metric",
        "n_frames", "n_valid_frames", "mean",
    ])


def scaffold_condition_comparisons(replicate: pd.DataFrame, system_repouts):
    """Compare cohorts across two configured conditions when unambiguous."""
    columns = [
        "system_1", "cohort_1", "system_2", "cohort_2",
        *KEYS, "metric", "mean_1", "mean_2",
        "difference_system2_minus_system1", "n_replicates_1",
        "n_replicates_2", "status",
    ]
    systems = list(system_repouts)
    if len(systems) != 2 or replicate.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for (*group, metric), subset in replicate.groupby([*KEYS, "metric"], sort=True):
        eligible = {}
        for system in systems:
            present = subset.loc[(subset.system == system) & subset["mean"].notna()]
            cohort_names = sorted(present.cohort.unique())
            eligible[system] = cohort_names
        cohorts = [eligible[system] for system in systems]
        if any(len(names) != 1 for names in cohorts):
            status = "missing_or_ambiguous_cohort"
        elif cohorts[0] == cohorts[1]:
            status = "no_chemical_contrast"
        else:
            status = "ok"
        values = []
        for system, candidate in zip(systems, cohorts):
            cohort = candidate[0] if len(candidate) == 1 else None
            valid = subset.loc[(subset.system == system) & (subset.cohort == cohort)]
            complete = (cohort is not None and valid.replicate.nunique() ==
                        len(system_repouts[system]) and valid["mean"].notna().all())
            values.append((cohort, valid["mean"].to_numpy(float) if complete
                           else np.array([], float)))
        if status == "ok" and any(not len(value) for _, value in values):
            status = "incomplete_replicas"
        rows.append({
            "system_1": systems[0], "cohort_1": values[0][0],
            "system_2": systems[1], "cohort_2": values[1][0],
            **dict(zip([*KEYS, "metric"], (*group, metric))),
            "mean_1": float(values[0][1].mean()) if len(values[0][1]) else np.nan,
            "mean_2": float(values[1][1].mean()) if len(values[1][1]) else np.nan,
            "difference_system2_minus_system1": (
                float(values[1][1].mean() - values[0][1].mean())
                if status == "ok" else np.nan
            ),
            "n_replicates_1": len(values[0][1]),
            "n_replicates_2": len(values[1][1]),
            "status": status,
        })
    return pd.DataFrame(rows, columns=columns)


def write_scaffold_condition_summary(cfg, root, system_repouts):
    if not cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get(
            "ceramide_scaffold_pooling", False):
        return 0
    replicas = summarize_scaffold_replicates(system_repouts)
    contrasts = scaffold_condition_comparisons(replicas, system_repouts)
    directory = Path(root) / "project_summary/statistics/ceramide_scaffold"
    atomic_parquet(directory / "replicate_means.parquet", replicas)
    atomic_parquet(directory / "condition_comparisons.parquet", contrasts)
    return len(contrasts)
