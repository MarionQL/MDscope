"""Deterministic coverage of saved between-system bootstrap intervals."""
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import atomic_parquet


# Identity columns must match the keys used to create the comparison rows.
SOURCES = (
    ("time_series", "time_series", "replicate_statistics", "system_comparisons",
     ("group", "metric"), ("mean", "median")),
    ("matched_time_series", "matched_time_series", "replicate_statistics", "system_comparisons",
     ("group", "metric"), ("mean", "median")),
    ("profiles", "profiles", "replicate_profiles", "system_comparisons",
     ("profile", "position_label"), ("saved",)),
    ("interactions", "interactions", "replicate_interaction_metrics",
     "system_interaction_comparisons", ("interaction_type", "interaction", "metric"), ("saved",)),
    ("coupling", "coupling", "replicate_correlations", "system_comparisons",
     ("mechanism", "x", "y"), ("saved",)),
    ("diffusion", "diffusion", "replicate_common_lag_slopes", "common_lag_comparisons",
     ("metric",), ("saved",)),
    ("solvation_rdfs", "solvation", "replicate_rdf_shells", "rdf_shell_comparisons",
     ("solvent", "target_class", "site_class", "normalization", "shell_cutoff_A"), ("saved",)),
)

COLUMNS = [
    "family", "observable", "system_1", "system_2", "estimator", "status",
    "source_table", "difference_system2_minus_system1", "ci_low", "ci_high", "ci_level",
    "p_value", "q_value",
]


def _identity(frame, keys):
    if "group" in keys:
        frame["group"] = frame["group"].fillna("")
    return frame[list(keys)].apply(
        lambda row: json.dumps(dict(zip(keys, row)), separators=(",", ":"), allow_nan=False),
        axis=1,
    )


def _interval_status(low, high):
    if not (np.isfinite(low) and np.isfinite(high)) or low > high:
        return "unavailable_ci"
    if low > 0:
        return "above_zero"
    if high < 0:
        return "below_zero"
    return "includes_zero"


def write_contrast_coverage(dest: Path, systems, families):
    """Write an exhaustive row audit of the registered statistics families."""
    rows = []
    pairs = list(combinations(systems, 2))
    for family, folder, replica_file, contrast_file, keys, estimators in SOURCES:
        path = dest / folder / f"{contrast_file}.parquet"
        source = str(path.relative_to(dest))
        if family not in families:
            continue
        enabled = families[family] == "enabled"
        if enabled and path.is_file():
            contrasts = pd.read_parquet(path)
            contrasts["observable"] = _identity(contrasts, keys)
        else:
            contrasts = pd.DataFrame(columns=["system_1", "system_2", "observable"])
        replica_path = dest / folder / f"{replica_file}.parquet"
        observed = {}
        lifetime_support = {}
        if enabled and replica_path.is_file():
            replicas = pd.read_parquet(replica_path)
            replicas["observable"] = _identity(replicas, keys)
            observed = replicas.groupby("observable").system.agg(set).to_dict()
            if family == "interactions":
                lifetime = replicas.loc[
                    replicas.metric.isin(["mean_lifetime_ps", "median_lifetime_ps"])
                ]
                lifetime_support = {key: set() for key in lifetime.observable}
                lifetime_support.update(
                    lifetime.loc[np.isfinite(lifetime.value)].groupby("observable").system.agg(set)
                )
        for first, second in pairs:
            base = dict(family=family, system_1=first, system_2=second, source_table=source)
            if not enabled:
                rows.append(dict(base, status="family_disabled"))
                continue
            subset = contrasts.loc[contrasts.system_1.eq(first) & contrasts.system_2.eq(second)]
            reported = set(subset.observable)
            expected = {key for key, present in observed.items() if first in present or second in present}
            for key in sorted(expected - reported):
                status = ("one_system_only" if not {first, second}.issubset(observed[key])
                          else "no_comparison")
                if key in lifetime_support and not {first, second}.issubset(lifetime_support[key]):
                    status = "missing_observed_episode_lifetime"
                for estimator in estimators:
                    rows.append(dict(base, observable=key, estimator=estimator, status=status))
            if not expected and subset.empty:
                rows.append(dict(base, status="no_observables"))
            for row in subset.to_dict("records"):
                for estimator in estimators:
                    prefix = "" if estimator == "saved" else estimator + "_"
                    low = row[prefix + "difference_ci_low"]
                    high = row[prefix + "difference_ci_high"]
                    rows.append(dict(
                        base, observable=row["observable"], estimator=estimator,
                        status=(row["ci_status"] if isinstance(row.get("ci_status"), str)
                                else _interval_status(low, high)),
                        difference_system2_minus_system1=row[prefix + "difference_system2_minus_system1"],
                        ci_low=low, ci_high=high, ci_level=row.get("ci_level"),
                        p_value=row.get(prefix + "p_value"), q_value=row.get(prefix + "q_value"),
                    ))
    coverage = pd.DataFrame(rows, columns=COLUMNS).sort_values(
        ["family", "system_1", "system_2", "observable", "estimator"],
        ignore_index=True,
    )
    atomic_parquet(dest / "coverage" / "contrast_coverage.parquet", coverage)
