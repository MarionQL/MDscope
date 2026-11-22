"""Select descriptive metrics and eligible inferential estimators."""
from __future__ import annotations

from ..contracts.metric_selection import (
    ESTIMATOR_FLAGS, INTERACTION_METRICS, SIMPLE_FAMILY_KEYS,
    _valid_estimator_flags as _valid_estimator_flags, validate_metric_selection as validate_metric_selection,
)

def time_series_selection(cfg):
    """None = legacy all metrics; set() = explicit none."""
    selection = cfg.get("statistics", {}).get("metric_selection")
    if selection is None or "time_series" not in selection:
        return None
    return {(item.get("group"), item["metric"]) for item in selection["time_series"]}


def selected_time_series(cfg, group, metric):
    wanted = time_series_selection(cfg)
    return wanted is None or (group, metric) in wanted


def interaction_selection(cfg):
    """None = all legacy interaction metrics; empty tuple = explicit none.

    Each item is (kind, optional exact canonical interaction label, optional
    metric). Omitting one of the latter two selects all corresponding entries.
    """
    selection = cfg.get("statistics", {}).get("metric_selection")
    if selection is None or "interactions" not in selection:
        return None
    return tuple((item["kind"], item.get("interaction"), item.get("metric"))
                 for item in selection["interactions"])


def selected_interaction_metrics(selectors, kind, interaction):
    """Return the selected statistic names, preserving canonical metric order."""
    if selectors is None:
        return INTERACTION_METRICS
    matches = (metric for group, identity, metric in selectors
               if group == kind and (identity is None or identity == interaction))
    accepted = set()
    for metric in matches:
        accepted.update(INTERACTION_METRICS if metric is None else (metric,))
    return frozenset(accepted)


def simple_family_selection(cfg, family):
    """Return exact persisted identity tuples, or None for legacy all.

    A present empty list means disabled. The caller must preflight requested
    identities against saved inputs before publishing statistical outputs.
    """
    keys = SIMPLE_FAMILY_KEYS[family]
    selection = cfg.get("statistics", {}).get("metric_selection")
    if selection is None or family not in selection:
        return None
    return {tuple(item[key] for key in keys) for item in selection[family]}


def studio_metric_selection_defaults():
    """New Studio configs opt out of all inference; old YAML stays legacy-all.

    The Studio authoring layer should write this explicit mapping to YAML,
    rather than changing DEFAULTS and silently invalidating existing projects.
    """
    return {family: [] for family in (
        "time_series", "profiles", "interactions", "solvation_rdfs",
        "diffusion", "coupling", "lipid_contacts",
    )}


def metric_estimator_enabled(cfg, family, identity, estimator):
    """Resolve exact-identity estimator toggles without changing legacy mode.

    For listed metrics, explicit ci/hypothesis_tests settings override
    inference (default True). Broad interaction selectors are OR-combined.
    """
    if estimator not in ESTIMATOR_FLAGS[family]:
        raise ValueError(f"{family} does not support estimator {estimator}")
    selectors = cfg.get("statistics", {}).get("metric_selection")
    if selectors is None or family not in selectors:
        return True
    if family == "interactions":
        kind, interaction, metric = identity
        matches = [
            row for row in selectors[family]
            if row["kind"] == kind
            and (row.get("interaction") is None or row["interaction"] == interaction)
            and (row.get("metric") is None or row["metric"] == metric)
        ]
    else:
        keys = ("group", "metric") if family == "time_series" else SIMPLE_FAMILY_KEYS[family]
        matches = [
            row for row in selectors[family]
            if tuple(row.get(key) for key in keys) == tuple(identity)
        ]
    return any(row.get(estimator, row.get("inference", True)) for row in matches)

