"""Select descriptive metrics and their eligible inferential estimators.

Missing selectors retain legacy family-level inference. Explicit empty lists
mean no statistics for that family; raw trajectory outputs are unaffected.
"""
from __future__ import annotations

# Canonical persisted identity fields for statistics-family selectors.
# Selectors identify whole scientific observables, never individual RDF bins,
# correlated frames, or implicitly inferred run IDs.
SIMPLE_FAMILY_KEYS = {
    "profiles": ("profile",),
    "solvation_rdfs": ("solvent", "target_class", "site_class"),
    "diffusion": ("metric",),
    "coupling": ("mechanism",),
}


# All families have descriptive point summaries. Whether to compute additional
# inferential products can be configured independently for each selected metric.
# The default remains True, so older selectors keep their original behavior.
# Capabilities reflect writer outputs, not GUI assumptions.
ESTIMATOR_FLAGS = {
    "time_series": frozenset({"ci", "hypothesis_tests"}),
    "profiles": frozenset({"ci", "hypothesis_tests"}),
    "interactions": frozenset({"ci", "hypothesis_tests"}),
    "solvation_rdfs": frozenset({"ci"}),
    "diffusion": frozenset({"ci"}),
    "coupling": frozenset({"ci"}),
}


def _valid_estimator_flags(family, entry):
    for key in ("inference", "ci", "hypothesis_tests"):
        if key in entry and not isinstance(entry[key], bool):
            return False
    if any(key in entry and key not in ESTIMATOR_FLAGS[family]
           for key in ("ci", "hypothesis_tests")):
        return False
    if entry.get("inference") is False and (
        entry.get("ci") is True or entry.get("hypothesis_tests") is True
    ):
        return False
    return True


INTERACTION_METRICS = frozenset({
    "occupancy", "event_rate_per_ns", "mean_lifetime_ps",
    "median_lifetime_ps", "total_lifetime_ps",
})


def validate_metric_selection(stats):
    selection = stats.get("metric_selection")
    if selection is None:
        return
    if not isinstance(selection, dict) or set(selection) - {"time_series", "interactions", *SIMPLE_FAMILY_KEYS}:
        raise ValueError("statistics.metric_selection contains unsupported family")
    for family, entries in selection.items():
        if not isinstance(entries, list):
            raise ValueError(f"statistics.metric_selection.{family} must be a list")
        seen = set()
        if family in SIMPLE_FAMILY_KEYS:
            keys = SIMPLE_FAMILY_KEYS[family]
            for entry in entries:
                if (not isinstance(entry, dict) or set(entry) - (set(keys) | {"inference", "ci", "hypothesis_tests"})
                        or any(key not in entry for key in keys)
                        or not _valid_estimator_flags(family, entry)
                        or any(not isinstance(entry.get(key), str) or not entry[key].strip()
                               for key in keys)):
                    raise ValueError(
                        f"statistics.metric_selection.{family} entries require "
                        f"nonempty {', '.join(keys)}"
                    )
                identity = tuple(entry[key] for key in keys)
                if identity in seen:
                    raise ValueError(f"Duplicate statistics.metric_selection.{family}: {identity!r}")
                seen.add(identity)
            continue
        if family == "interactions":
            for entry in entries:
                if (not isinstance(entry, dict) or set(entry) - {"kind", "interaction", "metric", "inference", "ci", "hypothesis_tests"}
                        or not _valid_estimator_flags(family, entry)
                        or not isinstance(entry.get("kind"), str) or not entry["kind"].strip()
                        or ("interaction" in entry and
                            (not isinstance(entry["interaction"], str) or not entry["interaction"].strip()))
                        or ("metric" in entry and
                            (not isinstance(entry["metric"], str)
                             or entry["metric"] not in INTERACTION_METRICS))):
                    raise ValueError(
                        "statistics.metric_selection.interactions entries require a "
                        "nonempty kind, optional exact interaction and a supported metric"
                    )
                key = (entry["kind"], entry.get("interaction"), entry.get("metric"))
                if key in seen:
                    raise ValueError(f"Duplicate statistics.metric_selection.interactions: {key!r}")
                seen.add(key)
            continue
        for entry in entries:
            if (not isinstance(entry, dict) or set(entry) - {"metric", "group", "inference", "ci", "hypothesis_tests"}
                    or not _valid_estimator_flags(family, entry)
                    or not isinstance(entry.get("metric"), str)
                    or not entry["metric"].strip()
                    or (entry.get("group") is not None
                        and (not isinstance(entry["group"], str) or not entry["group"].strip()))):
                raise ValueError(
                    "statistics.metric_selection.time_series entries require a nonempty "
                    "metric and optional nonempty group"
                )
            key = (entry.get("group"), entry["metric"])
            if key in seen:
                raise ValueError(f"Duplicate statistics.metric_selection.time_series: {key!r}")
            seen.add(key)


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
        "diffusion", "coupling",
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
