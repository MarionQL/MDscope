"""Explicit inference selection for persisted series and interactions (not raw analyses).

An absent selector preserves legacy family behavior. A present empty list
means no time-series inference. Identities are the existing saved
(group, metric) keys; no new metric naming convention is invented.
"""
from __future__ import annotations

INTERACTION_METRICS = frozenset({
    "occupancy", "event_rate_per_ns", "mean_lifetime_ps",
    "median_lifetime_ps", "total_lifetime_ps",
})


def validate_metric_selection(stats):
    selection = stats.get("metric_selection")
    if selection is None:
        return
    if not isinstance(selection, dict) or set(selection) - {"time_series", "interactions"}:
        raise ValueError("statistics.metric_selection supports time_series and interactions only")
    for family, entries in selection.items():
        if not isinstance(entries, list):
            raise ValueError(f"statistics.metric_selection.{family} must be a list")
        seen = set()
        if family == "interactions":
            for entry in entries:
                if (not isinstance(entry, dict) or set(entry) - {"kind", "interaction", "metric"}
                        or not isinstance(entry.get("kind"), str) or not entry["kind"].strip()
                        or ("interaction" in entry and
                            (not isinstance(entry["interaction"], str) or not entry["interaction"].strip()))
                        or ("metric" in entry and
                            entry["metric"] not in INTERACTION_METRICS)):
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
            if (not isinstance(entry, dict) or set(entry) - {"metric", "group"}
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
