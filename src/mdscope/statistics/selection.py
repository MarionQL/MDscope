"""Explicit inference selection for persisted time-series (not raw analyses).

An absent selector preserves legacy family behavior. A present empty list
means no time-series inference. Identities are the existing saved
(group, metric) keys; no new metric naming convention is invented.
"""
from __future__ import annotations


def validate_metric_selection(stats):
    selection = stats.get("metric_selection")
    if selection is None:
        return
    if not isinstance(selection, dict) or set(selection) - {"time_series"}:
        raise ValueError("statistics.metric_selection currently supports time_series only")
    for family, entries in selection.items():
        if not isinstance(entries, list):
            raise ValueError(f"statistics.metric_selection.{family} must be a list")
        seen = set()
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
