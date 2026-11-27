"""Validation rules for plots configuration."""
import math


def validate_plots(cfg):
    from ...contracts.plot_validation import FIGURE_FAMILIES
    allowed = set(FIGURE_FAMILIES)
    plots = cfg.get("plots", {})
    if not isinstance(plots.get("table_reports_pdf", False), bool):
        raise ValueError("plots.table_reports_pdf must be boolean")
    window_ns = plots.get("time_series_rolling_window_ns")
    if window_ns is not None and (
        type(window_ns) not in (int, float)
        or not math.isfinite(window_ns) or window_ns <= 0
    ):
        raise ValueError("plots.time_series_rolling_window_ns must be null or a positive finite number")
    families = plots.get("families", {})
    if (not isinstance(families, dict) or set(families) - allowed
            or any(not isinstance(value, bool) for value in families.values())):
        raise ValueError("plots.families must map known families to booleans")
    from ...contracts.plot_validation import validate_figure_options, validate_plot_ids
    presets = plots.get("presets")
    if presets is not None and (not isinstance(presets, list) or
            any(value not in {"primary", "secondary", "tertiary"} for value in presets) or
            len(set(presets)) != len(presets)):
        raise ValueError("plots.presets must be null or a unique list of primary, secondary, tertiary")
    validate_plot_ids(plots.get("include"))
    validate_figure_options(plots.get("options", {}))
    excluded = plots.get("interaction_exclude_kinds", [])
    if not isinstance(excluded, list) or any(not isinstance(value, str) for value in excluded):
        raise ValueError("plots.interaction_exclude_kinds must be a list")
