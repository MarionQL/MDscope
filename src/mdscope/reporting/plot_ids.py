"""Stable opt-in plot IDs for renderers that support individual selection.

No hidden broad matching: an include-list selects exact figure variants.
Unimplemented renderer families are rejected rather than silently ignored.
"""
from __future__ import annotations

from pathlib import Path

DYNAMIC_VARIANTS = {
    "time_series": frozenset({
        "timeseries", "distribution", "autocorrelation",
        "bootstrap_table", "comparison_table",
    }),
    "profiles": frozenset({
        "profile", "distribution", "system_difference",
        "comparison_table",
    }),
    "interactions": frozenset({
        "occupancy_heatmap", "lifetime_distribution",
        "bootstrap_table", "comparison_table",
    }),
}
STATIC_IDS = frozenset({
    "solvation:shell_rdf",
    "solvation:water_and_ions",
    "solvation:contact_survival",
    "diffusion:diffusion_coefficients",
    "lipid_enrichment:lipid_enrichment_heatmap",
    "lipid_enrichment:lipid_area_density_heatmap",
})


def validate_plot_ids(include):
    if include is None:
        return
    if not isinstance(include, list) or any(not isinstance(x, str) for x in include):
        raise ValueError("plots.include must be a list of named plot ID strings")
    if len(set(include)) != len(include):
        raise ValueError("plots.include must not contain repeated plot IDs")
    for identifier in include:
        if identifier in STATIC_IDS:
            continue
        parts = identifier.split(":")
        if len(parts) == 3 and parts[0] in {"time_series", "profiles"}:
            if parts[1] and parts[2] in DYNAMIC_VARIANTS[parts[0]]:
                continue
        if len(parts) in {3, 4} and parts[0] == "interactions":
            if (parts[1] and (len(parts) == 3 and parts[2] in {
                    "occupancy_heatmap", "lifetime_distribution"
                } or len(parts) == 4 and parts[2] and parts[3] in {
                    "bootstrap_table", "comparison_table"
                })):
                continue
        raise ValueError(
            f"plots.include has unsupported plot ID {identifier!r}; "
            "use a supported family:metric:variant or static plot ID"
        )


def id_from_figure_path(path: Path) -> str | None:
    """Recover stable IDs for the supported persisted output filenames."""
    path = Path(path)
    if path.suffix != ".png":
        return None
    parts = path.parts
    name = path.stem
    if "time_series" in parts or "profiles" in parts:
        family = "time_series" if "time_series" in parts else "profiles"
        if "__" in name:
            stem, variant = name.rsplit("__", 1)
            if variant in DYNAMIC_VARIANTS[family]:
                return f"{family}:{stem}:{variant}"
    if "interactions" in parts and len(parts) >= 2:
        kind = path.parent.name
        if name in {"occupancy_heatmap", "lifetime_distribution"}:
            return f"interactions:{kind}:{name}"
        if "__" in name:
            metric, variant = name.rsplit("__", 1)
            if variant in {"bootstrap_table", "comparison_table"}:
                return f"interactions:{kind}:{metric}:{variant}"
    for family in ("solvation", "diffusion", "lipid_enrichment"):
        identifier = f"{family}:{name}"
        if identifier in STATIC_IDS:
            return identifier
    return None
