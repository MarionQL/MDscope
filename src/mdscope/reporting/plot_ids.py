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
        "comparison_table", "bootstrap_table",
    }),
    "interactions": frozenset({
        "occupancy_heatmap", "lifetime_distribution",
        "bootstrap_table", "comparison_table",
    }),
}
SDC4_REPLICA_STEMS = frozenset({
    "lipid_region_moiety_contact_fraction",
    "lipid_region_mean_contact_fraction",
    "lipid_region_mean_local_enrichment",
    "lipid_region_mean_matched_local_contact_fraction",
    "lipid_region_mean_matched_local_enrichment",
    "lipid_density_maps", "lipid_enrichment_maps",
    "nxip_height_occlusion_joint", "vvhplvpl_proline_ramachandran",
})

STATIC_IDS = frozenset({
    "shared_states:state_populations",
    "shared_states:state_sampling_diagnostics",
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
        if (len(parts) == 4 and parts[0] == "sdc4" and parts[1] and parts[2]
                and parts[3] in SDC4_REPLICA_STEMS):
            continue
        if (len(parts) == 2 and parts[0] == "lipid_contacts" and
                parts[1].startswith(("lipid_contact_matrix_", "lipid_contact_difference_"))
                and parts[1].endswith(("_neighbors_per_molecule", "_enrichment"))):
            continue
        if len(parts) == 3 and parts[0] in {"time_series", "profiles"}:
            if parts[1] and parts[2] in DYNAMIC_VARIANTS[parts[0]]:
                continue
        if (len(parts) == 4 and parts[0] == "profiles" and parts[1]
                and parts[2] and parts[3] == "distance_heatmap"):
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


def validate_figure_options(options):
    """Only presentation resolution is configurable per ID for now."""
    if not isinstance(options, dict):
        raise ValueError("plots.options must map named plot IDs to option mappings")
    for identifier, spec in options.items():
        validate_plot_ids([identifier])
        if not isinstance(spec, dict) or set(spec) - {"dpi"}:
            raise ValueError(f"plots.options.{identifier} supports only dpi")
        if "dpi" in spec and (type(spec["dpi"]) is not int or
                              not 50 <= spec["dpi"] <= 600):
            raise ValueError(f"plots.options.{identifier}.dpi must be an integer from 50 to 600")


def figure_dpi(cfg, plot_id, default):
    if plot_id is None:
        return default
    return cfg.get("plots", {}).get("options", {}).get(plot_id, {}).get("dpi", default)


def id_from_figure_path(path: Path) -> str | None:
    """Recover stable IDs for the supported persisted output filenames."""
    path = Path(path)
    if path.suffix != ".png":
        return None
    parts = path.parts
    name = path.stem
    if (path.parent.name.startswith("residue_membrane_distance__")
            and name.endswith("__distance_heatmap")):
        system = name.removesuffix("__distance_heatmap")
        return f"profiles:{path.parent.name}:{system}:distance_heatmap"
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
    if "shared_states" in parts:
        identifier = f"shared_states:{name}"
        return identifier if identifier in STATIC_IDS else None
    if "lipid_contacts" in parts:
        identifier = f"lipid_contacts:{name}"
        if name.startswith(("lipid_contact_matrix_", "lipid_contact_difference_")):
            return identifier
    if (len(parts) >= 5 and parts[-2] == "sdc4_mechanisms"
            and parts[-3] == "extensions" and name in SDC4_REPLICA_STEMS):
        return f"sdc4:{parts[-5]}:{parts[-4]}:{name}"
    for family in ("solvation", "diffusion", "lipid_enrichment"):
        identifier = f"{family}:{name}"
        if identifier in STATIC_IDS:
            return identifier
    return None
