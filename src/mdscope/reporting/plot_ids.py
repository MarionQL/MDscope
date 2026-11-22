"""Stable opt-in plot IDs for renderers that support individual selection.

No hidden broad matching: an include-list selects exact figure variants.
Unimplemented renderer families are rejected rather than silently ignored.
"""
from __future__ import annotations

from pathlib import Path


def safe_plot_stem(value) -> str:
    """Match the filename normalization used by production figure renderers."""
    return "".join(ch if ch.isalnum() or ch in "-_." else "_"
                   for ch in str(value)).strip("_")


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
SDC4_CONDITION_STEMS = frozenset({
    "lipid_region_moiety_contact_fraction",
    "lipid_region_mean_contact_fraction",
    "lipid_region_mean_local_enrichment",
    "lipid_region_mean_matched_local_contact_fraction",
    "lipid_region_mean_matched_local_enrichment",
    "lipid_density_maps", "lipid_enrichment_maps",
    "nxip_height_occlusion_joint", "vvhplvpl_proline_ramachandran",
})

STATIC_IDS = frozenset({
    "profiles:residue_response:heatmap",
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
        if (len(parts) == 3 and parts[0] == "sdc4" and parts[1]
                and parts[2] in SDC4_CONDITION_STEMS):
            continue
        if (len(parts) == 2 and parts[0] == "lipid_contacts" and
                parts[1].startswith(("lipid_contact_matrix_", "lipid_contact_difference_",
                                     "lipid_contact_class_matrix_", "lipid_contact_class_difference_"))
                and parts[1].endswith(("_pair_count_undirected",
                                       "_neighbors_per_molecule", "_enrichment"))):
            continue
        if (len(parts) == 3 and parts[0] in {"time_series", "profiles"}
                and parts[1] and parts[2] in DYNAMIC_VARIANTS[parts[0]]):
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


def describe_plot(identifier: str) -> dict:
    """Describe prerequisites for one supported plot without loading data."""
    validate_plot_ids([identifier])
    parts = identifier.split(":")
    family = parts[0]
    source = []
    destination = None
    unit_notes = "Metric-dependent units shown on output axes"
    if family == "time_series":
        stem, variant = parts[1:]
        source = ["{system}/{replicate}/summary/timeseries.parquet",
                  "{system}/{replicate}/summary/group_timeseries.parquet"]
        if variant in {"autocorrelation", "bootstrap_table", "comparison_table"}:
            source.append("project_summary/statistics/time_series/")
        destination = f"project_summary/figures/time_series/{stem}__{variant}.png"
    elif family == "profiles":
        if len(parts) == 4 and parts[-1] == "distance_heatmap":
            source = [f"{{system}}/{{replicate}}/profiles/residue_membrane_distance/"
                      f"{parts[1].removeprefix('residue_membrane_distance__')}/"
                      "distance_matrix_A.npy"]
            destination = f"project_summary/figures/profiles/{parts[1]}/{parts[2]}__distance_heatmap.png"
            unit_notes = "Distance in Å; time in ns"
        else:
            stem, variant = parts[1:]
            source = ["project_summary/statistics/profiles/replicate_profiles.parquet",
                      "project_summary/statistics/profiles/system_profiles.parquet"]
            destination = f"project_summary/figures/profiles/{stem}__{variant}.png"
    elif family == "interactions":
        kind = parts[1]
        if len(parts) == 3:
            source = [f"{{system}}/{{replicate}}/interactions/{kind}/"
                      f"{'occupancy' if parts[2] == 'occupancy_heatmap' else 'lifetimes'}.parquet"]
            destination = f"project_summary/figures/interactions/{kind}/{parts[2]}.png"
            unit_notes = "Occupancy fraction or observed duration in ns"
        else:
            source = ["project_summary/statistics/interactions/"]
            destination = f"project_summary/figures/interactions/{kind}/{parts[2]}__{parts[3]}.png"
    elif family == "solvation":
        stem = parts[1]
        source = {
            "shell_rdf": ["project_summary/statistics/solvation/system_rdf.parquet"],
            "water_and_ions": ["{system}/{replicate}/extensions/solvation/species_timeseries.parquet"],
            "contact_survival": ["{system}/{replicate}/extensions/solvation/survival.parquet"],
        }[stem]
        destination = f"project_summary/figures/solvation/{stem}.png"
        unit_notes = "RDF g(r), coordination number, or censored survival probability"
    elif family == "diffusion":
        source = ["project_summary/diffusion/replicate_coefficients.parquet"]
        destination = "project_summary/diffusion/diffusion_coefficients.png"
        unit_notes = "nm²/ns"
    elif family == "lipid_enrichment":
        source = ["project_summary/statistics/time_series/replicate_statistics.parquet"]
        destination = f"project_summary/figures/membrane/{parts[1]}.png"
        unit_notes = "Local/leaflet or area-corrected ratio (dimensionless)"
    elif family == "lipid_contacts":
        source = ["{system}/{replicate}/extensions/lipid_contacts/matrix.parquet",
                  "project_summary/statistics/lipid_contacts/system_statistics.parquet",
                  "project_summary/statistics/lipid_contacts/system_comparisons.parquet"]
        destination = f"project_summary/figures/lipid_contacts/{parts[1]}.png"
        unit_notes = ("Unique undirected molecule pairs/frame, neighbors per reference "
                      "molecule, or random-mixing enrichment; per leaflet, pointwise bootstrap")
    elif family == "shared_states":
        stem = parts[1]
        source = ["project_summary/shared_states/" + (
            "replicate_populations.parquet" if stem == "state_populations"
            else "sampling_diagnostics.parquet")]
        destination = f"project_summary/shared_states/{stem}.png"
    elif family == "sdc4":
        stem = parts[2]
        filename = (
            "lipid_region_specificity.parquet" if stem.startswith("lipid_region_") else
            "lipid_density_maps.npz" if stem.startswith("lipid_") else
            "timeseries.parquet"
        )
        source = [f"{parts[1]}/{{replicate}}/extensions/sdc4_mechanisms/{filename}"]
        destination = f"project_summary/figures/sdc4/{parts[1]}/{stem}.png"
    variant_labels = {
        "timeseries": "Value over time", "distribution": "Replica distributions",
        "autocorrelation": "Autocorrelation", "bootstrap_table": "Replicate and condition statistics",
        "comparison_table": "Between-condition comparisons", "profile": "Condition profile",
        "system_difference": "Between-condition profile difference", "distance_heatmap": "Residue distance heatmap",
        "occupancy_heatmap": "Interaction occupancy", "lifetime_distribution": "Contact duration distribution",
    }
    label = variant_labels.get(parts[-1], parts[-1].replace("_", " ").capitalize())
    if len(parts) >= 3:
        label += " · " + parts[1].replace("__group_", " / group ").replace("_", " ")
    return {
        "label": label,
        "id": identifier, "family": family,
        "renderer": {
            "time_series": "reporting.plotting.plot_time_series",
            "profiles": "reporting.plotting.plot_profiles",
            "interactions": "reporting.plotting.plot_interactions",
            "solvation": "reporting.solvation_figures.plot_solvation",
            "diffusion": "reporting.membrane_figures.plot_project_diffusion",
            "lipid_enrichment": "reporting.membrane_figures.plot_lipid_enrichment_matrix",
            "lipid_contacts": "reporting.membrane_figures.plot_lipid_contacts",
            "shared_states": "reporting.sdc4_figures.plot_sdc4_mechanisms",
            "sdc4": "reporting.sdc4_figures.plot_sdc4_mechanisms",
        }[family],
        "source_candidates": source,
        "destination": destination,
        "units_note": unit_notes,
        "status": "conditional_on_saved_data",
        "available_overrides": ["dpi"],
    }


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
        if name.startswith(("lipid_contact_matrix_", "lipid_contact_difference_",
                            "lipid_contact_class_matrix_", "lipid_contact_class_difference_")):
            return identifier
    if (len(parts) >= 5 and parts[-3] == "sdc4" and
            parts[-4] == "figures" and name in SDC4_CONDITION_STEMS):
        return f"sdc4:{parts[-2]}:{name}"
    for family in ("solvation", "diffusion", "lipid_enrichment"):
        identifier = f"{family}:{name}"
        if identifier in STATIC_IDS:
            return identifier
    return None
