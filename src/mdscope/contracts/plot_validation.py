"""Dependency-free plot ID and options contract."""
# Figure families accepted by plots.families and catalogued by planning.
FIGURE_FAMILIES = (
    "time_series", "profiles", "interactions", "solvation", "sdc4",
    "shared_states", "diffusion", "lipid_enrichment", "lipid_contacts", "other",
)
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



def safe_plot_stem(value) -> str:
    """Match the filename normalization used by production figure renderers."""
    return "".join(ch if ch.isalnum() or ch in "-_." else "_"
                   for ch in str(value)).strip("_")


def describe_plot(identifier: str, *, table_reports_pdf=False) -> dict:
    """Describe prerequisites for one supported plot without loading data."""
    validate_plot_ids([identifier])
    parts = identifier.split(":")
    family = parts[0]
    source = []
    destination = None
    requires_statistics = False
    unit_notes = "Metric-dependent units shown on output axes"
    if family == "time_series":
        stem, variant = parts[1:]
        source = ["{system}/{replicate}/summary/timeseries.parquet",
                  "{system}/{replicate}/summary/group_timeseries.parquet"]
        requires_statistics = variant in {"autocorrelation", "bootstrap_table", "comparison_table"}
        if requires_statistics:
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
            requires_statistics = True
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
            requires_statistics = True
            source = ["project_summary/statistics/interactions/"]
            destination = f"project_summary/figures/interactions/{kind}/{parts[2]}__{parts[3]}.png"
    elif family == "solvation":
        stem = parts[1]
        requires_statistics = stem == "shell_rdf"
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
        requires_statistics = True
        source = ["project_summary/statistics/time_series/replicate_statistics.parquet"]
        destination = f"project_summary/figures/membrane/{parts[1]}.png"
        unit_notes = "Local/leaflet or area-corrected ratio (dimensionless)"
    elif family == "lipid_contacts":
        # Species neighbor/enrichment matrices have a descriptive saved-data
        # fallback. Class, pair-count and comparison figures require inference.
        requires_statistics = not (
            parts[1].startswith("lipid_contact_matrix_")
            and parts[1].endswith(("_neighbors_per_molecule", "_enrichment")))
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
    table_source = parts[-1] in {'bootstrap_table', 'comparison_table'}
    return {
        **({'consolidated_destination': 'project_summary/figures/statistical_reports/{registered_analysis}/statistics_tables.pdf',
            'csv_companion': 'statistics_tables.csv',
            'selection_semantics': 'Exact source rows within one analysis report; legacy PNG destination applies only with table_reports_pdf: false'} if table_source and table_reports_pdf else {}),
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
        "requires_statistics": requires_statistics,
        "destination": destination,
        "units_note": unit_notes,
        "status": "conditional_on_saved_data",
        "available_overrides": ["dpi"],
    }

