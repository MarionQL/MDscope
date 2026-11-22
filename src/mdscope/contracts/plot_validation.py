"""Dependency-free plot ID and options contract."""
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


