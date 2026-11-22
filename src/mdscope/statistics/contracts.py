"""Read-only writer-verified statistics columns; never inferred numerical values.

These are guaranteed base schemas for nonempty products. Estimator provenance,
coverage and optional chemistry/reference columns may extend them. Identities are
not independent samples: condition inference uses simulation replicas.
"""

import re

# Residue-resolved protein observables with comparable saved labels. Shared by
# YAML validation and the saved-data screen so eligibility has one definition.
RESIDUE_RESPONSE_PROFILE = re.compile(
    r"^(?:rmsf$|dssp_(?:helix|sheet|coil)_occupancy$|"
    r"sasa_mean(?:__|$)|sasa_buried_nm2__|"
    r"residue_occlusion__|residue_water_coordination__all$|"
    r"residue_membrane_distance__(?!.*__group_))"
)

BASE_COLUMNS = {
    "time_series": ("replicate_statistics.parquet", ["system", "replicate", "group", "metric", "mean", "median", "mean_ci_low", "mean_ci_high", "median_ci_low", "median_ci_high", "temporal_sampling_status", "bootstrap_draws", "ci_level"], "replicate / group / metric"),
    "profiles": ("system_profiles.parquet", ["system", "profile", "position", "position_label", "mean", "n_replicates", "mean_ci_low", "mean_ci_high", "bootstrap_draws", "ci_level"], "condition / profile / aligned position"),
    "interactions": ("replicate_interaction_metrics.parquet", ["system", "replicate", "interaction_type", "interaction", "metric", "value", "ci_low", "ci_high", "bootstrap_draws", "ci_level"], "replicate / interaction / metric"),
    "solvation_rdfs": ("replicate_rdf.parquet", ["system", "replicate", "solvent", "target_class", "site_class", "r_inner_A", "r_outer_A", "g_bulk"], "replicate / chemical identity / radial bin; descriptive g(r)"),
    "diffusion": ("replicate_common_lag_slopes.parquet", ["metric", "system", "replicate", "value"], "replicate / common-lag diffusion metric"),
    "coupling": ("replicate_correlations.parquet", ["system", "replicate", "mechanism", "x", "y", "pearson_block_r", "block_frames", "n_blocks", "interpretation"], "replicate / mechanism pair; exploratory association"),
    "lipid_contacts": ("replicate_statistics.parquet", ["system", "replicate", "level", "leaflet", "species_a", "species_b", "metric", "mean", "mean_ci_low", "mean_ci_high", "bootstrap_draws", "ci_level"], "replicate / lipid species or class / leaflet / pair / metric"),
}


def statistics_output_contract(family):
    """Return explicit base columns and limitations, without data or inference."""
    filename, columns, grain = BASE_COLUMNS[family]
    return {"filename": filename, "columns": list(columns),
            "columns_status": "verified_base_plus_conditional_provenance",
            "row_grain": grain, "independent_units": "simulation_replicates",
            "sample_rows": [{column: "illustrative placeholder" for column in columns} for _ in range(2)],
            "limitations": "Conditional file; unavailable CI/p/q remain nullable. No illustrative values are estimates. Full chemistry/provenance columns depend on configuration."}
