"""Dependency-free statistics metric selection and eligible profile contracts."""
import re

# Residue-resolved protein observables with comparable saved labels. Shared by
# YAML validation and the saved-data screen so eligibility has one definition.
RESIDUE_RESPONSE_PROFILE = re.compile(
    r"^(?:rmsf$|dssp_(?:helix|sheet|coil)_occupancy$|"
    r"sasa_mean(?:__|$)|sasa_buried_nm2__|"
    r"residue_occlusion__|residue_water_coordination__all$|"
    r"residue_membrane_distance__(?!.*__group_))"
)


# Every statistics family and its directory under project_summary/statistics/.
# The one family list for YAML validation, the dispatcher and read-only planning.
STATISTICS_FAMILIES = {
    "time_series": "time_series",
    "profiles": "profiles",
    "interactions": "interactions",
    "solvation_rdfs": "solvation",
    "diffusion": "diffusion",
    "coupling": "coupling",
    "lipid_contacts": "lipid_contacts",
}


SIMPLE_FAMILY_KEYS = {
    "profiles": ("profile",),
    "solvation_rdfs": ("solvent", "target_class", "site_class"),
    "diffusion": ("metric",),
    "coupling": ("mechanism",),
    "lipid_contacts": ("level", "leaflet", "species_a", "species_b", "metric"),
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
    "lipid_contacts": frozenset({"ci"}),
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




# File names below come from the current statistics writers. They are possible
# products, not completeness requirements: some are only written with data.
# Minimal writer-verified statistical identity and supported estimator catalog.
# Figure patterns are conditional, never guaranteed concrete figure IDs.
STATISTICS_METRIC_CONTRACTS = {
    "time_series": (["group", "metric"], "per-replica observable",
                    ["temporal_summary", "block_bootstrap", "hierarchical_ci",
                     "replicate_permutation", "bh_fdr"], ["time_series:{metric}:{variant}"]),
    "profiles": (["profile"], "per-replica position within profile",
                 ["position_summary", "temporal_block_bootstrap",
                  "hierarchical_ci", "replicate_permutation", "bh_fdr"],
                 ["profiles:{profile}:{variant}"]),
    "interactions": (["kind", "interaction", "metric"], "per-replica interaction metric",
                     ["occupancy_and_event_summary", "block_bootstrap",
                      "hierarchical_ci", "replicate_permutation", "bh_fdr"],
                     ["interactions:{kind}:{variant}",
                      "interactions:{kind}:{metric}:{variant}"]),
    "solvation_rdfs": (["solvent", "target_class", "site_class"],
                       "per-replica radial bin",
                       ["descriptive_mean_sd", "shell_replicate_contrast"],
                       ["solvation:shell_rdf"]),
    "diffusion": (["metric"], "per-replica shared-lag slope",
                  ["shared_lag_fit", "replicate_contrast"],
                  ["diffusion:diffusion_coefficients"]),
    "coupling": (["mechanism"], "per-replica mechanism pair",
                 ["block_correlation", "replicate_contrast"], []),
    "lipid_contacts": (["level", "leaflet", "species_a", "species_b", "metric"],
                       "one molecular pair class per leaflet and independent replicate",
                       ["temporal_block_bootstrap", "hierarchical_ci"],
                       ["lipid_contacts:lipid_contact_{matrix,difference}_*"]),
}

STATISTICS_CANDIDATES = {
    "time_series": ("replicate_statistics.parquet", "autocorrelation.parquet",
                    "system_statistics.parquet", "system_comparisons.parquet",
                    "hs_contact_episodes.parquet"),
    "profiles": ("replicate_profiles.parquet", "residue_regions.parquet",
                 "replicate_profile_statistics.parquet", "system_profiles.parquet",
                 "system_profile_statistics.parquet", "system_comparisons.parquet"),
    "interactions": ("replicate_interaction_metrics.parquet",
                     "system_interaction_metrics.parquet",
                     "system_interaction_comparisons.parquet"),
    "solvation_rdfs": ("replicate_rdf.parquet", "system_rdf.parquet",
                      "replicate_rdf_shells.parquet", "rdf_shell_comparisons.parquet"),
    "diffusion": ("replicate_common_lag_slopes.parquet",
                  "common_lag_comparisons.parquet"),
    "coupling": ("replicate_correlations.parquet", "system_comparisons.parquet"),
    "lipid_contacts": ("replicate_statistics.parquet", "system_statistics.parquet",
                       "system_comparisons.parquet"),
}


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
