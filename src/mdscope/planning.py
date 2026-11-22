"""Read-only, conservative inventory of configured MDscope work.

The artifact inventory reuses registry.analysis_output_plan and
registry.required_artifacts.  Entries not registered there are explicitly
reported as *unmapped*, rather than guessed from file naming conventions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import (
    analysis_runs,
    iter_replicates,
    plot_enabled,
    project_root,
    statistics_family_config,
)
from .core.glycans import glycan_definition
from .registry import (
    analysis_output_plan,
    artifact_schema,
    configured_profile_variables,
    density_output_contracts,
    required_artifacts,
    scalar_contracts,
    state_feature_catalog,
    structural_profile_contracts,
)
from .reporting.plot_ids import (
    SDC4_CONDITION_STEMS,
    STATIC_IDS,
    describe_plot,
    safe_plot_stem,
)
from .contracts.metric_selection import STATISTICS_FAMILIES
from .contracts.plot_validation import FIGURE_FAMILIES
from .contracts.state_feature_spec import SDC4_FEATURES
from .statistics.selection import ESTIMATOR_FLAGS

# Interpretation of coordinates, not a claim about PBC or bulk density.
# See docs/developer/scientific-contracts.md for the actual estimator contracts.
FITTED_ANALYSES = frozenset({"rmsd", "rmsf"})
SAVED_DATA_ANALYSES = frozenset({"shared_states", "coupling_screen"})
CONTROL_ANALYSES = frozenset({"alignment", "lifetimes"})


# Per-analysis semantics describe *eligibility*, never promise supported CIs.
# Output file names and companions remain owned by registry.py.
ANALYSIS_METADATA = {
    "alignment": ("Coordinate alignment control", (), ()),
    "rmsd": ("Fitted-structure displacement", ("time_series",), ("time_series",)),
    "rmsf": ("Fitted-structure atom/residue fluctuation", ("profiles",), ("profiles",)),
    "radius_of_gyration": ("Structure radius of gyration", ("time_series",), ("time_series",)),
    "density_1d": ("Selected-atom density along one axis", ("profiles",), ("profiles",)),
    "membrane": ("Bilayer area, thickness and lipid MSD", ("time_series", "diffusion"), ("time_series", "diffusion")),
    "membrane_order": ("Lipid-tail C-H orientational order", ("profiles", "time_series"), ("profiles",)),
    "residue_membrane_distance": ("Residue-to-membrane minimum distance", ("profiles",), ("profiles",)),
    "dssp": ("Per-residue secondary structure", ("profiles",), ("profiles",)),
    "sasa": ("Surface accessibility with explicit occluders", ("profiles",), ("profiles",)),
    "salt_bridges": ("Protein salt-bridge contact events", ("interactions",), ("interactions",)),
    "hydrogen_bonds": ("Hydrogen-bond contact events", ("interactions",), ("interactions",)),
    "ligand_contacts": ("Ligand–protein contact events", ("interactions",), ("interactions",)),
    "lipid_hydrogen_bonds": ("Lipid-associated hydrogen-bond events", ("interactions",), ("interactions",)),
    "lifetimes": ("Contact persistence postprocessing", ("interactions",), ("interactions",)),
    "bilayer_normal": ("Local lipid normal orientations", ("time_series",), ("time_series",)),
    "lipid_enrichment": ("Protein footprint / leaflet lipid enrichment", ("time_series",), ("lipid_enrichment",)),
    "lipid_contacts": ("Per-molecule lipid neighborhood contacts", ("time_series", "lipid_contacts"), ("lipid_contacts",)),
    "local_environment": ("Protein–lipid/glycan neighborhood observables", ("time_series",), ("time_series",)),
    "coupling_screen": ("Saved-feature coupling diagnostics", ("coupling",), ()),
    "sdc4_mechanisms": ("SDC4-specific motif, domain and membrane mechanisms", ("time_series", "profiles"), ("sdc4",)),
    "potassium_bridge_audit": ("SDC4 acidic-patch–potassium–HS audit", (), ()),
    "solvation": ("Site-resolved water and ion coordination", ("solvation_rdfs", "time_series"), ("solvation",)),
    "hydration": ("Unique-water hydration per sugar/residue/group", ("profiles", "time_series"), ("profiles",)),
    "shared_states": ("Replicate-balanced shared-state inference", (), ("shared_states",)),
}
REQUIRED_TOPOLOGY = {
    "rmsd": ("alignment atom set",),
    "rmsf": ("alignment atom set",),
    "membrane": ("lipid headgroup selection and valid box",),
    "membrane_order": ("lipid-tail carbon chains and C-H sites; bonds verify tail naming where available",),
    "sasa": ("target/occluder atoms and radii; optional FreeSASA backend",),
    "solvation": ("target chemical sites and water/ion identities",),
    "sdc4_mechanisms": ("SDC4 domain identities; topology bonds for hs_auto",),
    "hydration": ("water oxygens, target identities; topology bonds, glycan residues and "
                  "attachment atoms for linear glycan chains",),
    "residue_membrane_distance": ("query residues and membrane atom selection",),
    "bilayer_normal": ("lipid anchor selections and sufficient neighbors",),
    "lipid_enrichment": ("footprint and lipid anchor selections",),
    "hydrogen_bonds": ("donor/hydrogen/acceptor topology chemistry",),
    "lipid_hydrogen_bonds": ("lipid molecules and H-bond donor chemistry",),
}


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


def _format(path: str) -> str:
    """Identify the *storage container*, not an invented column schema."""
    return {
        ".parquet": "Parquet",
        ".json": "JSON",
        ".npz": "NumPy NPZ",
        ".npy": "NumPy NPY",
    }.get(Path(path).suffix.lower(), "other")


def _relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _selection_inputs(spec: dict) -> dict[str, Any]:
    """Expose user selections without reinterpreting MDAnalysis expressions."""
    fields = (
        "selection", "target_selection", "occluder_selection",
        "membrane_selection", "membrane_center_selection",
        "headgroup_selection", "lipid_selection", "protein_selection",
        "water_selection", "water_oxygen_selection", "ligand_selection",
        "alignment_selection", "features",
    )
    result = {name: spec[name] for name in fields if spec.get(name) is not None}
    for name in (
        "selections", "species_selections", "targets", "tm_chains",
        "ectodomains", "hs_chains", "lipid_anchors",
        "ion_selections", "motifs",
    ):
        if spec.get(name):
            result[name] = spec[name]
    return result


def _source(name: str, cfg: dict) -> str:
    if name in FITTED_ANALYSES:
        # RMSD/RMSF use the source-frame features when structural fitting is
        # disabled. Calling them "fitted" in the plan would misstate provenance.
        return ("fitted" if cfg["analyses"]["alignment"].get("enabled", True)
                else "source")
    if name in SAVED_DATA_ANALYSES:
        return "saved_replica_products"
    if name in CONTROL_ANALYSES:
        return "control"
    if name == "membrane" and cfg["analyses"]["membrane"].get("msd", {}).get("multi_origin"):
        return "source_or_separately_synchronized_msd_source"
    return "source"


def _inference_design(cfg, reps):
    """Schema-only permutation resolution based on configured replicate labels."""
    from collections import Counter

    from .statistics.permutation_resolution import permutation_floor
    counts = Counter(rep.system for rep in reps)
    labels = {}
    for rep in reps:
        labels.setdefault(rep.system, set()).add(rep.label)
    tests = cfg.get("statistics", {}).get("hypothesis_tests", {}) or {}
    design = tests.get("design", "independent")
    pairs = []
    for index, left in enumerate(counts):
        for right in list(counts)[index + 1:]:
            n = len(labels[left] & labels[right]) if design == "paired" else None
            value, kind = permutation_floor(counts[left], counts[right], n_matched=n, tests=tests)
            pairs.append({
                "system_1": left, "system_2": right,
                "n_replicates_system_1": counts[left],
                "n_replicates_system_2": counts[right],
                "n_matched_pairs": n,
                "resolution_kind": "monte_carlo" if kind == "monte_carlo_resolution" else kind,
                "min_attainable_p": value if kind == "exact" else None,
                "monte_carlo_resolution": value if kind == "monte_carlo_resolution" else None,
            })
    return {"design": design, "comparisons": pairs}


def build_plan(cfg: dict, config_path: str | Path, *,
               validation_mode: str = "schema_only",
               topology_counts: dict | None = None) -> dict[str, Any]:
    """Describe enabled work using real registered primary/companion outputs.

    No disk access and no trajectory or bootstrap computations occur here.
    Topology counts, when supplied, must come from genuine topology preflight.
    """
    from .config.runs import structural_selection_requirements
    config_path = Path(config_path).resolve()
    reps = list(iter_replicates(cfg, config_path))
    template_root = Path("__replica__")
    entries = analysis_output_plan(cfg, template_root)

    outputs = []
    registered_by_name: dict[str, list[str]] = {}
    for entry in entries:
        required = required_artifacts(entry, cfg)
        relative_paths = [_relative(Path(p), template_root) for p in required]
        registered_by_name.setdefault(entry.name, []).extend(relative_paths)
        for index, path in enumerate(relative_paths):
            outputs.append({
                "analysis": entry.name,
                "run_id": entry.run_id or "default",
                "path_within_replica": path,
                "format": _format(path),
                "role": "primary" if index == 0 else "required_companion",
                "status": "registered_required",
                **artifact_schema(entry, Path(path)),
                "rows": None,  # No trajectory sampling has been inspected.
            })

    selected = []
    state_sources = state_feature_catalog(cfg)["sources"]
    order_spec = cfg["analyses"]["membrane_order"]
    for name in cfg["analyses"]:
        for run in analysis_runs(cfg, name):
            status = ("registered" if name in registered_by_name
                      else "saved_data_stage" if name in SAVED_DATA_ANALYSES
                      else "control" if name in CONTROL_ANALYSES
                      else "primary_artifact_not_registered")
            meta = ANALYSIS_METADATA.get(name)
            candidate_statistics = meta[1] if meta else ()
            candidate_figures = meta[2] if meta else ()
            selected.append({
                "id": f"{name}:{run.run_id}",
                "analysis": name,
                "run_id": run.run_id,
                "label": run.label,
                "purpose": meta[0] if meta else "Not cataloged",
                "coordinate_source": _source(name, cfg),
                "selection_inputs": _selection_inputs(run.spec),
                **({"state_model_inputs": {
                    "preset": run.spec.get("preset", "sdc4"),
                    "features": (run.spec.get("features") if run.spec.get("features") is not None else SDC4_FEATURES),
                    "feature_sources": run.spec.get("feature_sources", {}),
                    "source_prerequisites": {
                        alias: state_sources[source["source"]]
                        for alias, source in run.spec.get("feature_sources", {}).items()},
                    "alignment": run.spec.get("alignment", "exact"),
                    "scaling": run.spec.get("scaling", "standardize"),
                }} if name == "shared_states" else {}),
                **({"glycan_definitions": {
                    glycan: glycan_definition(cfg["selections"]["glycans"], glycan)
                    for glycan in run.spec.get("glycans") or []}}
                   if name == "hydration" else {}),
                # Effective naming shared by membrane_order and SDC4 local order.
                **({"tail_naming": {
                    "tail_preset": order_spec.get("tail_preset"), "tails": order_spec["tails"],
                    "lipid_overrides": order_spec.get("lipid_overrides") or {}}}
                   if name == "membrane_order" or (
                       name == "sdc4_mechanisms" and (run.spec.get("local_order") or {}).get("enabled"))
                   else {}),
                "required_topology": list(REQUIRED_TOPOLOGY.get(name, ())),
                "topology_requirements_status": ("partial_contract" if name in REQUIRED_TOPOLOGY
                                                 else "not_registered"),
                "statistics_eligibility": [
                    {"family": family,
                     "status": ("enabled_conditional" if statistics_family_config(cfg, family)
                                is not None else "disabled")}
                    for family in candidate_statistics
                ],
                "plot_family_eligibility": [
                    {"family": family,
                     "status": "enabled_conditional" if plot_enabled(cfg, family) else "disabled"}
                    for family in candidate_figures
                ],
                "artifact_catalog_status": status,
                "descriptive_figures": analysis_figure_candidates(cfg, name),
            })

    # Hydration can be implicitly scheduled by enabled SDC4 sugar-water metrics.
    if "hydration" in registered_by_name and not any(
            item["analysis"] == "hydration" for item in selected):
        meta = ANALYSIS_METADATA["hydration"]
        selected.append({
            "id": "hydration:default",
            "analysis": "hydration", "run_id": "default",
            "label": "Hydration (required by SDC4)",
            "purpose": meta[0],
            "coordinate_source": "source",
            "selection_inputs": {},
            "required_topology": list(REQUIRED_TOPOLOGY["hydration"]),
            "topology_requirements_status": "partial_contract",
            "statistics_eligibility": [
                {"family": family,
                 "status": ("enabled_conditional" if statistics_family_config(cfg, family)
                            is not None else "disabled")}
                for family in meta[1]
            ],
            "plot_family_eligibility": [
                {"family": family,
                 "status": "enabled_conditional" if plot_enabled(cfg, family) else "disabled"}
                for family in meta[2]
            ],
            "artifact_catalog_status": "registered_implicit",
            "descriptive_figures": analysis_figure_candidates(cfg, "hydration"),
        })

    from .statistics.contracts import statistics_output_contract
    project_statistics = []
    for family in STATISTICS_FAMILIES:
        project_statistics.append({
            "family": family,
            "status": ("enabled_conditional" if statistics_family_config(cfg, family) is not None
                       else "disabled"),
            "output_directory": f"project_summary/statistics/{STATISTICS_FAMILIES[family]}/",
            "files": [f"project_summary/statistics/{STATISTICS_FAMILIES[family]}/{file}"
                      for file in STATISTICS_CANDIDATES[family]],
            "files_status": "conditional_candidates_from_statistics_writers",
            "independent_condition_units": "simulation_replicates",
            "candidate_source_runs": [
                {
                    "analysis": item["analysis"],
                    "run_id": item["run_id"],
                    "selection_inputs": item["selection_inputs"],
                    "status": "candidate; exact metric identity requires saved data",
                }
                for item in selected
                if family in (ANALYSIS_METADATA.get(item["analysis"], ("", (), ()))[1])
            ],
            "output_schema": statistics_output_contract(family),
            "metric_contract": {
                "identity_fields": STATISTICS_METRIC_CONTRACTS[family][0],
                "row_grain": STATISTICS_METRIC_CONTRACTS[family][1],
                "supported_estimators": STATISTICS_METRIC_CONTRACTS[family][2],
                "figure_id_patterns": STATISTICS_METRIC_CONTRACTS[family][3],
                "figure_ids_status": "conditional; exact dynamic IDs require saved data",
                "output_columns_status": "identity fields only; complete writer schema not catalogued",
                "raw_analysis_independent": True,
                "figure_selection_independent": True,
                "descriptive_selection": "included_identity",
                "inferential_selection": "per_entry_inference_boolean_default_true",
                "independent_estimator_toggles": sorted(ESTIMATOR_FLAGS[family]),
                "estimator_flag_default": "inference:true (legacy); specific flags override",
                "unsupported_estimator_policy": "reject_before_execution",
                "unavailable_estimate_policy": "nullable_ci_p_and_q_not_fabricated",
                "eligibility_requirements": {
                    "replicate_independence": "simulation_replicates",
                    "temporal": ("correlated_frames; full_grid_and_valid_blocks_required"
                                 if family in {"time_series", "profiles", "interactions"}
                                 else "not_a_framewise_bootstrap"),
                    "between_condition_ci": (
                        "complete_eligible_replicates; paired/unpaired design honored"
                    ),
                    "dynamic_ids": "requires_saved_data; schema_only_candidates",
                },
                "configured_statistical_figure_ids": [
                    identifier for identifier in (cfg.get("plots", {}).get("include") or [])
                    if identifier.split(":", 1)[0] == (
                        "solvation" if family == "solvation_rdfs" else family
                    )
                ],
            },
            "metric_selection": (
                cfg.get("statistics", {}).get("metric_selection", {}).get(family)
                if family in {"time_series", "interactions", "profiles",
                               "solvation_rdfs", "diffusion", "coupling"} and
                cfg.get("statistics", {}).get("metric_selection") is not None
                else None
            ),
            "metric_identity": (
                ["group", "metric"] if family == "time_series"
                else ["kind", "interaction", "metric"] if family == "interactions"
                else ["profile"] if family == "profiles"
                else ["solvent", "target_class", "site_class"] if family == "solvation_rdfs"
                else ["metric"] if family == "diffusion"
                else ["mechanism"] if family == "coupling"
                else None
            ),
            "metric_identity_status": (
                "saved_data_dependent; no dynamic identities invented"
                if family in {"time_series", "interactions", "profiles",
                               "solvation_rdfs", "diffusion", "coupling"}
                else "not_yet_registered"
            ),
        })
    if (cfg.get("statistics", {}).get("matched_window_ps") is not None and
            statistics_family_config(cfg, "time_series") is not None):
        project_statistics.append({
            "family": "matched_time_series",
            "status": "enabled_conditional",
            "output_directory": "project_summary/statistics/matched_time_series/",
            "files": None,
            "files_status": "depends_on_common_matched_window",
            "independent_condition_units": "simulation_replicates",
        })

    plots = cfg.get("plots", {})
    figure_families = []
    include = plots.get("include")
    for family in FIGURE_FAMILIES:
        # Figure descriptions come from the production plot-ID registry.
        enabled = (plot_enabled(cfg, family) and
                   (family != "solvation" or
                    plots.get("families", {}).get("solvation", False)))
        if include is None:
            ids = sorted(id for id in STATIC_IDS if id.startswith(f"{family}:"))
            if family == "profiles":
                # The renderer creates one distance heatmap per configured run
                # and system. No topology/trajectory inspection is required.
                ids.extend(
                    f"profiles:residue_membrane_distance__{safe_plot_stem(run.run_id)}:"
                    f"{safe_plot_stem(system)}:distance_heatmap"
                    for run in analysis_runs(cfg, "residue_membrane_distance")
                    for system in cfg["systems"]
                )
            if (family == "sdc4" and
                    cfg["analyses"]["sdc4_mechanisms"].get("enabled", False)):
                ids.extend(
                    f"sdc4:{safe_plot_stem(system)}:{stem}"
                    for system in cfg["systems"]
                    for stem in sorted(SDC4_CONDITION_STEMS)
                )
        else:
            ids = [identifier for identifier in include
                   if identifier.startswith(f"{family}:")]
        figures = []
        for identifier in ids:
            descriptor = describe_plot(identifier)
            descriptor["requires_statistics"] = any(
                path.startswith("project_summary/statistics/")
                for path in descriptor["source_candidates"]
            )
            descriptor["availability"] = (
                "conditional_on_saved_data" if enabled else "disabled_by_configuration"
            )
            if identifier.startswith("sdc4:"):
                _, system, _ = identifier.split(":")
                configured = {
                    safe_plot_stem(name): name for name in cfg["systems"]
                }
                descriptor["analysis"] = "sdc4_mechanisms"
                descriptor["system"] = configured.get(system)
                if (descriptor["system"] is None or not
                        cfg["analyses"]["sdc4_mechanisms"].get("enabled", False)):
                    descriptor["availability"] = "unavailable_unconfigured_analysis_or_system"
            if identifier.startswith("profiles:") and identifier.endswith(":distance_heatmap"):
                _, stem, system, _ = identifier.split(":")
                configured_runs = {
                    f"residue_membrane_distance__{safe_plot_stem(run.run_id)}": run.run_id
                    for run in analysis_runs(cfg, "residue_membrane_distance")
                }
                matched_systems = {
                    safe_plot_stem(name): name for name in cfg["systems"]
                }
                descriptor["analysis"] = "residue_membrane_distance"
                descriptor["run_id"] = configured_runs.get(stem)
                descriptor["system"] = matched_systems.get(system)
                if descriptor["run_id"] is None or descriptor["system"] is None:
                    descriptor["availability"] = "unavailable_unconfigured_run_or_system"
            figures.append(descriptor)
        figure_families.append({
            "family": family,
            "status": "enabled_conditional" if enabled else "disabled",
            "specific_figures": figures,
            "specific_figures_status": (
                "static and configured distance-heatmap candidates; other dynamic IDs require saved data"
                if include is None else
                "configured_ids_only; saved metric identity not validated"
            ),
            "source": "saved_products_only",
        })

    analyses = cfg["analyses"]
    membrane_msd = analyses["membrane"].get("msd", {})
    additional_passes = []
    if analyses["membrane"].get("enabled") and membrane_msd.get("enabled") and membrane_msd.get("multi_origin"):
        additional_passes.append("multi_origin_msd")
    if analyses["potassium_bridge_audit"].get("enabled"):
        additional_passes.append("potassium_bridge_audit")
    if analyses["dssp"].get("enabled") and analyses["dssp"].get("backend") == "parallel":
        additional_passes.append("parallel_dssp")

    from .config import residue_response_settings
    response = residue_response_settings(cfg)
    response_status = ("disabled_statistics" if not cfg["statistics"].get("enabled", True)
                       or statistics_family_config(cfg, "profiles") is None
                       else "disabled" if not response["enabled"]
                       else "empty_selection" if response["profiles"] == [] or response["chains"] == []
                       else "conditional_saved_profile_contrasts")
    return {
        "residue_response": {
            "settings": response, "status": response_status,
            "exploratory": True, "trajectory_reads": False,
            "requires": ["replicate_profiles.parquet", "system_comparisons.parquet"],
            "manifest": "project_summary/statistics/profiles/residue_response_manifest.json",
        },
        "schema_version": 1,
        "inference_design": _inference_design(cfg, reps),
        "validation": {
            "mode": validation_mode,
            "topology_checked": topology_counts is not None,
            "trajectory_frames_checked": False,
        },
        "project": cfg["project"],
        "resolved_output_directory": str(project_root(cfg, config_path)),
        "systems": [
            {"id": system, "label": spec.get("label", system),
             "replicates": sum(rep.system == system for rep in reps)}
            for system, spec in cfg["systems"].items()
        ],
        "replicates": [
            {"system": rep.system, "id": rep.label, "index": rep.index,
             "topology": str(rep.topology), "trajectory": str(rep.trajectory),
             "msd_trajectory": str(rep.msd_trajectory) if rep.msd_trajectory else None}
            for rep in reps
        ],
        "selections": cfg["selections"],
        "structural_selection_requirements": structural_selection_requirements(cfg),
        "topology_counts": topology_counts,
        "analyses": selected,
        "registered_replica_outputs": outputs,
        "observable_contracts": scalar_contracts(cfg),
        "structural_profile_outputs": structural_profile_contracts(cfg) + density_output_contracts(cfg),
        "profile_variables": configured_profile_variables(cfg),
        "conditional_replica_outputs": [
            {"path_within_replica": "summary/timeseries.parquet",
             "format": "Parquet", "status": "written_if_nonempty",
             "columns": ["time_ps"], "columns_status": "partial"},
            {"path_within_replica": "summary/group_timeseries.parquet",
             "format": "Parquet", "status": "written_if_nonempty",
             "columns": ["group", "time_ps"], "columns_status": "partial"},
        ],
        "statistics_families": project_statistics,
        "figure_families": figure_families,
        "additional_trajectory_passes": additional_passes,
        "coverage": {
            "registered_artifact_contracts": True,
            "all_generated_files_cataloged": False,
            "parquet_and_npz_schemas_cataloged": False,
            "partial_verified_column_schemas": True,
            "individual_figure_ids_cataloged": False,
            "unregistered_analysis_names": sorted({
                item["analysis"] for item in selected
                if item["artifact_catalog_status"] == "primary_artifact_not_registered"
            }),
        },
        "notes": [
            "Registered artifact paths are per-replica relative paths, not completed files.",
            "Some optional raw/cache/profiles and specialized statistics products are not yet cataloged.",
            "Fixed-column schemas are described for verified producer paths only; dynamic columns remain unregistered.",
            "Enabled statistics/figures are conditional on supported saved inputs, not guaranteed results.",
            "Statistics file lists are writer-produced candidates; empty/unsupported families may write none.",
            "Coordinate-frame labels describe input streams, not PBC or reference-model suitability.",
            "Missing/unavailable numerical values must not be silently treated as measured zeros.",
        ],
    }


# Writer-owned saved identity columns. This mapping is for reading existing
# statistics, not an alternate identity or output schema.
SAVED_STATISTICS_IDENTITY = {
    "time_series": ("replicate_statistics.parquet", ("group", "metric")),
    "profiles": ("replicate_profile_statistics.parquet", ("profile",)),
    "interactions": (
        "replicate_interaction_metrics.parquet",
        ("interaction_type", "interaction", "metric"),
    ),
    "solvation_rdfs": (
        "replicate_rdf.parquet", ("solvent", "target_class", "site_class"),
    ),
    "diffusion": ("replicate_common_lag_slopes.parquet", ("metric",)),
    "coupling": ("replicate_correlations.parquet", ("mechanism",)),
}


def inspect_saved_statistics(root: str | Path, *, max_identities: int = 256) -> dict:
    """Inspect completed statistics read-only; bound high-cardinality discovery.

    Schemas come from actual Parquet metadata, not an independently maintained
    GUI column list. Only identity columns are scanned, stopping at the cap.
    """
    import json

    import pyarrow.parquet as pq

    if max_identities < 1:
        raise ValueError("max_identities must be positive")
    dest = Path(root) / "project_summary" / "statistics"
    from .statistics.manifest import read_statistics_manifest
    manifest_file = dest / "manifest.json"
    facts = read_statistics_manifest(manifest_file)
    if facts.state == "missing":
        raise FileNotFoundError(f"No saved statistics manifest: {manifest_file}")
    if facts.state == "invalid":
        raise ValueError(f"Invalid saved statistics manifest {manifest_file}: {facts.error}")
    if not facts.complete:
        raise ValueError("Statistics are incomplete; cannot claim an output catalog")
    result = {}
    for family, (identity_file, columns) in SAVED_STATISTICS_IDENTITY.items():
        folder = dest / STATISTICS_FAMILIES[family]
        if facts.family_disabled(family):
            result[family] = {'status':'disabled','files':{},'identities':[],
                              'identities_truncated':False}
            continue
        files = {}
        for path in sorted(folder.glob("*.parquet")):
            info = pq.ParquetFile(path)
            schema = info.schema_arrow
            files[path.name] = {
                "columns": schema.names,
                "column_types": {
                    field.name: str(field.type) for field in schema
                },
                "rows": info.metadata.num_rows,
                "path": str(path.relative_to(Path(root))),
                "row_grain": (
                    STATISTICS_METRIC_CONTRACTS[family][1]
                    if path.name == identity_file else None
                ),
                "row_grain_status": (
                    "family_identity_source"
                    if path.name == identity_file else "not_verified_per_file"
                ),
            }
        source = folder / identity_file
        if not source.is_file():
            result[family] = {"status": "not_materialized",
                              "files": files, "identities": [],
                              "identities_truncated": False}
            continue
        reader = pq.ParquetFile(source)
        available = set(reader.schema_arrow.names)
        if not set(columns) <= available:
            raise ValueError(f"Saved statistics {source} lacks identity columns")
        identities = []
        seen = set()
        truncated = False
        for batch in reader.iter_batches(batch_size=4096, columns=list(columns)):
            for row in zip(*(batch.column(i).to_pylist()
                             for i in range(len(columns)))):
                identity = tuple(None if item is None else str(item) for item in row)
                if identity in seen:
                    continue
                if len(identities) >= max_identities:
                    truncated = True
                    break
                seen.add(identity)
                identities.append({
                    **dict(zip(columns, identity)),
                    "id": (f"statistics:{family}:"
                           + json.dumps(identity, separators=(",", ":"), ensure_ascii=False)),
                })
            if truncated:
                break
        result[family] = {
            "status": "saved",
            "identity_columns": list(columns),
            "identity_id_format": "statistics:<family>:<JSON array of persisted identity values>",
            "source_analysis_run_status": "not inferred from statistics rows; see plan candidate_source_runs",
            "identities": identities,
            "identities_truncated": truncated,
            "files": files,
        }
    return {"status": "complete", "families": result,
            "identity_limit": max_identities}


def inspect_saved_plot_ids(cfg: dict, root: str | Path, config_path: str | Path,
                           *, max_identities: int = 256) -> dict:
    """Discover actual saved figure stems without trajectory reads or rendering.

    Column names are obtained from Parquet metadata; grouped series and
    profiles read only their identity column. Figure contracts remain owned
    by reporting.plot_ids, not by this inspector.
    """
    import pyarrow.parquet as pq

    from .reporting.plot_ids import DYNAMIC_VARIANTS

    if max_identities < 1:
        raise ValueError("max_identities must be positive")
    root = Path(root)
    found = set()
    truncated = False

    def add(family, stem, variants):
        nonlocal truncated
        for variant in sorted(variants):
            identifier = f"{family}:{stem}:{variant}"
            if identifier in found:
                continue
            if len(found) >= max_identities:
                truncated = True
                return
            found.add(identifier)

    replicas = list(iter_replicates(cfg, config_path))
    for rep in replicas:
        folder = root / rep.system / rep.label
        for filename, grouped in (("timeseries.parquet", False),
                                  ("group_timeseries.parquet", True)):
            path = folder / "summary" / filename
            if not path.is_file():
                continue
            source = pq.ParquetFile(path)
            names = source.schema_arrow.names
            metrics = [name for name in names if name not in {"time_ps", "group"}]
            if grouped:
                if "group" not in names:
                    raise ValueError(f"Saved group time series lacks group column: {path}")
                groups = set()
                for batch in source.iter_batches(columns=["group"]):
                    groups.update(str(name) for name in batch.column(0).to_pylist()
                                  if name is not None)
                for group in sorted(groups):
                    for metric in metrics:
                        add("time_series", safe_plot_stem(
                            f"{metric}__group_{group}"), DYNAMIC_VARIANTS["time_series"])
            else:
                for metric in metrics:
                    add("time_series", safe_plot_stem(metric),
                        DYNAMIC_VARIANTS["time_series"])
        interactions = folder / "interactions"
        if interactions.is_dir():
            for kind in sorted(interactions.iterdir()):
                if not kind.is_dir():
                    continue
                if (kind / "occupancy.parquet").is_file():
                    add("interactions", safe_plot_stem(kind.name), ("occupancy_heatmap",))
                if (kind / "lifetimes.parquet").is_file():
                    add("interactions", safe_plot_stem(kind.name), ("lifetime_distribution",))

    from .config import profile_dir
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        for system in cfg["systems"]:
            if any(
                all((directory / filename).is_file()
                    for filename in ("distance_matrix_A.npy", "time_ps.npy",
                                     "residues.parquet"))
                for rep in replicas if rep.system == system
                for directory in [profile_dir(
                    root / rep.system / rep.label,
                    "residue_membrane_distance", run.run_id)]
            ):
                identifier = (
                    f"profiles:residue_membrane_distance__{safe_plot_stem(run.run_id)}:"
                    f"{safe_plot_stem(system)}:distance_heatmap"
                )
                if len(found) < max_identities:
                    found.add(identifier)
                else:
                    truncated = True

    stats_root = root / "project_summary" / "statistics"
    manifest_file = stats_root / "manifest.json"
    from .statistics.manifest import read_statistics_manifest
    facts = read_statistics_manifest(manifest_file)
    if facts.state == "invalid":
        raise ValueError(f"Invalid saved statistics manifest {manifest_file}: {facts.error}")
    if facts.complete:
        profile_path = stats_root / "profiles" / "replicate_profiles.parquet"
        if not facts.family_disabled('profiles') and profile_path.is_file():
            source = pq.ParquetFile(profile_path)
            if "profile" not in source.schema_arrow.names:
                raise ValueError(f"Saved profile table lacks profile column: {profile_path}")
            profiles = set()
            for batch in source.iter_batches(columns=["profile"]):
                profiles.update(str(name) for name in batch.column(0).to_pylist()
                                if name is not None)
            for profile in sorted(profiles):
                add("profiles", safe_plot_stem(profile), DYNAMIC_VARIANTS["profiles"])
        interactions_dir = stats_root / "interactions"
        for filename, variant in (
            ("system_interaction_metrics.parquet", "bootstrap_table"),
            ("system_interaction_comparisons.parquet", "comparison_table"),
        ):
            path = interactions_dir / filename
            if facts.family_disabled('interactions') or not path.is_file():
                continue
            source = pq.ParquetFile(path)
            fields = ("interaction_type", "metric")
            if not set(fields).issubset(source.schema_arrow.names):
                raise ValueError(f"Saved interaction statistics lack identity columns: {path}")
            pairs = set()
            for batch in source.iter_batches(columns=list(fields)):
                pairs.update(
                    (str(kind), str(metric))
                    for kind, metric in zip(batch.column(0).to_pylist(),
                                            batch.column(1).to_pylist())
                    if kind is not None and metric is not None
                )
            for kind, metric in sorted(pairs):
                add("interactions", f"{safe_plot_stem(kind)}:{safe_plot_stem(metric)}",
                    (variant,))

    chosen = cfg.get("plots", {}).get("include")
    candidates = []
    for identifier in sorted(found):
        meta = describe_plot(identifier)
        meta["requires_statistics"] = any(
            source.startswith("project_summary/statistics/")
            for source in meta["source_candidates"]
        )
        meta["availability"] = (
            "conditional_on_saved_data" if plot_enabled(cfg, meta["family"])
            else "disabled_by_configuration"
        )
        if meta['requires_statistics']:
            family = {'solvation':'solvation_rdfs'}.get(meta['family'], meta['family'])
            if facts.family_disabled(family):
                meta['availability'] = 'unavailable_disabled_statistics_family'
            elif not facts.complete:
                meta['availability'] = 'unavailable_' + facts.state + '_statistics'
        candidates.append(meta)
    missing = [] if chosen is None else sorted(
        set(chosen) - found - set(STATIC_IDS)
    )
    return {
        "status": "saved_identity_discovery",
        "figures": candidates,
        "unconfirmed_requested_ids": missing,
        "unconfirmed_status": "not_observed_in_inspected_saved_products",
        "truncated": truncated,
        "max_identities": max_identities,
    }


def analysis_figure_candidates(cfg, name):
    """Writer-supported descriptive figure candidates, without saved-data reads."""
    identifiers = []
    if name in {"hydrogen_bonds", "lipid_hydrogen_bonds", "salt_bridges"}:
        identifiers.append(f"interactions:{name}:occupancy_heatmap")
        if cfg["analyses"].get("lifetimes", {}).get("enabled", True):
            identifiers.append(f"interactions:{name}:lifetime_distribution")
    if name == "residue_membrane_distance":
        identifiers.extend(f"profiles:residue_membrane_distance__{safe_plot_stem(run.run_id)}:{safe_plot_stem(system)}:distance_heatmap"
                           for run in analysis_runs(cfg, name) for system in cfg["systems"])
    if name == "solvation":
        identifiers.append("solvation:water_and_ions")
        from .registry import fused_output
        if any(path.name == "survival.parquet" for path in required_artifacts(fused_output(Path(), "solvation"), cfg)):
            identifiers.append("solvation:contact_survival")
    if name == "shared_states":
        identifiers += ["shared_states:state_populations", "shared_states:state_sampling_diagnostics"]
    if name == "sdc4_mechanisms":
        identifiers.extend(f"sdc4:{safe_plot_stem(system)}:{stem}" for system in cfg["systems"] for stem in sorted(SDC4_CONDITION_STEMS))
    if name == "lipid_contacts":
        identifiers.extend(f"lipid_contacts:lipid_contact_matrix_{system}_{leaflet}_{metric}"
                           for system in cfg["systems"] for leaflet in ("upper", "lower")
                           for metric in ("neighbors_per_molecule", "enrichment"))
    if name == "membrane" and cfg["analyses"]["membrane"].get("msd", {}).get("enabled"):
        identifiers.append("diffusion:diffusion_coefficients")
    return [dict(describe_plot(identifier), requires_statistics=False,
                 availability="conditional_on_saved_data") for identifier in identifiers]


def statistical_figure_candidates(cfg):
    """Exact configured identities only; unknown interaction metrics stay unresolved."""
    selection = cfg.get("statistics", {}).get("metric_selection") or {}
    identifiers = []
    for entry in selection.get("profiles", []) or []:
        stem = safe_plot_stem(entry["profile"])
        identifiers.extend(f"profiles:{stem}:{variant}" for variant in
                           ("profile", "distribution", "system_difference", "comparison_table", "bootstrap_table"))
    for entry in selection.get("interactions", []) or []:
        if entry.get("metric"):
            stem = f"{safe_plot_stem(entry['kind'])}:{safe_plot_stem(entry['metric'])}"
            identifiers.extend(f"interactions:{stem}:{variant}" for variant in
                               ("bootstrap_table", "comparison_table"))
    if selection.get("solvation_rdfs"):
        identifiers.append("solvation:shell_rdf")
    if any(item.get("metric", "").startswith("lipid_footprint_") for item in selection.get("time_series", []) or []):
        identifiers.extend(["lipid_enrichment:lipid_enrichment_heatmap", "lipid_enrichment:lipid_area_density_heatmap"])
    return [dict(describe_plot(identifier), requires_statistics=True,
                 availability="conditional_on_saved_data") for identifier in dict.fromkeys(identifiers)]


def authoring_catalog(cfg, *, include_inactive=False):
    """One catalog for authoring reconciliation and UI figure eligibility."""
    from copy import deepcopy

    from .registry import configured_profile_variables, scalar_contracts
    from .reporting.plot_ids import DYNAMIC_VARIANTS
    configured = deepcopy(cfg)
    from .config import normalize_analysis_booleans
    normalize_analysis_booleans(configured)
    if include_inactive:
        for spec in configured["analyses"].values():
            spec["enabled"] = True
            for key, run in (spec.get("runs") or {}).items():
                spec["runs"][key] = {**run, "enabled": True} if isinstance(run, dict) else {"enabled": True}
        # Ownership includes supported but currently disabled suboptions.
        membrane = configured["analyses"]["membrane"]
        for field in ("area_per_lipid", "thickness", "msd"):
            membrane.setdefault(field, {})["enabled"] = True
        for name in ("rmsd", "radius_of_gyration"):
            configured["analyses"][name].pop("groups", None)
        distance = configured["analyses"]["residue_membrane_distance"]
        for spec in [distance, *(distance.get("runs") or {}).values()]:
            if isinstance(spec, dict):
                spec["groups"] = True
        density = configured["analyses"]["density_1d"]
        for spec in [density, *(density.get("runs") or {}).values()]:
            if isinstance(spec, dict):
                for key, target in (spec.get("selections") or {}).items():
                    spec["selections"][key] = {**({"selection": target} if isinstance(target, str) else target), "groups": True}
    scalars = scalar_contracts(configured)
    profiles = configured_profile_variables(configured, include_optional=include_inactive)
    figures = {name: {} for name in configured["analyses"]}
    for variable in scalars:
        for figure in variable["figures"]:
            owner = "statistics" if figure["requires_statistics"] else variable["analysis"]
            figures.setdefault(owner, {})[figure["id"]] = figure
    for name, spec in configured["analyses"].items():
        if spec.get("enabled"):
            for figure in analysis_figure_candidates(configured, name):
                figures[name][figure["id"]] = figure
    for variable in profiles:
        for variant in sorted(DYNAMIC_VARIANTS["profiles"]):
            identifier = f"profiles:{safe_plot_stem(variable['profile'])}:{variant}"
            figures.setdefault("statistics", {})[identifier] = dict(describe_plot(identifier), requires_statistics=True)
    for figure in statistical_figure_candidates(configured):
        figures.setdefault("statistics", {})[figure["id"]] = figure
    if not include_inactive:
        selected = configured.get("statistics", {}).get("metric_selection")
        scalar_selection = None if selected is None else selected.get("time_series")
        statistics = {}
        for variable in scalars:
            if scalar_selection is None or any(item["metric"] == variable["metric"] and item.get("group") == variable["group"] for item in scalar_selection):
                for figure in variable["figures"]:
                    if figure["requires_statistics"]:
                        statistics[figure["id"]] = figure
        for figure in statistical_figure_candidates(configured):
            statistics[figure["id"]] = figure
        figures["statistics"] = statistics
    return {"scalars": scalars, "profiles": profiles,
            "motifs": configured["selections"].get("motifs", {}),
            "glycans": configured["selections"].get("glycans", {}),
            "state_features": state_feature_catalog(configured),
            "figures": {owner: list(items.values()) for owner, items in figures.items()}}
