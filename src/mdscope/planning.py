"""Read-only, conservative inventory of configured MDscope work.

The artifact inventory reuses registry.analysis_output_plan and
registry.required_artifacts.  Entries not registered there are explicitly
reported as *unmapped*, rather than guessed from file naming conventions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import analysis_runs, iter_replicates, plot_enabled, statistics_family_config
from .registry import analysis_output_plan, required_artifacts, artifact_schema
from .statistics_families import STATISTICS_FAMILIES


# Interpretation of coordinates, not a claim about PBC or bulk density.
# See docs/developer/scientific-contracts.md for the actual estimator contracts.
FITTED_ANALYSES = frozenset({"rmsd", "rmsf"})
SAVED_DATA_ANALYSES = frozenset({"shared_states", "coupling_screen"})
CONTROL_ANALYSES = frozenset({"alignment", "lifetimes"})
FIGURE_FAMILIES = (
    "time_series", "profiles", "interactions", "solvation", "sdc4",
    "shared_states", "diffusion", "lipid_enrichment", "lipid_contacts", "other",
)
STATISTICS_DIRS = {
    "time_series": "time_series",
    "profiles": "profiles",
    "interactions": "interactions",
    "solvation_rdfs": "solvation",
    "diffusion": "diffusion",
    "coupling": "coupling",
}


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
    "lipid_contacts": ("Per-molecule lipid neighborhood contacts", ("time_series",), ("lipid_contacts",)),
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
    "membrane_order": ("lipid-tail C-H sites; bond connectivity where available",),
    "sasa": ("target/occluder atoms and radii; optional FreeSASA backend",),
    "solvation": ("target chemical sites and water/ion identities",),
    "sdc4_mechanisms": ("SDC4 domain identities; topology bonds for hs_auto",),
    "hydration": ("water oxygens, target identities, topology bonds for glycans",),
    "residue_membrane_distance": ("query residues and membrane atom selection",),
    "bilayer_normal": ("lipid anchor selections and sufficient neighbors",),
    "lipid_enrichment": ("footprint and lipid anchor selections",),
    "hydrogen_bonds": ("donor/hydrogen/acceptor topology chemistry",),
    "lipid_hydrogen_bonds": ("lipid molecules and H-bond donor chemistry",),
}


# File names below come from the current statistics writers. They are possible
# products, not completeness requirements: some are only written with data.
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


def build_plan(cfg: dict, config_path: str | Path, *,
               validation_mode: str = "schema_only",
               topology_counts: dict | None = None) -> dict[str, Any]:
    """Describe enabled work using real registered primary/companion outputs.

    No disk access and no trajectory or bootstrap computations occur here.
    Topology counts, when supplied, must come from genuine topology preflight.
    """
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
        })

    project_statistics = []
    for family in STATISTICS_FAMILIES:
        project_statistics.append({
            "family": family,
            "status": ("enabled_conditional" if statistics_family_config(cfg, family) is not None
                       else "disabled"),
            "output_directory": f"project_summary/statistics/{STATISTICS_DIRS[family]}/",
            "files": [f"project_summary/statistics/{STATISTICS_DIRS[family]}/{file}"
                      for file in STATISTICS_CANDIDATES[family]],
            "files_status": "conditional_candidates_from_statistics_writers",
            "independent_condition_units": "simulation_replicates",
            "metric_selection": (
                cfg.get("statistics", {}).get("metric_selection", {}).get("time_series")
                if family == "time_series" and
                cfg.get("statistics", {}).get("metric_selection") is not None
                else None
            ),
            "metric_identity": (
                ["group", "metric"] if family == "time_series" else None
            ),
            "metric_identity_status": (
                "saved_data_dependent; no dynamic identities invented"
                if family == "time_series" else "not_yet_registered"
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
    for family in FIGURE_FAMILIES:
        # Solvation renderer has an additional explicit family opt-in.
        enabled = (plot_enabled(cfg, family) and
                   (family != "solvation" or
                    plots.get("families", {}).get("solvation", False)))
        figure_families.append({
            "family": family,
            "status": "enabled_conditional" if enabled else "disabled",
            "specific_figures": None,
            "specific_figures_status": "not_registered",
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

    return {
        "schema_version": 1,
        "validation": {
            "mode": validation_mode,
            "topology_checked": topology_counts is not None,
            "trajectory_frames_checked": False,
        },
        "project": cfg["project"],
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
        "topology_counts": topology_counts,
        "analyses": selected,
        "registered_replica_outputs": outputs,
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
