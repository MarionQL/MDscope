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
from .statistics import STATISTICS_FAMILIES


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
        return "fitted"
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
            selected.append({
                "analysis": name,
                "run_id": run.run_id,
                "label": run.label,
                "coordinate_source": _source(name, cfg),
                "selection_inputs": _selection_inputs(run.spec),
                "artifact_catalog_status": status,
            })

    # Hydration can be implicitly scheduled by enabled SDC4 sugar-water metrics.
    if "hydration" in registered_by_name and not any(
            item["analysis"] == "hydration" for item in selected):
        selected.append({
            "analysis": "hydration", "run_id": "default",
            "label": "Hydration (required by SDC4)",
            "coordinate_source": "source",
            "selection_inputs": {},
            "artifact_catalog_status": "registered_implicit",
        })

    project_statistics = []
    for family in STATISTICS_FAMILIES:
        project_statistics.append({
            "family": family,
            "status": ("enabled_conditional" if statistics_family_config(cfg, family) is not None
                       else "disabled"),
            "output_directory": f"project_summary/statistics/{STATISTICS_DIRS[family]}/",
            "files": None,
            "files_status": "depends_on_saved_data",
            "independent_condition_units": "simulation_replicates",
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
            "Coordinate-frame labels describe input streams, not PBC or reference-model suitability.",
            "Missing/unavailable numerical values must not be silently treated as measured zeros.",
        ],
    }
