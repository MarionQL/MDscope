"""
Project orchestration: run replicates, consolidate outputs, then compute
statistics and figures.
"""

from __future__ import annotations

import glob
import json
import os
import tempfile
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
import pandas as pd
import yaml

from .. import __version__
from ..config import analysis_enabled, iter_replicates, load_config, project_root
from ..contracts.outputs import EXTENSION_ANALYSES, MEMBRANE_ANALYSES
from ..core.persistence import (
    atomic_json,
    atomic_parquet,
    ensure_dir,
    publish_directory,
)
from ..core.progress import stage_progress
from ..core.provenance import config_hash, environment_metadata
from ..errors import ConfigError, ReplicateIndexError
from ..interactions import postprocess_interactions
from ..mechanisms.residue_occlusion import (
    OCCLUSION_FILE,
    residue_environment_pending,
)
from ..observables.density_1d import missing_membrane_density_runs
from ..observables.hydration import temporal_hs_hydration_profiles
from ..observables.lipid_enrichment import (
    footprint_pending,
    summarize_lipid_composition,
)
from ..observables.membrane_order import membrane_order_timeseries
from ..observables.msd import (
    summarize_diffusion,
    summarize_project_diffusion,
    summarize_project_multi_origin,
)
from ..observables.multi_origin_pipeline import run_multi_origin_diffusion
from ..observables.structure_profiles import postprocess_structure_profiles
from ..registry import analysis_output_plan, fused_output, fused_output_plan
from ..statistics.dispatch import analyze_project
from ..statistics.shared_states import fit_shared_states, publish_shared_states
from ..systems.sdc4.mechanisms import complete_residue_occlusion
from ..systems.sdc4.mechanisms import enabled as sdc4_enabled
from .fused_all import run_fused_only
from .persisted_inputs import (
    check_artifact as _check_artifact,
)
from .persisted_inputs import (
    configuration_with_saved_sources,
    inspect_results_directory,
)
from .persisted_inputs import (
    replicate_state_prefix as _replicate_state_prefix,
)
from .persisted_inputs import (
    validate_analysis_output as _validate_analysis_output,
)
from .persisted_inputs import (
    validate_replica_products as _validate_replica_products,
)
from .preflight import validate_selections
from .replica_tasks import additional_time_series, specialized_tasks
from .streaming import (
    complete_membrane_density,
    run_streaming_pass,
    validate_stream_inventory,
)


class RunState:
    # Small completion markers let interrupted long trajectories resume
    # without repeating finished stages.
    def __init__(self, root: Path):
        self.root = ensure_dir(root / ".state")

    def marker(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def complete(self, key: str, required_artifacts: tuple[Path, ...] | None = None) -> bool:
        path = self.marker(key)
        if not path.exists():
            return False
        try:
            marker = json.loads(path.read_text())
            if not isinstance(marker, dict) or marker.get("status") != "complete":
                return False
            if required_artifacts is None:
                return True
            if "artifacts" not in marker:
                # Older markers can still establish completion if the actual
                # expected outputs are readable.
                for artifact in required_artifacts:
                    _check_artifact(artifact)
                return True
            recorded = marker["artifacts"]
            if not isinstance(recorded, dict):
                return False
            for artifact in required_artifacts:
                saved_size = recorded.get(str(artifact))
                if type(saved_size) is not int or _check_artifact(artifact) != saved_size:
                    return False
            return set(recorded) == {str(p) for p in required_artifacts}
        except Exception:
            # Corrupted markers, filesystem failures and invalid Parquet
            # metadata must invalidate reuse, not crash a restart.
            return False

    def mark(self, key: str, extra=None):
        atomic_json(
            self.marker(key),
            {**(extra or {}), "status": "complete"}
        )


def _create_once(path: Path, payload: str) -> bool:
    """Atomically publish a new file without ever replacing another task's file."""
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}-", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            return False
        return True
    finally:
        os.unlink(temporary)


def initialize(config_path: str | Path, cfg: dict):
    root = ensure_dir(project_root(cfg, config_path))
    resolved = root / "resolved_config.yaml"
    serialized = yaml.safe_dump(cfg, sort_keys=False)
    if not _create_once(resolved, serialized):
        try:
            existing = yaml.safe_load(resolved.read_text(encoding="utf-8"))
        except (OSError, yaml.YAMLError) as exc:
            raise RuntimeError(f"Cannot verify saved configuration {resolved}") from exc
        if config_hash(existing) != config_hash(cfg):
            raise ConfigError(
                f"Existing project {root} has config hash {config_hash(existing)}; "
                f"new run has {config_hash(cfg)}. Use a separate output directory."
            )
    metadata = environment_metadata()
    metadata.update({
        "config_path": str(Path(config_path).resolve()),
        "mdscope_version": __version__,
        "config_hash": config_hash(cfg),
    })
    _create_once(root / "run_metadata.json", json.dumps(metadata, indent=2) + "\n")
    return root, RunState(root)


def _load_parts(pattern: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(pattern)))
    return pd.concat(
        [pd.read_parquet(path) for path in files],
        ignore_index=True
    ) if files else pd.DataFrame()


def _merge_timeseries(tables: list[pd.DataFrame], group: bool = False) -> pd.DataFrame:
    # Merge once after each analysis has written its compact table; this
    # keeps analysis functions independent without reopening
    # trajectories.
    tables = [table for table in tables if table is not None and not table.empty]
    if not tables:
        return pd.DataFrame()
    keys = ["group", "time_ps"] if group else ["time_ps"]
    for source_index, table in enumerate(tables):
        if table.duplicated(keys).any():
            raise ValueError(
                f"Duplicate time-series keys in source {source_index}: {keys}"
            )
    out = tables[0].sort_values(keys)
    for table in tables[1:]:
        table = table.sort_values(keys)
        duplicate_metrics = [
            col for col in table.columns if col in out.columns and col not in keys
        ]
        if duplicate_metrics:
            raise ValueError(
                "Conflicting observable columns across time-series sources: "
                + ", ".join(sorted(duplicate_metrics))
            )
        if len(table.columns) > len(keys):
            out = out.merge(table, on=keys, how="outer")
    return out.sort_values(keys).reset_index(drop=True)


def summarize_replicate(repout: Path, cfg: dict):
    """
    Consolidate raw outputs; scientific statistics are handled uniformly
    later.
    """
    raw, summary = repout / "raw", ensure_dir(repout / "summary")
    base = _load_parts(raw / "timeseries-part-*.parquet")
    if not base.empty:
        if base.duplicated("time_ps").any():
            raise ValueError("Duplicate times in base replicate time series")
        base = base.sort_values("time_ps")
    group_base = _load_parts(raw / "group_timeseries-part-*.parquet")
    if not group_base.empty:
        if group_base.duplicated(["group", "time_ps"]).any():
            raise ValueError("Duplicate group/time keys in base replicate time series")
        group_base = group_base.sort_values(["group", "time_ps"])

    profile_ts, profile_group_ts = postprocess_structure_profiles(repout, cfg)
    time_sources = [base] + profile_ts + membrane_order_timeseries(repout)
    time_sources.extend(additional_time_series(repout,cfg))
    available_times = [
        table[
            "time_ps"
        ].to_numpy(
            float
        ) for table in time_sources if table is not None and not table.empty and "time_ps" in table
    ]
    if not available_times:
        for kind in ("hydrogen_bonds","lipid_hydrogen_bonds"):
            counts=repout/"interactions"/kind/"count_timeseries.parquet"
            if counts.exists():
                available_times.append(pd.read_parquet(counts)["time_ps"].to_numpy(float))
    all_times = (
        np.sort(
            np.unique(np.concatenate(available_times))
        ) if available_times else np.array([], float)
    )
    interaction_ts = postprocess_interactions(repout, cfg, all_times)

    whole = _merge_timeseries(time_sources + interaction_ts)
    grouped = _merge_timeseries(
        ([group_base] if not group_base.empty else []) + profile_group_ts,
        group=True
    )
    written = []
    if not whole.empty:
        path = summary / "timeseries.parquet"
        atomic_parquet(path, whole)
        written.append(path)
        if cfg.get("analyses",{}).get("membrane",{}).get("enabled",False) and cfg.get("analyses",{}).get("membrane",{}).get("msd",{}).get("enabled",False):
            summarize_diffusion(repout,cfg)
    if not grouped.empty:
        path = summary / "group_timeseries.parquet"
        atomic_parquet(path, grouped)
        written.append(path)
    return tuple(written)


def _output_pending(entry, cfg: dict) -> bool:
    """Missing product, or a lipid footprint whose committed reference.json differs."""
    if entry.name == "lipid_enrichment":
        return footprint_pending(entry.path, cfg["analyses"]["lipid_enrichment"])
    return not entry.path.is_file()


def _plan_fused_replay(cfg: dict, repout: Path):
    """Return configured output identities and the instances needing a frame pass.

    Called after existing-product preflight and the optional primary stream;
    validation failures are not treated as missing or recoverable outputs.
    """
    expected = fused_output_plan(cfg, repout)
    pending = tuple(entry for entry in expected if _output_pending(entry, cfg))
    return expected, pending


def run_replicate(rep, cfg: dict, root: Path, state: RunState):
    prefix = _replicate_state_prefix(rep)
    repout = ensure_dir(root / rep.system / rep.label)
    metadata = environment_metadata()
    metadata.update({
        "system": rep.system, "replicate": rep.label,
        "mdscope_version": __version__, "config_hash": config_hash(cfg),
    })
    _create_once(repout / "run_metadata.json", json.dumps(metadata, indent=2) + "\n")
    stage_timings = []
    # Resolve the registered inventory once for both preflight and the final
    # audit: no independent re-selection after the trajectory has run.
    registered = analysis_output_plan(cfg, repout)
    verified_artifacts = {}
    # Reject damaged existing products before a reader can overwrite them.
    # The new occlusion companion can be calculated independently if absent.
    for entry in registered:
        if entry.name == "lipid_enrichment" and not (entry.path.parent / "reference.json").is_file():
            continue  # An interrupted publication has no committed footprint.
        if entry.path.exists():
            _validate_analysis_output(entry, cfg, allow_missing_occlusion=True,
                                      verified=verified_artifacts)
    streaming_analyses = [
        "rmsd",
        "rmsf",
        "radius_of_gyration",
        "density_1d",
        "membrane",
        "salt_bridges",
        "ligand_contacts"
    ]
    streaming_required = any(analysis_enabled(
        cfg,
        name
    ) for name in streaming_analyses) or bool(cfg.get(
        "cache",
        {}
    ).get("aligned_coordinates", False))
    if streaming_required:
        key = prefix + "-stream"
        if state.complete(key):
            validate_stream_inventory(rep, repout, completed=True)
        else:
            t0 = perf_counter()
            state.mark(key, run_streaming_pass(rep, cfg, root))
            stage_timings.append({"stage": "streaming", "seconds": perf_counter() - t0})

    # One registry-derived decision serves observer-only runs and recovery.
    # Successful observers publish their markers before a later observer can
    # fail; incomplete siblings never receive completion markers.
    expected, pending = _plan_fused_replay(cfg, repout)
    pending_set = set(pending)
    pending_names = {entry.name for entry in pending}
    published = set()
    membrane_recovered = bool(pending_names.intersection(MEMBRANE_ANALYSES))
    extension_recovered = bool(pending_names.intersection((*EXTENSION_ANALYSES, "hydration")))
    distance_recovered = "residue_membrane_distance" in pending_names

    def mark_fused(entry):
        if not entry.path.is_file():
            raise RuntimeError(f"Missing fused {entry.name} output: {entry.path}")
        if entry in pending_set and entry not in published:
            # Validate the committed companions before marking this output.
            # SDC4 residue occlusion remains a separate completion stage.
            _validate_analysis_output(entry, cfg, allow_missing_occlusion=True,
                                      verified=verified_artifacts)
        key = prefix + "-" + entry.checkpoint_key
        if not state.complete(key):
            if entry.name in MEMBRANE_ANALYSES:
                source = ("streaming" if streaming_required and not membrane_recovered
                          else "fused-membrane-fallback")
            elif entry.name in (*EXTENSION_ANALYSES, "hydration"):
                source = ("streaming" if streaming_required and not extension_recovered
                          else "fused-extensions-fallback")
            elif entry.name == "residue_membrane_distance":
                source = ("streaming" if streaming_required and not distance_recovered
                          else "fused-distances-fallback")
            else:
                source = "streaming" if streaming_required else "fused-only"
            state.mark(key, {"fused_into": source})
        if entry in pending_set:
            published.add(entry)

    def mark_finalized_observer():
        # Observers finalize independently inside the shared source-frame pass.
        for entry in pending:
            if entry not in published and entry.path.is_file():
                mark_fused(entry)

    if pending:
        t0 = perf_counter()
        run_fused_only(
            rep, cfg, repout, only=pending,
            validate_summary_grid=streaming_required,
            on_finalized=mark_finalized_observer,
        )
        stage_timings.append({
            "stage": "fused-recovery" if streaming_required else "fused-only",
            "seconds": perf_counter() - t0,
        })

    for entry in expected:
        mark_fused(entry)

    # Specialized tasks follow the fused readiness rule: existing products
    # were validated before any reader ran, so only an absent registered
    # product schedules work. A stage marker alone neither proves nor
    # disproves completion.
    for task in specialized_tasks(cfg):
        if not any(_output_pending(entry, cfg) for entry in registered
                   if entry.name in task.produces):
            continue
        t0 = perf_counter()
        result = task.run(rep, cfg, repout)
        state.mark(prefix + "-" + task.key, result if isinstance(result, dict) else {})
        stage_timings.append({"stage": task.key, "seconds": perf_counter() - t0})
    msd_cfg = cfg.get("analyses", {}).get("membrane", {}).get("msd", {})
    if msd_cfg.get("enabled") and msd_cfg.get("multi_origin", False):
        targets = tuple(repout / "summary" / name for name in (
            "multi_origin_diffusion.parquet", "multi_origin_msd.parquet",
        ))
        if not state.complete(prefix + "-multi-origin-msd", required_artifacts=targets):
            t0 = perf_counter()
            result = run_multi_origin_diffusion(rep, cfg, repout)
            state.mark(prefix + "-multi-origin-msd", {
                **result, "artifacts": {str(p): _check_artifact(p) for p in targets},
            })
            stage_timings.append({"stage": "multi-origin-msd", "seconds": perf_counter() - t0})
    summary_inputs_changed = bool(stage_timings)
    if sdc4_enabled(cfg) and cfg["analyses"]["sdc4_mechanisms"].get("residue_occlusion"):
        entry = fused_output(repout, "sdc4_mechanisms")
        if residue_environment_pending(entry.path.parent / OCCLUSION_FILE, cfg["analyses"]["sdc4_mechanisms"]):
            t0 = perf_counter()
            complete_residue_occlusion(rep, cfg, repout)
            stage_timings.append({"stage": "residue-occlusion", "seconds": perf_counter() - t0})

    for entry in registered:
        _validate_analysis_output(entry, cfg, verified=verified_artifacts)
    whole = repout / "summary" / "timeseries.parquet"
    grouped = repout / "summary" / "group_timeseries.parquet"
    summary_paths = tuple(
        path for path, required in (
            (whole, streaming_required),
            (grouped, any((repout / "raw").glob("group_timeseries-part-*.parquet"))),
        )
        if required or path.is_file()
    )
    # A summary marker is reusable only when its recorded outputs still
    # exist and have valid Parquet footers. Do not infer completion from
    # the marker or from one plausible file alone.
    summary_reusable = state.complete(prefix + "-summary", required_artifacts=summary_paths)
    if summary_inputs_changed or not summary_reusable:
        t0 = perf_counter()
        summary_paths = summarize_replicate(repout, cfg)
        if streaming_required and repout / "summary" / "timeseries.parquet" not in summary_paths:
            raise RuntimeError(f"Streaming summary produced no time series: {repout}")
        state.mark(prefix + "-summary", {
            "artifacts": {str(p): _check_artifact(p) for p in summary_paths},
        })
        stage_timings.append({"stage": "summary", "seconds": perf_counter() - t0})
    density_runs = missing_membrane_density_runs(repout, cfg)
    if density_runs:
        t0 = perf_counter()
        complete_membrane_density(rep, cfg, repout, density_runs)
        stage_timings.append({"stage": "membrane-density", "seconds": perf_counter() - t0})
    if cfg["analyses"].get("sdc4_mechanisms", {}).get("hs_sugar_water_cutoff_A") is not None:
        times = pd.read_parquet(whole, columns=["time_ps"]).time_ps.to_numpy(float)
        temporal_hs_hydration_profiles(repout, cfg, times)
    if stage_timings:
        atomic_json(repout / "stage_timings.json", {
            "system": rep.system, "replicate": rep.label,
            "stages": stage_timings,
        })
        for item in stage_timings:
            print(f"TIMING {rep.system}/{rep.label} {item['stage']}: {item['seconds']:.3f} s", flush=True)


def finalize_project(config_path: str | Path, cfg=None, root=None, state=None, phase="statistics", *, results_dir=None):
    if phase not in {"statistics", "plotting"}:
        raise ValueError(f"Unknown finalization phase: {phase}")
    cfg = cfg or load_config(config_path)
    verbose = cfg.get("progress", {}).get("enabled", True)
    if results_dir is not None:
        if root is not None:
            raise ValueError("Specify either root or results_dir, not both")
        root = Path(results_dir).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"Existing results directory not found: {root}")
        # Never call initialize() here: it overwrites the saved resolved config
        # that is required to establish the provenance of existing products.
        inventory = inspect_results_directory(cfg, Path(config_path).resolve().parent, root)
        # The inventory includes saved-only analyses discovered from the
        # persisted configuration as well as currently scheduled analyses.
        # Require complete coverage before offline finalization can consume
        # either class of product.
        # Derive expected coverage from the current configuration, not the
        # inventory itself: a missing inventory row must never become an
        # implicitly successful preflight.
        from ..registry import saved_output_plan
        replicates_for_preflight = list(iter_replicates(cfg, config_path))
        expected = {
            (entry.name, entry.run_id or "default", rep.system, rep.label)
            for rep in replicates_for_preflight
            for entry in saved_output_plan(
                cfg, root / rep.system / rep.label)
        }
        by_key = {
            (row["analysis"], row["run_id"], row["condition"], row["replica"]): row
            for row in inventory
        }
        invalid = [
            (*key, by_key[key]["status"] if key in by_key else "missing")
            for key in sorted(expected)
            if key not in by_key or by_key[key]["status"] != "verified_existing"
        ]
        if invalid:
            raise RuntimeError(
                "Cannot finalize: required persisted replica data are missing or invalid. "
                "Existing results are not verified for every scheduled analysis/run/"
                "condition/replica; refusing downstream finalization: " +
                "; ".join("/".join(map(str, row)) +
                          (": " + by_key[row[:4]].get("reason", "") if row[:4] in by_key else "")
                          for row in invalid[:20]))
        if phase == "statistics":
            from .persisted_inputs import unavailable_saved_statistics
            if unavailable_saved_statistics(cfg, inventory, allow_scheduled=False):
                raise RuntimeError(
                    "Existing results lack the saved coordinate cache required for exact "
                    "temporal RMSF statistics. Existing outputs were preserved; restore "
                    "the saved cache or remove the unavailable exact request.")
        cfg = configuration_with_saved_sources(cfg, inventory, statistics=phase == "statistics")
    initialize_after_validation = root is None and phase == "statistics"
    if root is None:
        root = project_root(cfg, config_path)
    replicates = list(iter_replicates(cfg, config_path))
    if phase == "statistics":
        with stage_progress("Statistics replica validation", enabled=verbose):
            _validate_replica_products(cfg, root, replicates)
    else:
        if not cfg.get("plots", {}).get("enabled", False):
            raise ValueError("Set plots.enabled: true to generate diagnostic plots")
        from ..planning import resolve_figure_presets
        cfg, preset_omissions = resolve_figure_presets(cfg, root, config_path)
        selected_figures = cfg.get("plots", {}).get("include")
        if selected_figures is not None:
            # With an exact allowlist, raw-data plots can be regenerated before
            # statistics. Check their actual source contracts, not plot names.
            from ..config import figure_enabled
            from ..contracts.plot_validation import describe_plot
            from ..reporting.plot_ids import requires_finalization
            requires_statistics = any(
                describe_plot(identifier)["requires_statistics"]
                for identifier in selected_figures if figure_enabled(cfg, identifier)
            )
        else:
            # Legacy family-wide plotting can dispatch statistical figures.
            requires_statistics = True
        needs_finalization = selected_figures is None or any(
            requires_finalization(identifier)
            for identifier in selected_figures if figure_enabled(cfg, identifier))
        if needs_finalization:
            from ..statistics.manifest import read_statistics_manifest
            manifest = root / "project_summary/statistics/manifest.json"
            facts = read_statistics_manifest(manifest)
            if requires_statistics and facts.state == "missing":
                raise RuntimeError(
                    f"Plotting requires completed statistics for selected statistical figures: "
                    f"missing {manifest}. Run 'mdscope statistics' first; existing plots were preserved."
                )
            if facts.state == "invalid":
                raise RuntimeError(
                    f"Plotting requires completed statistics for selected statistical figures: "
                    f"invalid {manifest} ({facts.error}). Inspect the saved statistics manifest "
                    "before retrying; existing plots were preserved."
                )
            if (requires_statistics and not facts.complete) or facts.state == "incomplete":
                raise RuntimeError(
                    f"Plotting requires completed statistics for selected statistical figures: "
                    f"{manifest}: {facts.error or ('status ' + repr(facts.status))}. Finish statistics first; "
                    "existing plots were preserved."
                )
        if not requires_statistics:
            # Even descriptive figures must not read incomplete, stale or
            # corrupted per-replica products.
            _validate_replica_products(cfg, root, replicates)
    system_repouts = {
        system: [
            root / rep.system / rep.label
            for rep in replicates
            if rep.system == system
        ]
        for system in cfg["systems"]
    }
    with ExitStack() as preparation:
        if phase == "statistics":
            atomic_json(root / ".state/finalization.json", {
                "status": "running", "config_hash": config_hash(cfg),
            })
        shared_generation = None
        shared_preparation_seconds = 0.
        if phase == "statistics" and cfg["analyses"].get("shared_states", {}).get("enabled"):
            # Prepare the complete model before archiving the live summary.
            directory = preparation.enter_context(TemporaryDirectory(prefix=".shared-state-preflight-", dir=ensure_dir(root)))
            started = perf_counter()
            shared_generation = fit_shared_states(cfg, directory, system_repouts)
            shared_preparation_seconds = perf_counter() - started
        if initialize_after_validation:
            with stage_progress("Statistics project initialization", enabled=verbose):
                root, _ = initialize(config_path, cfg)
        generation_root = root
        if phase == "statistics":
            generation_root = Path(preparation.enter_context(
                TemporaryDirectory(prefix=".finalization-", dir=ensure_dir(root))))
        project_summary = generation_root / "project_summary"
        # Finalization always provides a project summary directory, even when
        # optional statistics and plots are disabled.
        ensure_dir(project_summary)
        # Finalization only consumes persisted tables/arrays: no
        # trajectories are reopened.
        timings_path = project_summary / "finalization_timings.json"
        if phase == "plotting" and timings_path.exists():
            try:
                finalization_timings = json.loads(timings_path.read_text()).get("stages", [])
            except (ValueError, OSError):
                finalization_timings = []
        else:
            finalization_timings = []
        def timed(stage, fn, *args, **kwargs):
            started = perf_counter()
            try:
                # Statistics reports its own active family and heartbeat.
                with stage_progress(f"Finalize {stage}", enabled=verbose and phase == "statistics"
                                    and stage != "statistics"):
                    return fn(*args, **kwargs)
            finally:
                seconds = perf_counter() - started
                if stage == "shared_states":
                    seconds += shared_preparation_seconds
                finalization_timings.append({"stage": stage, "seconds": seconds})
                atomic_json(project_summary / "finalization_timings.json", {
                    "stages": finalization_timings,
                    "total_seconds": sum(row["seconds"] for row in finalization_timings),
                })
                print(f"TIMING finalize {stage}: {seconds:.3f} s", flush=True)

        if phase != "plotting":
            if cfg.get("analyses", {}).get("lipid_enrichment", {}).get("enabled"):
                with stage_progress("Statistics lipid composition prepass", enabled=verbose):
                    for paths in system_repouts.values():
                        for repout in paths:
                            summarize_lipid_composition(repout)
            timed("statistics", analyze_project, cfg, root, system_repouts, output_root=generation_root)
            timed("shared_states", publish_shared_states, shared_generation, generation_root)
            timed("diffusion_summary", summarize_project_diffusion, generation_root, system_repouts)
            timed("multi_origin_summary", summarize_project_multi_origin, generation_root, system_repouts, cfg)
            if cfg.get("statistics", {}).get("enabled", True):
                from ..reporting.report import generate_report
                timed("findings", generate_report, project_summary / "statistics", project_root=root)
            with publish_directory(project_summary, root / "project_summary",
                                   root / ".archives/project_summary", "previous"):
                atomic_json(root / ".state/finalization.json", {
                    "status": "complete", "config_hash": config_hash(cfg),
                })
            if cfg.get("statistics", {}).get("enabled", True):
                print(f"FINDINGS: {root / 'project_summary/statistics/findings.md'}", flush=True)
            return root

        from ..reporting.figure_inventory import (
            figure_generation,
        )
        from ..reporting.membrane_figures import (
            plot_lipid_contacts,
            plot_lipid_enrichment_matrix,
            plot_project_diffusion,
        )
        from ..reporting.plotting import plot_all
        from ..reporting.sdc4_figures import plot_sdc4_mechanisms
        from ..reporting.solvation_figures import plot_solvation

        with figure_generation(cfg, root, system_repouts) as figure_root:
            if cfg.get('plots', {}).get('presets') is not None:
                atomic_json(figure_root / 'project_summary/figures/preset_selection.json', {
                    'presets': cfg['plots']['presets'], 'selected_ids': selected_figures,
                    'omitted': preset_omissions})
            timed("standard_plots", plot_all, cfg, figure_root, system_repouts)
            timed("lipid_enrichment_plots", plot_lipid_enrichment_matrix, cfg, figure_root)
            timed("sdc4_plots", plot_sdc4_mechanisms, cfg, figure_root, system_repouts)
            timed("lipid_contact_plots", plot_lipid_contacts, cfg, figure_root, system_repouts)
            timed("diffusion_plot", plot_project_diffusion, cfg, figure_root)
            timed("solvation_plots", plot_solvation, cfg, figure_root, system_repouts)
        return root


def run_one(config_path: str | Path, index: int):
    """
    Run one zero-based systemxreplicate entry; useful for external job
    arrays.
    """
    cfg = load_config(config_path)
    replicates = list(iter_replicates(cfg, config_path))
    if index < 0 or index >= len(replicates):
        raise ReplicateIndexError(
            f"replicate index {index} is out of range 0..{len(replicates) - 1}"
        )
    rep = replicates[index]
    validate_selections(cfg, Path(config_path).resolve().parent, systems={rep.system})
    root, state = initialize(config_path, cfg)
    run_replicate(rep, cfg, root, state)
    return root / rep.system / rep.label


def run_all(config_path: str | Path):
    cfg = load_config(config_path)
    validate_selections(cfg, Path(config_path).resolve().parent)
    root, state = initialize(config_path, cfg)
    for rep in iter_replicates(cfg, config_path):
        run_replicate(rep, cfg, root, state)
    finalize_project(config_path, cfg, root, state)
    return root
