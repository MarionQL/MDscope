"""
Project orchestration: run replicates, consolidate outputs, then compute
statistics and figures.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import tempfile
from time import perf_counter
from contextlib import ExitStack
from tempfile import TemporaryDirectory
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from ..core.persistence import archive_output as _archive_previous_output
from ..core.validation import validate_selections

from .. import __version__
from ..errors import ConfigError, ReplicateIndexError
from ..config import (
    analysis_enabled,
    analysis_runs,
    atomic_json,
    atomic_parquet,
    atomic_text,
    config_hash,
    ensure_dir,
    environment_metadata,
    iter_replicates,
    load_config,
)
from ..interactions import run_hydrogen_bonds
from ..interactions import postprocess_interactions
from ..observables.membrane_order import membrane_order_timeseries, run_membrane_order
from ..observables.residue_membrane_distance import run_residue_membrane_distance
from ..observables.structure_profiles import postprocess_structure_profiles, run_structure_profiles
from ..statistics import analyze_project
from .persisted_inputs import (
    check_artifact as _check_artifact,
    validate_analysis_output as _validate_analysis_output,
    replicate_state_prefix as _replicate_state_prefix,
    validate_replica_products as _validate_replica_products,
    inspect_results_directory, configuration_with_saved_sources,
)
from .progress import stage_progress
from .streaming import run_streaming_pass, validate_stream_inventory, complete_membrane_density
from ..observables.density_1d import missing_membrane_density_runs
from .fused_membrane import enabled_names as fused_enabled_names, run_fused_membrane
from .fused_extensions import enabled_names as extension_enabled_names, run_fused_extensions
from .fused_distance import enabled_names as distance_enabled_names, run_fused_distances
from .fused_all import run_fused_only
from ..mechanisms.sdc4_mechanisms import enabled as sdc4_enabled, run_sdc4_mechanisms, complete_residue_occlusion
from ..mechanisms.residue_occlusion import (
    OCCLUSION_FILE, residue_environment_pending,
)
from ..statistics.shared_states import fit_shared_states, publish_shared_states
from ..mechanisms.hs_hydration import temporal_hs_hydration_profiles
from ..registry import fused_output, fused_output_plan, analysis_output_plan, required_artifacts
from .replica_tasks import specialized_tasks, additional_time_series
from ..observables.lipid_enrichment import footprint_pending, summarize_lipid_composition
from ..observables.msd import summarize_diffusion, summarize_project_diffusion, summarize_project_multi_origin
from ..observables.multi_origin_pipeline import run_multi_origin_diffusion


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


def project_root(cfg: dict, config_path: str | Path) -> Path:
    path = Path(cfg["project"]["output"])
    return (
        path
        if path.is_absolute()
        else Path(config_path).resolve().parent / path
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
    # Reject damaged existing products before a reader can overwrite them.
    # The new occlusion companion can be calculated independently if absent.
    for entry in analysis_output_plan(cfg, repout):
        if entry.name == "lipid_enrichment" and not (entry.path.parent / "reference.json").is_file():
            continue  # An interrupted publication has no committed footprint.
        if entry.path.exists():
            _validate_analysis_output(entry, cfg, allow_missing_occlusion=True)
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

    # Observer-only mode replays every pending instance in one frame pass.
    # Preserve run IDs: multiple SASA and distance runs have distinct outputs.
    # Lipid footprint needs its committed reference.json as well as a table.
    if not streaming_required:
        pending = tuple(
            entry for entry in fused_output_plan(cfg, repout)
            if (footprint_pending(entry.path, cfg["analyses"]["lipid_enrichment"])
                if entry.name == "lipid_enrichment" else not entry.path.is_file())
        )
        if pending:
            t0 = perf_counter()
            run_fused_only(rep, cfg, repout, only=pending)
            stage_timings.append({"stage": "fused-only", "seconds": perf_counter() - t0})

    # These analyses are observers of the existing streaming frames, not
    # separate trajectory readers. For geometry-only systems, one shared
    # fallback pass is sufficient. Keep legacy task markers for resume.
    fused_names = fused_enabled_names(cfg)
    if fused_names:
        expected = {name: fused_output(repout, name) for name in fused_names}
        # A completed streaming marker guarantees finalized observer outputs.
        # If outputs are missing (or there was no protein streaming pass),
        # rebuild them with exactly one independent shared pass.
        pending_names = tuple(
            name for name in fused_names
            if (footprint_pending(expected[name].path, cfg["analyses"][name])
                if name == "lipid_enrichment" else not expected[name].path.is_file())
        )
        need_fallback = bool(pending_names)
        if need_fallback:
            t0 = perf_counter()
            run_fused_membrane(rep, cfg, repout, names=pending_names)
            stage_timings.append({"stage": "fused-membrane-fallback", "seconds": perf_counter() - t0})
        for name in fused_names:
            if not expected[name].path.is_file():
                raise RuntimeError(f"Missing fused membrane output: {expected[name].path}")
            key = prefix + "-" + expected[name].checkpoint_key
            if not state.complete(key):
                state.mark(key, {"fused_into": "streaming" if streaming_required and not need_fallback
                                 else "fused-membrane-fallback"})

    extension_names = extension_enabled_names(cfg)
    if extension_names:
        extension_outputs = {name: fused_output(repout, name)
                             for name in extension_names if name != "sasa"}
        if "sasa" in extension_names:
            from ..config import analysis_runs
            for run in analysis_runs(cfg, "sasa"):
                extension_outputs["sasa" if run.run_id == "default" else f"sasa-{run.run_id}"] = (
                    fused_output(repout, "sasa") if run.run_id == "default"
                    else fused_output(repout, "sasa", run_id=run.run_id)
                )
        missing = [entry for entry in extension_outputs.values()
                   if not entry.path.is_file()]
        if missing:
            t0 = perf_counter()
            run_fused_extensions(
                rep, cfg, repout,
                names=tuple(dict.fromkeys(entry.name for entry in missing)),
                sasa_run_ids=tuple(
                    entry.run_id if entry.run_id is not None else "default"
                    for entry in missing if entry.name == "sasa"
                ),
            )
            stage_timings.append({"stage": "fused-extensions-fallback", "seconds": perf_counter() - t0})
        for name in extension_outputs:
            if not extension_outputs[name].path.is_file():
                raise RuntimeError(f"Missing fused output: {extension_outputs[name].path}")
            key = prefix + "-" + extension_outputs[name].checkpoint_key
            if not state.complete(key):
                state.mark(key, {"fused_into": "streaming" if streaming_required and not missing
                                 else "fused-extensions-fallback"})

    distance_names = distance_enabled_names(cfg)
    if distance_names:
        expected = {
            run_id: fused_output(repout, "residue_membrane_distance", run_id=run_id)
            for run_id in distance_names
        }
        missing = [name for name, entry in expected.items() if not entry.path.is_file()]
        if missing:
            t0 = perf_counter()
            run_fused_distances(rep, cfg, repout, run_ids=tuple(missing))
            stage_timings.append({"stage": "fused-distances-fallback", "seconds": perf_counter() - t0})
        for run_id, entry in expected.items():
            if not entry.path.is_file():
                raise RuntimeError(f"Missing fused residue-distance output: {entry.path}")
            key = prefix + "-" + entry.checkpoint_key
            if not state.complete(key):
                state.mark(key, {"fused_into": "streaming" if streaming_required and not missing
                                 else "fused-distances-fallback"})

    if sdc4_enabled(cfg):
        output = fused_output(repout, "sdc4_mechanisms").path
        if not output.is_file():
            t0 = perf_counter()
            run_sdc4_mechanisms(rep, cfg, repout)
            stage_timings.append({"stage": "sdc4-mechanisms-fallback",
                                  "seconds": perf_counter() - t0})
        if not output.is_file():
            raise RuntimeError(f"Missing SDC4 mechanism output: {output}")
        key = prefix + "-" + fused_output(repout, "sdc4_mechanisms").checkpoint_key
        if not state.complete(key):
            state.mark(key, {"fused_into": "streaming" if streaming_required
                             else "fused-only"})

    for task in specialized_tasks(cfg):
        key=prefix+"-"+task.key
        if task.key == "potassium_bridge_audit":
            # The complete audit directory is published atomically. Reuse it
            # even if the process stopped before recording the stage marker.
            if fused_output(repout, task.key).path.is_file():
                continue
        if task.key == "potassium_bridge_audit" or not state.complete(key):
            t0 = perf_counter()
            result=task.run(rep,cfg,repout)
            state.mark(key,result if isinstance(result,dict) else {})
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

    for entry in analysis_output_plan(cfg, repout):
        _validate_analysis_output(entry, cfg)
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
        selected_figures = cfg.get("plots", {}).get("include")
        if selected_figures is not None:
            # With an exact allowlist, raw-data plots can be regenerated before
            # statistics. Check their actual source contracts, not plot names.
            from ..reporting.plot_ids import describe_plot
            from ..config import figure_enabled
            requires_statistics = any(
                any(source.startswith("project_summary/statistics/")
                    for source in describe_plot(identifier)["source_candidates"])
                for identifier in selected_figures if figure_enabled(cfg, identifier)
            )
        else:
            # Legacy family-wide plotting can dispatch statistical figures.
            requires_statistics = True
        if requires_statistics:
            from ..statistics_manifest import read_statistics_manifest
            manifest = root / "project_summary/statistics/manifest.json"
            facts = read_statistics_manifest(manifest)
            if facts.state == "missing":
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
            if not facts.complete:
                raise RuntimeError(
                    f"Plotting requires completed statistics for selected statistical figures: "
                    f"{manifest} has status {facts.status!r}. Finish statistics first; "
                    "existing plots were preserved."
                )
        else:
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
        project_summary = root / "project_summary"
        if phase != "plotting" and project_summary.exists():
            archived=_archive_previous_output(project_summary,root/".archives"/"project_summary","previous")
            print(f"Archived previous project summary: {archived}")
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
        def timed(stage, fn, *args):
            started = perf_counter()
            try:
                # Statistics reports its own active family and heartbeat.
                with stage_progress(f"Finalize {stage}", enabled=verbose and phase == "statistics"
                                    and stage != "statistics"):
                    return fn(*args)
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
            timed("statistics", analyze_project, cfg, root, system_repouts)
            timed("shared_states", publish_shared_states, shared_generation, root)
            timed("diffusion_summary", summarize_project_diffusion, root, system_repouts)
            timed("multi_origin_summary", summarize_project_multi_origin, root, system_repouts, cfg)
            if cfg.get("statistics", {}).get("enabled", True):
                from ..reporting.report import generate_report
                report = timed("findings", generate_report, project_summary / "statistics")
                print(f"FINDINGS: {report}", flush=True)
            return root

        from ..reporting.figure_inventory import figure_snapshot, record_generated_figures
        from ..reporting.plotting import plot_all
        from ..reporting.membrane_figures import (
            plot_lipid_enrichment_matrix, plot_lipid_contacts, plot_project_diffusion,
        )
        from ..reporting.sdc4_figures import plot_sdc4_mechanisms
        from ..reporting.solvation_figures import plot_solvation

        figure_before = figure_snapshot(root, system_repouts)
        timed("standard_plots", plot_all, cfg, root, system_repouts)
        timed("lipid_enrichment_plots", plot_lipid_enrichment_matrix, cfg, root)
        timed("sdc4_plots", plot_sdc4_mechanisms, cfg, root, system_repouts)
        timed("lipid_contact_plots", plot_lipid_contacts, cfg, root, system_repouts)
        timed("diffusion_plot", plot_project_diffusion, cfg, root)
        timed("solvation_plots", plot_solvation, cfg, root, system_repouts)
        record_generated_figures(cfg, root, system_repouts, figure_before)
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

