"""
Project orchestration: run replicates, consolidate outputs, then compute
statistics and figures.
"""

from __future__ import annotations

import glob
import hashlib
import json
import shutil
from time import perf_counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import __version__
from .config import (
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
from .interactions import run_hydrogen_bonds
from .interactions import postprocess_interactions
from .membrane_order import membrane_order_timeseries, run_membrane_order
from .plotting import plot_all
from .membrane_figures import plot_lipid_enrichment_matrix
from .profiles import run_residue_membrane_distance
from .profiles import postprocess_structure_profiles, run_structure_profiles
from .statistics import analyze_project
from .streaming import run_streaming_pass
from .fused_membrane import enabled_names as fused_enabled_names, run_fused_membrane
from .fused_extensions import enabled_names as extension_enabled_names, run_fused_extensions
from .fused_distance import enabled_names as distance_enabled_names, run_fused_distances
from .fused_all import run_fused_only
from .sdc4_mechanisms import enabled as sdc4_enabled, run_sdc4_mechanisms
from .shared_states import fit_shared_states
from .sdc4_figures import plot_sdc4_mechanisms
from .report import generate_report
from .registry import specialized_tasks,additional_time_series,fused_output,fused_output_plan
from .msd import summarize_diffusion, summarize_project_diffusion, summarize_project_multi_origin
from .observables.multi_origin_pipeline import run_multi_origin_diffusion
from .coupling_screen import screen_couplings


class RunState:
    # Small completion markers let interrupted long trajectories resume
    # without repeating finished stages.
    def __init__(self, root: Path, config_hash_value: str):
        self.root = ensure_dir(root / ".state")
        self.config_hash = config_hash_value

    def marker(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def complete(self, key: str) -> bool:
        path = self.marker(key)
        if not path.exists():
            return False
        try:
            return json.loads(path.read_text()).get("config_hash") == self.config_hash
        except Exception:
            return False

    def mark(self, key: str, extra=None):
        atomic_json(
            self.marker(key),
            {"status": "complete", "config_hash": self.config_hash, **(extra or {})}
        )


def project_root(cfg: dict, config_path: str | Path) -> Path:
    path = Path(cfg["project"]["output"])
    return (
        path
        if path.is_absolute()
        else Path(config_path).resolve().parent / path
    )


def package_code_hash() -> str:
    """Hash package source so restart markers invalidate after code edits."""
    digest = hashlib.sha256()
    package_dir = Path(__file__).resolve().parent
    for path in sorted(package_dir.rglob("*.py")):
        digest.update(path.relative_to(package_dir).as_posix().encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()[:16]


def initialize(config_path: str | Path, cfg: dict):
    root = ensure_dir(project_root(cfg, config_path))
    code_hash = package_code_hash()
    hash_value = config_hash({
        "config": cfg,
        "mdscope_version": __version__,
        "code_hash": code_hash,
    })
    atomic_text(
        root / "resolved_config.yaml",
        yaml.safe_dump(cfg, sort_keys=False),
    )
    metadata = environment_metadata()
    metadata.update({
        "config_hash": hash_value,
        "code_hash": code_hash,
        "config_path": str(Path(config_path).resolve()),
        "mdscope_version": __version__,
    })
    atomic_json(root / "run_metadata.json", metadata)
    return root, RunState(root, hash_value)


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
    tables = [table.copy() for table in tables if table is not None and not table.empty]
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
    if not whole.empty:
        atomic_parquet(summary / "timeseries.parquet", whole)
        if cfg.get("analyses",{}).get("membrane",{}).get("enabled",False) and cfg.get("analyses",{}).get("membrane",{}).get("msd",{}).get("enabled",False):
            summarize_diffusion(repout,cfg)
    if not grouped.empty:
        atomic_parquet(summary / "group_timeseries.parquet", grouped)


def _replicate_state_prefix(rep) -> str:
    identity_hash = config_hash({
        "system": rep.system,
        "replicate": rep.label,
    })
    return f"rep-{identity_hash}"


def _archive_previous_output(source: Path, archive_root: Path, label: str) -> Path:
    """Archive stale results rather than silently deleting expensive MD analyses."""
    ensure_dir(archive_root)
    candidate=archive_root / label
    index=1
    while candidate.exists():
        candidate=archive_root / f"{label}-{index}"
        index+=1
    source.rename(candidate)
    return candidate


def _prepare_replicate_output(
    rep,
    root: Path,
    state: RunState,
) -> Path:
    repout = root / rep.system / rep.label
    hash_file = repout / ".config_hash"
    if repout.exists():
        previous_hash = None
        try:
            previous_hash = hash_file.read_text().strip()
        except OSError:
            pass
        if previous_hash != state.config_hash:
            prior = previous_hash or "unknown_hash"
            archived=_archive_previous_output(
                repout,root/".archives"/"replicates"/rep.system/rep.label,prior)
            print(f"Archived previous replicate output: {archived}")
    repout = ensure_dir(repout)
    atomic_text(hash_file, state.config_hash + "\n")
    return repout


def run_replicate(rep, cfg: dict, root: Path, state: RunState):
    prefix = _replicate_state_prefix(rep)
    repout = _prepare_replicate_output(rep, root, state)
    stage_timings = []
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
    if streaming_required and not state.complete(prefix + "-stream"):
        t0 = perf_counter()
        state.mark(prefix + "-stream", run_streaming_pass(rep, cfg, root, state.config_hash))
        stage_timings.append({"stage": "streaming", "seconds": perf_counter() - t0})

    # Without protein streaming, run *all* optional observer families in
    # one shared frame loop instead of three separate fallback traversals.
    if not streaming_required:
        pending_outputs = fused_output_plan(cfg, repout)
        if any(not entry.path.is_file() for entry in pending_outputs):
            t0 = perf_counter()
            run_fused_only(rep, cfg, repout)
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
        need_fallback = any(not expected[name].path.is_file() for name in fused_names)
        if need_fallback:
            t0 = perf_counter()
            run_fused_membrane(rep, cfg, repout)
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
        extension_outputs = {name: fused_output(repout, name) for name in extension_names}
        missing = [name for name in extension_names if not extension_outputs[name].path.is_file()]
        if missing:
            t0 = perf_counter()
            run_fused_extensions(rep, cfg, repout)
            stage_timings.append({"stage": "fused-extensions-fallback", "seconds": perf_counter() - t0})
        for name in extension_names:
            if not extension_outputs[name].path.is_file():
                raise RuntimeError(f"Missing fused output: {extension_outputs[name].path}")
            key = prefix + "-" + extension_outputs[name].checkpoint_key
            if key.endswith("-profiles") and not all(
                item in extension_names for item in ("dssp", "sasa")
                if cfg["analyses"][item].get("enabled", False)
            ):
                continue
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
            run_fused_distances(rep, cfg, repout)
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
        key = prefix + "-sdc4-mechanisms"
        if not state.complete(key):
            state.mark(key, {"fused_into": "streaming" if streaming_required
                             else "fused-only"})

    for task in specialized_tasks(cfg):
        if task.key.startswith("residue-membrane-") and task.key[len("residue-membrane-"):] in distance_names:
            continue
        if task.key in {"membrane-order", "bilayer_normal", "lipid_enrichment"}:
            continue
        if task.key in {"hydrogen-bonds", "lipid_hydrogen_bonds", "local_environment"} and {
            "hydrogen-bonds": "hydrogen_bonds",
            "lipid_hydrogen_bonds": "lipid_hydrogen_bonds",
            "local_environment": "local_environment",
        }[task.key] in extension_names:
            continue
        if task.key == "profiles" and all(
            name in extension_names
            for name in ("dssp", "sasa")
            if cfg["analyses"][name].get("enabled", False)
        ):
            continue
        key=prefix+"-"+task.key
        if not state.complete(key):
            t0 = perf_counter()
            result=task.run(rep,cfg,repout)
            state.mark(key,result if isinstance(result,dict) else {})
            stage_timings.append({"stage": task.key, "seconds": perf_counter() - t0})
    msd_cfg = cfg.get("analyses", {}).get("membrane", {}).get("msd", {})
    if msd_cfg.get("enabled") and msd_cfg.get("multi_origin", False):
        target = repout / "summary" / "multi_origin_diffusion.parquet"
        if not state.complete(prefix + "-multi-origin-msd") or not target.is_file():
            t0 = perf_counter()
            result = run_multi_origin_diffusion(rep, cfg, repout)
            state.mark(prefix + "-multi-origin-msd", result)
            stage_timings.append({"stage": "multi-origin-msd", "seconds": perf_counter() - t0})
    if not state.complete(prefix + "-summary"):
        t0 = perf_counter()
        summarize_replicate(repout, cfg)
        state.mark(prefix + "-summary")
        stage_timings.append({"stage": "summary", "seconds": perf_counter() - t0})
    if stage_timings:
        atomic_json(repout / "stage_timings.json", {
            "system": rep.system, "replicate": rep.label,
            "stages": stage_timings,
        })
        for item in stage_timings:
            print(f"TIMING {rep.system}/{rep.label} {item['stage']}: {item['seconds']:.3f} s", flush=True)


def finalize_project(config_path: str | Path, cfg=None, root=None, state=None):
    cfg = cfg or load_config(config_path)
    if root is None or state is None:
        root, state = initialize(config_path, cfg)
    replicates = list(iter_replicates(cfg, config_path))
    incomplete = [
        f"{rep.system}/{rep.label}"
        for rep in replicates
        if not state.complete(_replicate_state_prefix(rep) + "-summary")
    ]
    if incomplete:
        raise RuntimeError(
            "Cannot finalize before all replicates finish the current "
            "configuration/code version. Incomplete: "
            + ", ".join(incomplete)
        )
    project_summary = root / "project_summary"
    if project_summary.exists():
        archived=_archive_previous_output(project_summary,root/".archives"/"project_summary",state.config_hash)
        print(f"Archived previous project summary: {archived}")
    # Even with statistics/plots disabled, a completed project needs a\n    # summary directory for optional report generation and downstream tools.\n    ensure_dir(project_summary)\n    system_repouts = {\n        system: [
            root / rep.system / rep.label
            for rep in replicates
            if rep.system == system
        ]
        for system in cfg["systems"]
    }
    # Finalization only consumes persisted tables/arrays: no
    # trajectories are reopened.
    analyze_project(cfg, root, system_repouts)
    plot_all(cfg, root, system_repouts)
    plot_lipid_enrichment_matrix(cfg,root)
    fit_shared_states(cfg, root, system_repouts)
    plot_sdc4_mechanisms(cfg, root, system_repouts)
    summarize_project_diffusion(root,system_repouts,cfg)
    summarize_project_multi_origin(root,system_repouts)
    if cfg.get('analyses',{}).get('coupling_screen',{}).get('enabled',False):
        screen_couplings(cfg,root,system_repouts)
    report_result = generate_report(cfg, root, system_repouts)
    if report_result:
        print(f"REPORT PDF: {report_result['pdf']}", flush=True)
        if report_result.get('bundle'):
            print(f"REPORT BUNDLE: {report_result['bundle']}", flush=True)
    return root


def run_one(config_path: str | Path, index: int):
    """
    Run one zero-based systemxreplicate entry; useful for external job
    arrays.
    """
    cfg = load_config(config_path)
    replicates = list(iter_replicates(cfg, config_path))
    if index < 0 or index >= len(replicates):
        raise IndexError(f"replicate index {index} is out of range 0..{len(replicates) - 1}")
    root, state = initialize(config_path, cfg)
    rep = replicates[index]
    run_replicate(rep, cfg, root, state)
    return root / rep.system / rep.label


def run_all(config_path: str | Path):
    cfg = load_config(config_path)
    root, state = initialize(config_path, cfg)
    for rep in iter_replicates(cfg, config_path):
        run_replicate(rep, cfg, root, state)
    finalize_project(config_path, cfg, root, state)
    return root
