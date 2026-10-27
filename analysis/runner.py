"""
Project orchestration: run replicates, consolidate outputs, then compute
statistics and figures.
"""

from __future__ import annotations

import glob
import hashlib
import json
import shutil
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
from .profiles import run_residue_membrane_distance
from .profiles import postprocess_structure_profiles, run_structure_profiles
from .statistics import analyze_project
from .timeseries import run_streaming_pass
from .local_environment import run_local_environment
from .msd import summarize_diffusion
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
    for path in sorted(package_dir.glob("*.py")):
        digest.update(path.name.encode("ascii"))
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
    out = tables[0].sort_values(keys).drop_duplicates(keys, keep="last")
    for table in tables[1:]:
        table = table.sort_values(keys).drop_duplicates(keys, keep="last")
        duplicate_metrics = [
            col for col in table.columns if col in out.columns and col not in keys
        ]
        table = table.drop(columns=duplicate_metrics)
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
        base = base.sort_values("time_ps").drop_duplicates("time_ps", keep="last")
    group_base = _load_parts(raw / "group_timeseries-part-*.parquet")
    if not group_base.empty:
        group_base = group_base.sort_values(["group", "time_ps"]).drop_duplicates(
            ["group", "time_ps"],
            keep="last"
        )

    profile_ts, profile_group_ts = postprocess_structure_profiles(repout, cfg)
    time_sources = [base] + profile_ts + membrane_order_timeseries(repout)
    local_path = repout / 'extensions' / 'local_environment' / 'timeseries.parquet'
    if local_path.exists():
        time_sources.append(pd.read_parquet(local_path))
    available_times = [
        table[
            "time_ps"
        ].to_numpy(
            float
        ) for table in time_sources if table is not None and not table.empty and "time_ps" in table
    ]
    if not available_times:
        hb_counts = repout / "interactions" / "hydrogen_bonds" / "count_timeseries.parquet"
        if hb_counts.exists():
            available_times.append(pd.read_parquet(hb_counts)["time_ps"].to_numpy(float))
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
        for col in whole.columns:
            if col != "time_ps" and (col.endswith("_count") or col.endswith("_residue_count")):
                whole[col] = pd.to_numeric(whole[col], errors="coerce").fillna(0)
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
            shutil.rmtree(repout)
    repout = ensure_dir(repout)
    atomic_text(hash_file, state.config_hash + "\n")
    return repout


def run_replicate(rep, cfg: dict, root: Path, state: RunState):
    prefix = _replicate_state_prefix(rep)
    repout = _prepare_replicate_output(rep, root, state)
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
        state.mark(prefix + "-stream", run_streaming_pass(rep, cfg, root, state.config_hash))

    for run in analysis_runs(cfg, "residue_membrane_distance"):
        key = f"{prefix}-residue-membrane-{run.run_id}"
        if not state.complete(key):
            state.mark(key, run_residue_membrane_distance(rep, cfg, repout, run.spec, run.run_id))

    if cfg["analyses"]["membrane_order"].get("enabled", False):
        key = prefix + "-membrane-order"
        if not state.complete(key):
            state.mark(key, run_membrane_order(rep, cfg, repout))

    if (cfg["analyses"]["dssp"].get(
        "enabled",
        False
    ) or cfg["analyses"]["sasa"].get(
        "enabled",
        False
    )) and not state.complete(prefix + "-profiles"):
        run_structure_profiles(rep, cfg, repout)
        state.mark(prefix + "-profiles")
    if cfg["analyses"]["hydrogen_bonds"].get(
        "enabled",
        False
    ) and not state.complete(prefix + "-hydrogen-bonds"):
        run_hydrogen_bonds(rep, cfg, repout)
        state.mark(prefix + "-hydrogen-bonds")
    for name, fn in [('local_environment',run_local_environment)]:
        if cfg.get('analyses',{}).get(name,{}).get('enabled',False):
            key=prefix+'-'+name
            if not state.complete(key):
                state.mark(key,fn(rep,cfg,repout))
    if not state.complete(prefix + "-summary"):
        summarize_replicate(repout, cfg)
        state.mark(prefix + "-summary")


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
        shutil.rmtree(project_summary)
    system_repouts = {
        system: [
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
    if cfg.get('analyses',{}).get('coupling_screen',{}).get('enabled',False):
        screen_couplings(cfg,root,system_repouts)
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
