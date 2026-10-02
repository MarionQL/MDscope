"""
Project orchestration: run replicates, consolidate outputs, then compute
statistics and figures.
"""

from __future__ import annotations

import glob
import json
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
from .plotting import plot_all
from .profiles import run_residue_membrane_distance
from .profiles import postprocess_structure_profiles, run_structure_profiles
from .statistics import analyze_project
from .timeseries import run_streaming_pass


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
    return path if path.is_absolute() else Path(config_path).resolve().parent / path


def initialize(config_path: str | Path, cfg: dict):
    root = ensure_dir(project_root(cfg, config_path))
    hash_value = config_hash({"config": cfg, "mdscope_version": __version__})
    atomic_text(root / "resolved_config.yaml", yaml.safe_dump(cfg, sort_keys=False))
    metadata = environment_metadata()
    metadata.update({
        "config_hash": hash_value,
        "config_path": str(Path(config_path).resolve()),
        "mdscope_version": __version__
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
        duplicate_metrics =[
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
    time_sources = [base] + profile_ts
    available_times =[
        table[
            "time_ps"
        ].to_numpy(
            float
        )for table in time_sources if table is not None and not table.empty and "time_ps" in table
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
    if not grouped.empty:
        atomic_parquet(summary / "group_timeseries.parquet", grouped)


def run_replicate(rep, cfg: dict, root: Path, state: RunState):
    prefix = f"{rep.system}-{rep.label}"
    repout = ensure_dir(root / rep.system / rep.label)
    streaming_analyses =[
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
    )for name in streaming_analyses) or bool(cfg.get(
        "cache",
        {}
    ).get("aligned_coordinates", False))
    if streaming_required and not state.complete(prefix + "-stream"):
        state.mark(prefix + "-stream", run_streaming_pass(rep, cfg, root, state.config_hash))

    for run in analysis_runs(cfg, "residue_membrane_distance"):
        key = f"{prefix}-residue-membrane-{run.run_id}"
        if not state.complete(key):
            state.mark(key, run_residue_membrane_distance(rep, cfg, repout, run.spec, run.run_id))

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
    if not state.complete(prefix + "-summary"):
        summarize_replicate(repout, cfg)
        state.mark(prefix + "-summary")


def finalize_project(config_path: str | Path, cfg=None, root=None, state=None):
    cfg = cfg or load_config(config_path)
    if root is None or state is None:
        root, state = initialize(config_path, cfg)
    replicates = list(iter_replicates(cfg, config_path))
    system_repouts ={
        system:[
            root / rep.system / rep.label for rep in replicates if rep.system == system
        ]for system in cfg["systems"]
    }
    # Finalization only consumes persisted tables/arrays: no
    # trajectories are reopened.
    analyze_project(cfg, root, system_repouts)
    plot_all(cfg, root, system_repouts)
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
