"""
Configuration defaults, YAML parsing, labels, paths, and safe output
helpers.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

DEFAULT_INSTANCE_ID = "default"
MULTI_RUN_ANALYSES = {"density_1d", "residue_membrane_distance"}

DEFAULTS = {
    "project": {"name": "mdscope_project", "output": "mdscope_results"},
    "progress": {"enabled": True},
    "trajectory": {"start": None, "stop": None, "stride": 1, "chunk_frames": 1000},
    "time": {"dt_ps": None, "origin_ps": 0.0},
    "selections": {
        "protein": "protein",
        "alignment": "protein and name CA",
        "features": "protein and name CA",
        "ligand": None,
        "groups": {},
    },
    "analyses": {
        "alignment": {"enabled": True, "selection": None},
        "rmsd": {"enabled": True},
        "rmsf": {"enabled": True},
        "radius_of_gyration": {"enabled": True, "selection": None},
        "density_1d": {
            "enabled": False,
            "axis": "z",
            "bin_width_A": 1.0,
            "center_selection": None,
            "selections": {},
        },
        "membrane": {
            "enabled": False,
            "axis": "z",
            "headgroup_selection": None,
            "area_per_lipid": {"enabled": True},
            "thickness": {"enabled": True},
            "msd":{
                "enabled": False,
                "selection": None,
                "include_3d": False,
                "slope_start_fraction": 0.2,
                "slope_end_fraction": 0.8
            },
        },
        "membrane_order": {
            "enabled": False,
            "axis": "z",
            "lipid_selection": None,
            "exclude_resnames": ["CHOL", "CHL1"],
            "hydrogen_cutoff_A": 1.25,
            "min_carbons_per_tail": 4,
            "per_lipid_type": True,
            "per_leaflet": True,
            "leaflet_anchor_selection": None,
            "tails": {
                "sn1": {"carbon_pattern": r"^C2(\d+)$"},
                "sn2": {"carbon_pattern": r"^C3(\d+)$"},
            },
            "lipid_overrides": {},
        },
        "residue_membrane_distance": {
            "enabled": False,
            "selection": None,
            "lipid_resnames":[
                "POPC",
                "POPE",
                "POPI",
                "DOPC",
                "DPPC",
                "PSM",
                "CHL1",
                "BGAL",
                "GM1a",
                "CER1"
            ],
            "membrane_selection": None,
            "exclude_hydrogens": True,
            "expand_residue_atoms": True,
            "initial_search_cutoff_A": 6.0,
            "search_growth_factor": 2.0,
            "fallback_block_atoms": 4096,
            "contact_cutoff_A": 5.0,
            "heatmap_max_A": None,
            "groups": False,
        },
        "dssp": {
            "enabled": True,
            "selection": None,
            "guess_hydrogens": True,
            "backend": "serial",
            "n_workers": None,
            "n_parts": None,
        },
        "bilayer_normal": {"enabled": False, "axis": "z", "species_selections": {},
                           "neighbor_radius_A": 15.0, "min_neighbors": 5,
                           "save_lipid_angles": False, "angle_chunk_frames": 100},
        "lipid_enrichment": {"enabled": False, "axis": "z", "species_selections": {},
                             "protein_selection": None,
                             "inner_radius_A": 0.0, "outer_radius_A": 15.0},
        "lipid_hydrogen_bonds": {"enabled": False, "lipid_selection": None,
                                 "donors_sel": None, "hydrogens_sel": None, "acceptors_sel": None,
                                 "d_h_cutoff_A": 1.2, "d_a_cutoff_A": 3.0,
                                 "angle_cutoff_deg": 150.0,"exclude_same_residue": True,
                                 "update_selections": False},
        "salt_bridges": {"enabled": True, "cutoff_A": 4.0},
        "ligand_contacts": {"enabled": False, "cutoff_A": 4.0, "save_atom_events": True},
        "hydrogen_bonds": {
            "enabled": False,
            "mode": "protein_protein",
            "d_h_cutoff_A": 1.2,
            "d_a_cutoff_A": 3.0,
            "angle_cutoff_deg": 150.0,
            "update_selections": False,
        },
        "sasa": {"enabled": False},
        "lifetimes": {"enabled": True, "gap_factor": 1.5},
    },
    "statistics": {
        "enabled": True,
        "bootstrap_draws": 2000,
        "ci": 0.95,
        "seed": 2026,
        "block": {
            "method": "moving",
            "length": "auto",
            "tau_multiplier": 2.0,
            "max_acf_lag": None,
        },
        "hierarchy": {
            "block_within_replicate": True,
            "resample_replicates": True,
        },
        "hypothesis_tests": {
            "enabled": True,
            "design": "independent",
            "method": "permutation",
            "permutation_draws": 10000,
            "exact_max_assignments": 100000,
            "multiple_testing": "fdr_bh",
        },
    },
    "cache": {"aligned_coordinates": True, "dtype": "float32", "keep": True},
    "plots": {
        "enabled": True,
        "min_occupancy": 0.05,
        "time_series": True,
        "distributions": True,
        "autocorrelation": True,
        "bootstrap_tables": True,
        "profiles": True,
        "interactions": True,
        "replicate_lightness_range": [0.28, 0.72],
        "font_sizes": {"base": 13, "title": 16, "suptitle": 18, "legend": 10, "ticks": 11},
        "dpi": 250,
    },
}


def merge(base: dict, update: dict) -> dict:
    out = deepcopy(base)
    for key, value in update.items():
        out[key] = merge(
            out[key],
            value
        ) if isinstance(value, dict) and isinstance(out.get(key), dict) else deepcopy(value)
    return out


@dataclass(frozen=True)
class AnalysisRun:
    name: str
    run_id: str
    label: str
    spec: dict[str, Any]


@dataclass(frozen=True)
class Replicate:
    system: str
    index: int
    topology: Path
    trajectory: Path
    label: str


def _validate_path_component(value: str, field: str) -> str:
    value = str(value)
    if not value or value in {".", ".."} or "/" in value or "\\" in value:
        raise ValueError(
            f"{field} must be a non-empty filesystem-safe name without / or \\"
        )
    return value


def analysis_runs(cfg: dict, name: str) -> list[AnalysisRun]:
    raw = cfg.get("analyses", {}).get(name, {})
    raw = {"enabled": raw} if isinstance(raw, bool) else raw
    if not isinstance(raw, dict):
        raise ValueError(f"analyses.{name} must be a mapping or boolean")
    runs = raw.get("runs")
    base = {k: deepcopy(v) for k, v in raw.items() if k != "runs"}
    if runs is None:
        return [
            AnalysisRun(
                name, DEFAULT_INSTANCE_ID, str(
                    base.get("label", name.replace("_", " ").title())
                ), base
            )
        ] if base.get("enabled", False) else []
    if name not in MULTI_RUN_ANALYSES:
        supported = ", ".join(sorted(MULTI_RUN_ANALYSES))
        raise ValueError(
            f"analyses.{name}.runs is not supported; named runs are only used "
            f"for {supported}"
        )
    if not isinstance(runs, dict) or not runs:
        raise ValueError(f"analyses.{name}.runs must be a non-empty mapping")
    enabled = bool(base.get("enabled", True))
    out = []
    for run_id, run_spec in runs.items():
        run_id = _validate_path_component(
            str(run_id),
            f"analyses.{name}.runs run ID",
        )
        if run_id == DEFAULT_INSTANCE_ID:
            raise ValueError(
                f"{DEFAULT_INSTANCE_ID!r} is reserved and cannot be used as a "
                "named run"
            )
        run_spec = {"enabled": run_spec} if isinstance(run_spec, bool) else run_spec
        if not isinstance(run_spec, dict):
            raise ValueError(f"analyses.{name}.runs.{run_id} must be a mapping or boolean")
        spec = merge(base, run_spec)
        spec["enabled"] = enabled and bool(run_spec.get("enabled", True))
        spec.pop("runs", None)
        if spec["enabled"]:
            out.append(AnalysisRun(
                name,
                run_id,
                str(run_spec.get("label", run_id.replace("_", " ").title())),
                spec
            ))
    return out


def analysis_enabled(cfg: dict, name: str) -> bool:
    return bool(analysis_runs(cfg, name))


def iter_replicates(cfg: dict, config_path: str | Path):
    base = Path(config_path).resolve().parent
    for system, spec in cfg["systems"].items():
        top = Path(spec["topology"])
        top = top if top.is_absolute() else (base / top).resolve()
        for i, rep in enumerate(spec["replicates"], 1):
            if isinstance(rep, dict):
                traj = Path(rep["trajectory"])
                label = str(rep.get("label", f"rep{i}"))
            else:
                traj = Path(rep)
                label = f"rep{i}"
            traj = traj if traj.is_absolute() else (base / traj).resolve()
            yield Replicate(system, i, top, traj, label)


def analysis_groups(cfg: dict) -> dict[str, dict]:
    out = {}
    for i, (group_id, raw) in enumerate(
        (cfg.get("selections", {}).get("groups", {}) or {}).items()
    ):
        spec = {"selection": raw} if isinstance(raw, str) else raw
        if not isinstance(spec, dict):
            raise ValueError(
                f"selections.groups.{group_id} must be a selection string "
                "or mapping"
            )
        group_id = _validate_path_component(
            str(group_id),
            "selections.groups group ID",
        )
        out[group_id] = {
            "selection": spec.get("selection"),
            "alignment_selection": spec.get("alignment_selection"),
            "label": str(spec.get("label", group_id)),
            "linestyle": spec.get("linestyle", ["-", "--", "-.", ":"][i % 4]),
        }
    return out


def system_label(cfg: dict, system: str) -> str:
    return str(cfg.get("systems", {}).get(system, {}).get("label", system))


def replicate_specs(cfg: dict, system: str) -> list[dict]:
    out = []
    for i, rep in enumerate(cfg.get("systems", {}).get(system, {}).get("replicates", []), 1):
        if isinstance(rep, dict):
            out.append({
                "id": str(rep.get("label", f"rep{i}")),
                "label": str(
                    rep.get("plot_label", rep.get("label", f"Replicate {i}"))
                ),
            })
        else:
            out.append({"id": f"rep{i}", "label": f"Replicate {i}"})
    return out


def replicate_label(cfg: dict, system: str, rep_id: str) -> str:
    return next(
        (x["label"] for x in replicate_specs(cfg, system) if x["id"] == str(rep_id)),
        str(rep_id)
    )


def group_label(cfg: dict, group_id: str) -> str:
    return analysis_groups(cfg).get(str(group_id), {}).get("label", str(group_id))


def profile_dir(repout: Path, analysis: str, run_id: str = DEFAULT_INSTANCE_ID) -> Path:
    return Path(repout) / "profiles" / analysis / run_id


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def atomic_json(path: str | Path, obj: Any) -> None:
    # Write to a sibling temporary file first so an interrupted job
    # never leaves a half-written state file.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp = tempfile.mkstemp(prefix=path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(obj, handle, indent=2, default=str)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_text(path: str | Path, text: str) -> None:
    # Atomic replacement matters because resolved configuration is also
    # used to audit/restart long runs.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp = tempfile.mkstemp(prefix=path.name, dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(text)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_parquet(path: str | Path, df: pd.DataFrame) -> None:
    # Use a unique temporary file so concurrent/aborted writers cannot
    # collide, then atomically replace the destination.
    path = Path(path)
    ensure_dir(path.parent)
    fd, tmp_name = tempfile.mkstemp(
        prefix=path.name + ".",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        try:
            df.to_parquet(tmp, index=False)
        except ImportError as exc:
            raise RuntimeError(
                "Parquet output requires pyarrow (recommended) or "
                "fastparquet."
            ) from exc
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def config_hash(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def environment_metadata() -> dict[str, Any]:
    import MDAnalysis as mda
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        commit = None
    return {
        "created_unix": time.time(),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "mdanalysis": mda.__version__,
        "git_commit": commit,
    }


def time_ps(ts, frame_index: int, cfg: dict) -> float:
    # When dt_ps is supplied it is interpreted as spacing between
    # analyzed frames, after trajectory.stride.
    spec = cfg.get("time", {})
    if spec.get("dt_ps") is not None:
        return float(spec.get("origin_ps", 0.0)) + frame_index * float(spec["dt_ps"])
    try:
        value = float(ts.time)
        if np.isfinite(value):
            return value
    except Exception:
        pass
    raise RuntimeError("Trajectory has no usable time. Set time.dt_ps in the YAML config.")


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path) as handle:
        raw = yaml.safe_load(handle) or {}
    cfg = merge(DEFAULTS, raw)
    raw_stats = raw.get("statistics", {}) or {}
    if "block_scale" in raw_stats and not (isinstance(
        raw_stats.get("block"),
        dict
    ) and "tau_multiplier" in raw_stats["block"]):
        cfg["statistics"]["block"]["tau_multiplier"] = raw_stats["block_scale"]
    if "max_acf_lag" in raw_stats and not (isinstance(
        raw_stats.get("block"),
        dict
    ) and "max_acf_lag" in raw_stats["block"]):
        cfg["statistics"]["block"]["max_acf_lag"] = raw_stats["max_acf_lag"]
    for name, raw_spec in (raw.get("analyses", {}) or {}).items():
        if isinstance(raw_spec, dict) and "runs" in raw_spec and "enabled" not in raw_spec:
            cfg.setdefault("analyses", {}).setdefault(name, {})["enabled"] = True
    for name, value in list(cfg.get("analyses", {}).items()):
        if isinstance(value, bool):
            cfg["analyses"][name] = merge(DEFAULTS["analyses"].get(name, {}), {"enabled": value})
    # Removed in 0.8: prose report and the parallel
    # distribution-statistics framework.
    # Preserve the useful part of the old autocorrelation block as a
    # plotting/max-lag alias.
    old_acf = cfg.get("analyses", {}).pop("autocorrelation", None)
    if isinstance(old_acf, dict):
        if old_acf.get("max_lag") is not None and cfg.get(
            "statistics",
            {}
        ).get("block", {}).get("max_acf_lag") is None:
            cfg.setdefault(
                "statistics",
                {}
            ).setdefault("block", {})["max_acf_lag"] = old_acf["max_lag"]
        if old_acf.get("enabled") is False:
            cfg.setdefault("plots", {})["autocorrelation"] = False
    cfg.pop("report", None)
    cfg.get("analyses", {}).pop("distribution_statistics", None)
    validate_config(cfg, Path(path).resolve().parent)
    return cfg


def validate_config(cfg: dict, base: Path = Path(".")) -> None:
    systems = cfg.get("systems")
    if not systems:
        raise ValueError("Config must define at least one system under systems:")
    if int(cfg["trajectory"]["stride"]) < 1:
        raise ValueError("trajectory.stride must be >= 1")
    if int(cfg["trajectory"].get("chunk_frames", 1000)) < 1:
        raise ValueError("trajectory.chunk_frames must be >= 1")
    time_spec = cfg.get("time", {}) or {}
    if time_spec.get("dt_ps") is not None:
        dt_ps = float(time_spec["dt_ps"])
        if not np.isfinite(dt_ps) or dt_ps <= 0:
            raise ValueError("time.dt_ps must be a finite value > 0")
    origin_ps = float(time_spec.get("origin_ps", 0.0))
    if not np.isfinite(origin_ps):
        raise ValueError("time.origin_ps must be finite")
    for name, spec in systems.items():
        if not isinstance(name, str):
            raise ValueError("System names must be strings")
        _validate_path_component(name, "system name")
        if not isinstance(spec, dict):
            raise ValueError(f"System {name!r} must be a mapping")
        if not spec.get("topology") or not spec.get("replicates"):
            raise ValueError(
                f"System {name!r} must define topology and at least one replicate"
            )
        top = Path(spec["topology"])
        top = top if top.is_absolute() else (base / top).resolve()
        if not top.exists():
            raise FileNotFoundError(top)
        labels = []
        for index, rep in enumerate(spec["replicates"], 1):
            if isinstance(rep, dict):
                if not rep.get("trajectory"):
                    raise ValueError(
                        f"System {name!r} replicate {index} is missing trajectory"
                    )
                trajectory_value = rep["trajectory"]
                label = str(rep.get("label", f"rep{index}"))
            else:
                trajectory_value = rep
                label = f"rep{index}"
            _validate_path_component(
                label,
                f"System {name!r} replicate {index} label",
            )
            labels.append(label)
            traj = Path(trajectory_value)
            traj = traj if traj.is_absolute() else (base / traj).resolve()
            if not traj.exists():
                raise FileNotFoundError(traj)
        duplicates = sorted({label for label in labels if labels.count(label) > 1})
        if duplicates:
            raise ValueError(
                f"System {name!r} has duplicate replicate labels: "
                + ", ".join(duplicates)
            )

    align = cfg.get("analyses", {}).get("alignment", {})
    if align.get(
        "enabled",
        True
    ) and not (align.get("selection") or cfg.get("selections", {}).get("alignment")):
        raise ValueError("alignment is enabled but no alignment selection was provided")
    for group_id, spec in analysis_groups(cfg).items():
        if not spec.get("selection"):
            raise ValueError(f"selections.groups.{group_id} is missing selection")

    for run in analysis_runs(cfg, "density_1d"):
        if str(run.spec.get("axis", "z")).lower() not in {"x", "y", "z"}:
            raise ValueError(
                f"density_1d[{run.run_id}].axis must be x, y, or z"
            )
        selections = run.spec.get("selections") or {}
        if float(run.spec.get("bin_width_A", 1.0)) <= 0 or not selections:
            raise ValueError(
                f"density_1d[{run.run_id}] requires positive bin_width_A and "
                "non-empty selections"
            )
        if not isinstance(selections, dict):
            raise ValueError(
                f"density_1d[{run.run_id}].selections must be a mapping"
            )
        for selection_id, selection_spec in selections.items():
            _validate_path_component(
                str(selection_id),
                f"density_1d[{run.run_id}] selection ID",
            )
            if isinstance(selection_spec, str):
                continue
            if not isinstance(selection_spec, dict):
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id} "
                    "must be a selection string or mapping"
                )
            if not selection_spec.get("selection"):
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id} "
                    "is missing selection"
                )
            mode = str(selection_spec.get("mode", "number")).lower()
            if mode not in {"number", "mass"}:
                raise ValueError(
                    f"density_1d[{run.run_id}].selections.{selection_id}.mode "
                    "must be number or mass"
                )

    mem = cfg.get("analyses", {}).get("membrane", {})
    if mem.get("enabled"):
        if str(mem.get("axis", "z")).lower() not in {
            "x",
            "y",
            "z"
        } or not mem.get("headgroup_selection"):
            raise ValueError(
                "enabled membrane analysis requires axis x/y/z and headgroup_selection"
            )
        frac = float((mem.get("msd", {}) or {}).get("slope_start_fraction", 0.2))
        end_frac = float((mem.get("msd", {}) or {}).get("slope_end_fraction", 0.8))
        if not 0 <= frac < end_frac <= 1:
            raise ValueError("membrane.msd slope fractions require 0 <= start < end <= 1")

    order = cfg.get("analyses", {}).get("membrane_order", {})
    if order.get("enabled"):
        if str(order.get("axis", "z")).lower() not in {"x", "y", "z"}:
            raise ValueError("membrane_order.axis must be x, y, or z")
        if float(order.get("hydrogen_cutoff_A", 1.25)) <= 0:
            raise ValueError("membrane_order.hydrogen_cutoff_A must be > 0")
        if int(order.get("min_carbons_per_tail", 4)) < 1:
            raise ValueError("membrane_order.min_carbons_per_tail must be >= 1")
        tails = order.get("tails", {}) or {}
        if not tails:
            raise ValueError("membrane_order.tails must define at least one tail")
        for tail, tail_spec in tails.items():
            if not isinstance(tail_spec, dict) or not (
                tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names")
            ):
                raise ValueError(
                    f"membrane_order.tails.{tail} requires carbon_pattern or carbon_names"
                )
            if tail_spec.get("carbon_names") and tail_spec.get("positions") is not None:
                if len(tail_spec["carbon_names"]) != len(tail_spec["positions"]):
                    raise ValueError(
                        f"membrane_order.tails.{tail}.positions must match carbon_names"
                    )

    if cfg["analyses"]["ligand_contacts"].get("enabled") and not cfg["selections"].get("ligand"):
        raise ValueError("ligand_contacts is enabled but selections.ligand is null")
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        if float(run.spec.get("initial_search_cutoff_A", 6.0)) <= 0:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "initial_search_cutoff_A must be > 0"
            )
        if float(run.spec.get("search_growth_factor", 2.0)) <= 1:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "search_growth_factor must be > 1"
            )
        if int(run.spec.get("fallback_block_atoms", 4096)) < 1:
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}]."
                "fallback_block_atoms must be >= 1"
            )
        if not (run.spec.get("membrane_selection") or run.spec.get("lipid_resnames")):
            raise ValueError(
                f"residue_membrane_distance[{run.run_id}] requires "
                "membrane_selection or lipid_resnames"
            )

    for analysis_name in ("lipid_enrichment","bilayer_normal"):
        ext = cfg["analyses"].get(analysis_name, {})
        if ext.get("enabled"):
            if not ext.get("species_selections") or not isinstance(ext.get("species_selections"), dict):
                raise ValueError(f"{analysis_name}.species_selections must define lipid species anchor selections")
            if str(ext.get("axis","z")).lower() not in {"x","y","z"}:
                raise ValueError(f"{analysis_name}.axis must be x, y, or z")
    enrich = cfg["analyses"].get("lipid_enrichment",{})
    if enrich.get("enabled") and not (0<=float(enrich.get("inner_radius_A",0))<float(enrich.get("outer_radius_A",15))):
        raise ValueError("lipid_enrichment radii require 0 <= inner < outer")
    normal = cfg["analyses"].get("bilayer_normal",{})
    if normal.get("enabled") and (float(normal.get("neighbor_radius_A",15))<=0 or int(normal.get("min_neighbors",5))<3):
        raise ValueError("bilayer_normal requires positive radius and >=3 neighbors")
    lipid_hb = cfg["analyses"].get("lipid_hydrogen_bonds",{})
    if lipid_hb.get("enabled") and not lipid_hb.get("lipid_selection"):
        raise ValueError("lipid_hydrogen_bonds.lipid_selection is required")

    stats = cfg.get("statistics", {})
    if int(stats.get("bootstrap_draws", 2000)) < 100:
        raise ValueError("statistics.bootstrap_draws must be >= 100")
    if not 0 < float(stats.get("ci", 0.95)) < 1:
        raise ValueError("statistics.ci must be between 0 and 1")
    block = stats.get("block", {}) or {}
    if str(block.get("method", "moving")) != "moving":
        raise ValueError("statistics.block.method currently supports only 'moving'")
    length = block.get("length", "auto")
    if length != "auto" and int(length) < 1:
        raise ValueError("statistics.block.length must be 'auto' or an integer >= 1")
    if float(block.get("tau_multiplier", 2.0)) <= 0:
        raise ValueError("statistics.block.tau_multiplier must be > 0")
    max_acf_lag = block.get("max_acf_lag")
    if max_acf_lag is not None and int(max_acf_lag) < 1:
        raise ValueError("statistics.block.max_acf_lag must be null or >= 1")
    tests = stats.get("hypothesis_tests", {}) or {}
    if tests.get("enabled", True):
        design = str(tests.get("design", "independent"))
        if design not in {"independent", "paired"}:
            raise ValueError(
                "statistics.hypothesis_tests.design must be independent or "
                "paired"
            )
        if design == "paired" and len(cfg.get("systems", {})) > 1:
            label_sets = []
            for system in cfg["systems"]:
                labels = {
                    spec["id"]
                    for spec in replicate_specs(cfg, system)
                }
                label_sets.append((system, labels))
            reference_system, reference_labels = label_sets[0]
            for system, labels in label_sets[1:]:
                if labels != reference_labels:
                    raise ValueError(
                        "paired hypothesis tests require identical replicate "
                        f"labels across systems; {reference_system!r} and "
                        f"{system!r} differ"
                    )
        if str(tests.get("method", "permutation")) != "permutation":
            raise ValueError(
                "statistics.hypothesis_tests.method currently supports only permutation"
            )
        if int(tests.get("permutation_draws", 10000)) < 100:
            raise ValueError(
                "statistics.hypothesis_tests.permutation_draws must be >= 100"
            )
        if int(tests.get("exact_max_assignments", 100000)) < 1:
            raise ValueError(
                "statistics.hypothesis_tests.exact_max_assignments must be >= 1"
            )
        if str(tests.get("multiple_testing", "fdr_bh")) not in {"fdr_bh", "none"}:
            raise ValueError("statistics.hypothesis_tests.multiple_testing must be fdr_bh or none")
