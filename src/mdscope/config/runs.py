"""Analysis instances, replica identity, labels and sampled trajectory time."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .defaults import DEFAULT_INSTANCE_ID, MULTI_RUN_ANALYSES, merge

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
    msd_trajectory: Path | None = None


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


def sdc4_enabled(cfg: dict) -> bool:
    """SDC4 feature activation without importing the execution pipeline."""
    return bool(cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False))


def project_root(cfg: dict, config_path: str | Path) -> Path:
    """Resolve project.output relative to the YAML, not the working directory."""
    path = Path(cfg["project"]["output"])
    return (path if path.is_absolute() else Path(config_path).resolve().parent / path).resolve()


def iter_replicates(cfg: dict, config_path: str | Path):
    base = Path(config_path).resolve().parent
    for system, spec in cfg["systems"].items():
        top = Path(spec["topology"])
        top = top if top.is_absolute() else (base / top).resolve()
        for i, rep in enumerate(spec["replicates"], 1):
            if isinstance(rep, dict):
                traj = Path(rep["trajectory"])
                msd_source = rep.get("msd_trajectory")
                label = str(rep.get("label", f"rep{i}"))
            else:
                traj = Path(rep)
                msd_source = None
                label = f"rep{i}"
            traj = traj if traj.is_absolute() else (base / traj).resolve()
            msd_path = Path(msd_source) if msd_source else None
            if msd_path is not None and not msd_path.is_absolute():
                msd_path = (base / msd_path).resolve()
            yield Replicate(system, i, top, traj, label, msd_path)


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


