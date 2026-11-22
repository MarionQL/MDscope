"""YAML loading and public configuration API, with domain validation separated."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from ..contracts.metric_selection import STATISTICS_FAMILIES
from ..errors import ConfigError, InputPathError
from .defaults import DEFAULTS, DEFAULT_INSTANCE_ID, MULTI_RUN_ANALYSES, TAIL_PRESETS, merge
from .runs import (
    AnalysisRun, Replicate, _validate_path_component,
    analysis_runs, analysis_enabled, analysis_groups,
    group_label, iter_replicates, profile_dir, project_root,
    replicate_label, replicate_specs, system_label, time_ps,
)
from .validation import validate_config, residue_response_settings

def load_config(path: str | Path, *, check_paths: bool = True) -> dict[str, Any]:
    """Load a user YAML file; any parse or validation failure is a ``ConfigError``."""
    try:
        with open(path) as handle:
            raw = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:
        raise InputPathError(f"Configuration file does not exist: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path}: invalid YAML: {exc}") from exc
    try:
        return resolve_config(raw, base=Path(path).resolve().parent, check_paths=check_paths)
    except ValueError as exc:
        raise ConfigError(f"{path}: {exc}") from exc


def resolve_config(raw: dict, *, base: Path = Path("."),
                   check_paths: bool = False) -> dict[str, Any]:
    """Resolve an in-memory YAML mapping through the same CLI validation path."""
    if not isinstance(raw, dict):
        raise ValueError("Configuration YAML must be a mapping")
    for section in ("project", "trajectory", "time", "selections", "analyses", "statistics", "plots"):
        if section in raw and not isinstance(raw[section], dict):
            raise ValueError(f"{section} must be a mapping")
    cfg = merge(DEFAULTS, raw)
    if "report" in raw:
        raise ValueError("Remove the retired report section. Statistics now write findings.md.")
    if {"block_scale", "max_acf_lag"} & raw.get("statistics", {}).keys():
        raise ValueError("Use statistics.block.tau_multiplier and statistics.block.max_acf_lag.")
    if {"autocorrelation", "distribution_statistics"} & raw.get("analyses", {}).keys():
        raise ValueError("Remove retired analyses.autocorrelation or analyses.distribution_statistics.")
    from ..config_schema import validate_known_keys
    validate_known_keys(raw, DEFAULTS)
    for name, raw_spec in (raw.get("analyses", {}) or {}).items():
        if isinstance(raw_spec, dict) and "runs" in raw_spec and "enabled" not in raw_spec:
            cfg.setdefault("analyses", {}).setdefault(name, {})["enabled"] = True
    normalize_analysis_booleans(cfg)
    _apply_tail_preset(cfg, raw)
    validate_config(cfg, Path(base).resolve(), check_paths=check_paths)
    return cfg


def _apply_tail_preset(cfg, raw):
    """Expand ``membrane_order.tail_preset`` once; explicit tails replace preset tails.

    Without a preset, user tails keep merging into the default tails as before.
    The preset name stays in the resolved configuration as provenance.
    """
    order = cfg["analyses"]["membrane_order"]
    name = order.get("tail_preset")
    if name is None:
        return
    if name not in TAIL_PRESETS:
        raise ValueError(
            f"membrane_order.tail_preset must be null or one of {sorted(TAIL_PRESETS)}; got {name!r}"
        )
    raw_order = (raw.get("analyses") or {}).get("membrane_order") or {}
    explicit = raw_order.get("tails") or {}
    if not isinstance(explicit, dict):
        raise ValueError("membrane_order.tails must be a mapping")
    order["tails"] = {**deepcopy(TAIL_PRESETS[name]), **deepcopy(explicit)}


def normalize_analysis_booleans(cfg):
    """Normalize supported boolean shorthand for configuration and draft catalogs."""
    for name, value in list(cfg.get("analyses", {}).items()):
        if isinstance(value, bool):
            cfg["analyses"][name] = merge(DEFAULTS["analyses"].get(name, {}), {"enabled": value})
        elif not isinstance(value, dict):
            raise ValueError(f"analyses.{name} must be a mapping or boolean")

def plot_enabled(cfg, family):
    """Resolve master/family switches; skip entirely unselected families."""
    options = cfg.get("plots", {})
    if (not options.get("enabled", True)
            or not options.get("families", {}).get(family, True)
            or not options.get(family, True)):
        return False
    selected = options.get("include")
    return selected is None or any(item.startswith(f"{family}:") for item in selected)


def figure_enabled(cfg, plot_id):
    """A specific figure is eligible only when all its switches allow it."""
    if not plot_enabled(cfg, plot_id.split(":", 1)[0]):
        return False
    selected = cfg.get("plots", {}).get("include")
    return selected is None or plot_id in selected


def statistics_family_config(cfg, family):
    """Resolve one statistics family's explicit overrides without mutating cfg."""
    from copy import deepcopy
    stats = cfg.get("statistics", {})
    option = stats.get("families", {}).get(family, {})
    option = {"enabled": option} if isinstance(option, bool) else option
    if not stats.get("enabled", True) or not option.get("enabled", True):
        return None
    if family in STATISTICS_FAMILIES:
        selection = stats.get("metric_selection")
        if selection is not None and selection.get(family) == []:
            return None
    if not option and stats.get("bootstrap", True):
        return cfg
    result = deepcopy(cfg)
    for key in ("bootstrap", "bootstrap_draws", "ci", "adaptive_hierarchy"):
        if key in option:
            result["statistics"][key] = option[key]
    if "hypothesis_tests" in option:
        result["statistics"].setdefault("hypothesis_tests", {})["enabled"] = option["hypothesis_tests"]
    if "hierarchy" in option:
        result["statistics"].setdefault("hierarchy", {}).update(option["hierarchy"])
    if not result["statistics"].get("bootstrap", True):
        result["statistics"]["bootstrap_draws"] = 0
    return result

# Export the historical public config API while modules are split.
__all__ = [
    "DEFAULTS", "DEFAULT_INSTANCE_ID", "MULTI_RUN_ANALYSES", "TAIL_PRESETS",
    "AnalysisRun", "Replicate", "_validate_path_component", "analysis_runs",
    "analysis_enabled", "analysis_groups", "group_label", "iter_replicates",
    "profile_dir", "project_root", "replicate_label", "replicate_specs",
    "system_label", "time_ps", "merge", "load_config", "resolve_config",
    "normalize_analysis_booleans", "validate_config", "plot_enabled",
    "figure_enabled", "residue_response_settings", "statistics_family_config",
]
