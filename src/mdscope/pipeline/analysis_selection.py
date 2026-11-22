"""Lightweight, authoritative enabled-analysis predicates used by producers and planning.

Keep scheduling decisions here so read-only output planning need not import
trajectory implementations. Production consumers use these same predicates.
"""
from __future__ import annotations

from ..config import analysis_runs
from ..contracts.outputs import EXTENSION_ANALYSES, MEMBRANE_ANALYSES

FUSED_MEMBRANE_ANALYSES = MEMBRANE_ANALYSES


def membrane_enabled_names(cfg):
    return tuple(
        name for name in FUSED_MEMBRANE_ANALYSES
        if cfg.get("analyses", {}).get(name, {}).get("enabled", False)
    )


def extensions_enabled_names(cfg):
    specs = cfg.get("analyses", {})
    names = []
    for name in EXTENSION_ANALYSES:
        if not specs.get(name, {}).get("enabled", False):
            continue
        # Parallel DSSP requires the specialized worker implementation.
        if name in {"dssp", "sasa"} and specs.get("dssp", {}).get("enabled", False) and specs["dssp"].get("backend", "serial") != "serial":
            continue
        names.append(name)
    if (specs.get("hydration", {}).get("enabled", False)
            or specs.get("sdc4_mechanisms", {}).get("hs_sugar_water_cutoff_A") is not None):
        names.append("hydration")
    return tuple(names)


def distance_enabled_names(cfg):
    return tuple(run.run_id for run in analysis_runs(cfg, "residue_membrane_distance"))


def sdc4_enabled(cfg):
    return bool(cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False))
