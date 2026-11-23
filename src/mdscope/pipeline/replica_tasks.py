"""Execution-only specialized replica tasks and saved scalar collection.

Output paths stay authoritative in mdscope.contracts.outputs and are also
reexported by the registry. The read-only planner does not import this module
or its executable analysis implementations.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..config import analysis_enabled
from ..contracts.outputs import FUSED_OUTPUTS
from ..observables.structure_profiles import run_structure_profiles
from ..systems.sdc4.potassium_bridge_audit import run_potassium_bridge_audit


@dataclass(frozen=True)
class AnalysisTask:
    key: str
    run: Callable
    # Registered analyses this task publishes; the registry owns their paths.
    produces: tuple[str, ...]


def specialized_tasks(cfg):
    """Only work not scheduled by fused observers, retaining legacy stage keys.

    Non-serial DSSP runs through the structural-profile worker; when DSSP uses
    that backend, SASA is also excluded from fused extensions and is handled by
    the same profiles task. Potassium bridge audit remains an independent task.
    """
    dssp = cfg.get("analyses", {}).get("dssp", {})
    if dssp.get("enabled", False) and dssp.get("backend", "serial") != "serial":
        yield AnalysisTask("profiles", run_structure_profiles, ("dssp", "sasa"))
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        yield AnalysisTask("potassium_bridge_audit", run_potassium_bridge_audit,
                           ("potassium_bridge_audit",))


def additional_time_series(repout: Path, cfg):
    """Read optional scalar products without constructing executable tasks.

    The ordering is a persisted summary contract: local environment, lipid
    footprint, bilayer normal, SDC4 and lipid contacts.
    """
    output = []
    for name in (
        "local_environment", "lipid_enrichment", "bilayer_normal",
        "sdc4_mechanisms", "lipid_contacts",
    ):
        if not analysis_enabled(cfg, name):
            continue
        path = Path(repout) / FUSED_OUTPUTS[name]
        if path.exists():
            output.append(pd.read_parquet(path))
    return output
