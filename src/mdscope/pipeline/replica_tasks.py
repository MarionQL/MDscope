"""Execution-only specialized replica tasks and saved scalar collection.

Output paths stay authoritative in mdscope.registry; the read-only planner
does not import this module or the callable analysis implementations.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd

from ..config import analysis_enabled, analysis_runs
from ..observables.residue_membrane_distance import run_residue_membrane_distance
from ..observables.membrane_order import run_membrane_order
from ..observables.structure_profiles import run_structure_profiles
from ..interactions.hydrogen_bonds import run_hydrogen_bonds
from ..mechanisms.local_environment import run_local_environment
from ..mechanisms.potassium_bridge_audit import run_potassium_bridge_audit
from ..observables.bilayer_normal import run_bilayer_normal
from ..interactions.lipid_hydrogen_bonds import run_lipid_hydrogen_bonds


@dataclass(frozen=True)
class AnalysisTask:
    key: str
    run: Callable

def specialized_tasks(cfg):
    """Yield uniquely checkpointed tasks; avoid re-reading trajectory unnecessarily."""
    for run in analysis_runs(cfg,"residue_membrane_distance"):
        yield AnalysisTask("residue-membrane-"+run.run_id,
            lambda rep,cfg,out,run=run:run_residue_membrane_distance(rep,cfg,out,run.spec,run.run_id))
    if cfg["analyses"]["membrane_order"].get("enabled",False):
        yield AnalysisTask("membrane-order",run_membrane_order)
    if cfg["analyses"]["dssp"].get("enabled",False) or cfg["analyses"]["sasa"].get("enabled",False):
        yield AnalysisTask("profiles",run_structure_profiles)
    if cfg["analyses"]["hydrogen_bonds"].get("enabled",False):
        yield AnalysisTask("hydrogen-bonds",run_hydrogen_bonds)
    if cfg.get("analyses",{}).get("local_environment",{}).get("enabled",False):
        yield AnalysisTask("local_environment", run_local_environment)

    if analysis_enabled(cfg, "potassium_bridge_audit"):
        yield AnalysisTask("potassium_bridge_audit", run_potassium_bridge_audit)

    for name,fn in (("bilayer_normal",run_bilayer_normal),
                    ("lipid_hydrogen_bonds",run_lipid_hydrogen_bonds)):
        if cfg.get("analyses",{}).get(name,{}).get("enabled",False):
            yield AnalysisTask(name, fn)

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

