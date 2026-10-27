"""Explicit, typed registry for specialized replicate analyses.

Adding a specialized analysis requires:
1. implement run(rep,cfg,repout) in its own module;
2. add its enabled predicate and output contract here.
Time series go to the same summary/stats/plot path as existing observables.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
import pandas as pd
from .config import analysis_runs
from .residue_membrane_distance import run_residue_membrane_distance
from .membrane_order import run_membrane_order
from .structure_profiles import run_structure_profiles
from .hydrogen_bonds import run_hydrogen_bonds
from .local_environment import run_local_environment

@dataclass(frozen=True)
class AnalysisTask:
    key: str
    run: Callable
    scalar_output: str|None=None

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
        yield AnalysisTask("local_environment",run_local_environment,
                           "extensions/local_environment/timeseries.parquet")

def additional_time_series(repout:Path):
    """Canonical scalar output merger for registered extensions."""
    table=Path(repout)/"extensions"/"local_environment"/"timeseries.parquet"
    return [pd.read_parquet(table)] if table.exists() else []
