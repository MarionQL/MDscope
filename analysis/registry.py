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
from .observables.residue_membrane_distance import run_residue_membrane_distance
from .observables.membrane_order import run_membrane_order
from .observables.structure_profiles import run_structure_profiles
from .interactions.hydrogen_bonds import run_hydrogen_bonds
from .mechanisms.local_environment import run_local_environment
from .observables.lipid_enrichment import run_lipid_enrichment
from .observables.bilayer_normal import run_bilayer_normal
from .interactions.lipid_hydrogen_bonds import run_lipid_hydrogen_bonds

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

    for name,fn in (("lipid_enrichment",run_lipid_enrichment),
                    ("bilayer_normal",run_bilayer_normal),
                    ("lipid_hydrogen_bonds",run_lipid_hydrogen_bonds)):
        if cfg.get("analyses",{}).get(name,{}).get("enabled",False):
            output=(f"extensions/{name}/timeseries.parquet"
                    if name!="lipid_hydrogen_bonds" else None)
            yield AnalysisTask(name,fn,output)

def additional_time_series(repout:Path,cfg):
    """Discover scalar series from registered analysis output contracts."""
    output=[]
    for task in specialized_tasks(cfg):
        if task.scalar_output is None:continue
        path=Path(repout)/task.scalar_output
        if path.exists():
            output.append(pd.read_parquet(path))
    if cfg.get("analyses", {}).get("sdc4_mechanisms", {}).get("enabled", False):
        path = Path(repout) / "extensions/sdc4_mechanisms/timeseries.parquet"
        if path.exists():
            output.append(pd.read_parquet(path))
    return output


# Canonical persisted outputs for fused observers. Keys deliberately retain the
# existing analysis names; marker aliases preserve completed RunState histories.
FUSED_OUTPUTS = {
    "membrane_order": ("profiles/membrane_order/timeseries.parquet", "membrane-order", "membrane"),
    "bilayer_normal": ("extensions/bilayer_normal/timeseries.parquet", "bilayer_normal", "membrane"),
    "lipid_enrichment": ("extensions/lipid_enrichment/timeseries.parquet", "lipid_enrichment", "membrane"),
    "hydrogen_bonds": ("interactions/hydrogen_bonds/events.parquet", "hydrogen-bonds", "extensions"),
    "lipid_hydrogen_bonds": ("interactions/lipid_hydrogen_bonds/events.parquet", "lipid_hydrogen_bonds", "extensions"),
    "dssp": ("profiles/dssp/assignments.parquet", "profiles", "extensions"),
    "sasa": ("profiles/sasa/residue_sasa.parquet", "profiles", "extensions"),
    "local_environment": ("extensions/local_environment/timeseries.parquet", "local_environment", "extensions"),
    "sdc4_mechanisms": ("extensions/sdc4_mechanisms/timeseries.parquet", "sdc4-mechanisms", "sdc4"),
}


@dataclass(frozen=True)
class FusedOutput:
    name: str
    path: Path
    checkpoint_key: str
    family: str


def fused_output(repout: Path, name: str, *, run_id: str | None = None) -> FusedOutput:
    """The persisted output and legacy completion marker of a fused analysis."""
    if name == "residue_membrane_distance":
        if run_id is None:
            raise ValueError("residue_membrane_distance requires a run_id")
        from .config import _validate_path_component
        run_id = _validate_path_component(str(run_id), "distance run ID")
        return FusedOutput(
            name, Path(repout) / "profiles" / name / run_id / "distance_matrix_A.npy",
            "residue-membrane-" + run_id, "distance",
        )
    rel, key, family = FUSED_OUTPUTS[name]
    return FusedOutput(name, Path(repout) / rel, key, family)


def fused_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Plan enabled fused outputs without creating readers or opening trajectories."""
    # Imports stay local: registry is also used by postprocessing.
    from .pipeline.fused_membrane import enabled_names as membrane_enabled
    from .pipeline.fused_extensions import enabled_names as extensions_enabled
    from .pipeline.fused_distance import enabled_names as distance_enabled
    from .mechanisms.sdc4_mechanisms import enabled as sdc4_enabled
    entries = [fused_output(repout, name) for name in membrane_enabled(cfg)]
    entries.extend(fused_output(repout, name) for name in extensions_enabled(cfg))
    entries.extend(fused_output(repout, "residue_membrane_distance", run_id=run_id)
                   for run_id in distance_enabled(cfg))
    if sdc4_enabled(cfg):
        entries.append(fused_output(repout, "sdc4_mechanisms"))
    return tuple(entries)
