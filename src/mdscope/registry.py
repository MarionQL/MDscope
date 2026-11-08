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
    if cfg.get('analyses', {}).get('lipid_contacts', {}).get('enabled', False):
        path = Path(repout)/'extensions/lipid_contacts/timeseries.parquet'
        if path.exists():
            output.append(pd.read_parquet(path))
    return output


# Persisted output paths are independent of scheduling and checkpoint naming.
FUSED_OUTPUTS = {
    "lipid_contacts": "extensions/lipid_contacts/timeseries.parquet",
    "membrane_order": "profiles/membrane_order/timeseries.parquet",
    "bilayer_normal": "extensions/bilayer_normal/timeseries.parquet",
    "lipid_enrichment": "extensions/lipid_enrichment/timeseries.parquet",
    "hydrogen_bonds": "interactions/hydrogen_bonds/events.parquet",
    "lipid_hydrogen_bonds": "interactions/lipid_hydrogen_bonds/events.parquet",
    "dssp": "profiles/dssp/assignments.parquet",
    "sasa": "profiles/sasa/residue_sasa.parquet",
    "local_environment": "extensions/local_environment/timeseries.parquet",
    "solvation": "extensions/solvation/counts.parquet",
    "sdc4_mechanisms": "extensions/sdc4_mechanisms/timeseries.parquet",
}


@dataclass(frozen=True)
class FusedOutput:
    name: str
    path: Path
    run_id: str | None = None

    @property
    def checkpoint_key(self) -> str:
        """Unambiguous current-format identifier, independent of old markers."""
        return self.name if self.run_id is None else f"{self.name}-{self.run_id}"


def fused_output(repout: Path, name: str, *, run_id: str | None = None) -> FusedOutput:
    """Return the expected artifact for a named analysis instance."""
    if name == "sasa" and run_id is not None:
        from .config import _validate_path_component
        run_id = _validate_path_component(str(run_id), "sasa run ID")
        return FusedOutput(name, Path(repout) / "profiles" / "sasa" /
                           run_id / "residue_sasa.parquet", run_id)
    if name == "residue_membrane_distance":
        if run_id is None:
            raise ValueError("residue_membrane_distance requires a run_id")
        from .config import _validate_path_component
        run_id = _validate_path_component(str(run_id), "distance run ID")
        return FusedOutput(
            name, Path(repout) / "profiles" / name / run_id / "distance_matrix_A.npy",
            run_id,
        )
    if run_id is not None:
        raise ValueError(f"{name} does not yet support independently named runs")
    return FusedOutput(name, Path(repout) / FUSED_OUTPUTS[name])


def fused_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Plan enabled fused outputs without creating readers or opening trajectories."""
    # Imports stay local: registry is also used by postprocessing.
    from .pipeline.fused_membrane import enabled_names as membrane_enabled
    from .pipeline.fused_extensions import enabled_names as extensions_enabled
    from .config import analysis_runs
    from .pipeline.fused_distance import enabled_names as distance_enabled
    from .mechanisms.sdc4_mechanisms import enabled as sdc4_enabled
    entries = [fused_output(repout, name) for name in membrane_enabled(cfg)]
    entries.extend(fused_output(repout, name) for name in extensions_enabled(cfg)
                   if name != "sasa")
    if "sasa" in extensions_enabled(cfg):
        entries.extend(fused_output(repout, "sasa", run_id=run.run_id)
                       if run.run_id != "default" else fused_output(repout, "sasa")
                       for run in analysis_runs(cfg, "sasa"))
    entries.extend(fused_output(repout, "residue_membrane_distance", run_id=run_id)
                   for run_id in distance_enabled(cfg))
    if sdc4_enabled(cfg):
        entries.append(fused_output(repout, "sdc4_mechanisms"))
    return tuple(entries)



def analysis_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """The same persisted contracts apply to fused and standalone scheduling."""
    entries = list(fused_output_plan(cfg, repout))
    names = {entry.name for entry in entries}
    if cfg['analyses'].get('dssp', {}).get('enabled') and 'dssp' not in names:
        entries.append(fused_output(repout, 'dssp'))
    if cfg['analyses'].get('sasa', {}).get('enabled') and 'sasa' not in names:
        entries.extend(fused_output(repout, 'sasa', run_id=run.run_id)
                       if run.run_id != 'default' else fused_output(repout, 'sasa')
                       for run in analysis_runs(cfg, 'sasa'))
    return tuple(entries)


def required_artifacts(entry: FusedOutput, cfg: dict) -> tuple[Path, ...]:
    """Required companions, without opening trajectories or inventing a manifest.

    Sparse event products that are legitimately absent after no events are not
    required. Large numerical tables are checked by the caller using metadata.
    """
    directory, name = entry.path.parent, entry.name
    companions = {
        'lipid_contacts': ('matrix.parquet', 'selection_audit.json'),
        'lipid_enrichment': ('counts.parquet', 'selection_audit.parquet'),
        'membrane_order': ('site_mapping.parquet', 'profile_index.parquet'),
        'hydrogen_bonds': ('count_timeseries.parquet',),
        'lipid_hydrogen_bonds': ('count_timeseries.parquet', 'molecule_audit.parquet'),
        'residue_membrane_distance': ('time_ps.npy', 'residues.parquet',
                                    'median_profile.parquet', 'selection_audit.parquet',
                                    'membrane_atom_indices.npy'),
        'solvation': ('rdf.parquet', 'target_residues.parquet', 'metadata.json',
                      'species_timeseries.parquet', 'target_summary.parquet',
                      'shell_candidates.parquet'),
        'sdc4_mechanisms': ('selection_audit.json', 'lipid_density_maps.npz',
                            'hs_conditioned_solvation.parquet'),
    }.get(name, ())
    paths = [entry.path, *(directory / filename for filename in companions)]
    spec = cfg['analyses'].get(name, {})
    if name == 'lipid_enrichment' and spec.get('area_reference'):
        paths.append(directory/'area_reference.json')
    if name == 'sasa':
        spec = next(run.spec for run in analysis_runs(cfg, name)
                    if run.run_id == (entry.run_id or 'default'))
        if spec.get('atom_sasa_selection'):
            paths.append(directory/'atom_sasa.parquet')
        if spec.get('reference_run') is not None:
            paths.extend(directory/f for f in ('residue_burial.parquet', 'residue_burial_summary.parquet'))
    if name == 'sdc4_mechanisms' and spec.get('motif_lipid_residence', False):
        paths.extend(directory/f for f in ('motif_lipid_contact_events.parquet',
            'motif_lipid_residence.parquet', 'motif_lipid_residence_intermittent.parquet',
            'lipid_region_specificity.parquet', 'lipid_residue_contacts.parquet'))
    if name == 'solvation':
        options = spec.get('options') or {}
        bridges = options.get('bridges', spec.get('preset', 'standard') != 'standard')
        if options.get('events', spec.get('preset', 'standard') != 'standard') or bridges or options.get('water_hbond_bridges'):
            paths.extend(directory/f for f in ('episodes.parquet', 'survival.parquet'))
        if options.get('bridge_target_pairs'):
            paths.append(directory/'bridge_pair_summary.parquet')
        if spec.get('aqueous_reference'):
            paths.append(directory/'aqueous_reference.parquet')
        if spec.get('grid_frames'):
            paths.append(directory/'water_density_3d.parquet')
    return tuple(paths)
