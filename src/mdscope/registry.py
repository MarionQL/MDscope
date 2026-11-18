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
from .config import analysis_enabled, analysis_runs
from .observables.residue_membrane_distance import run_residue_membrane_distance
from .observables.membrane_order import run_membrane_order
from .observables.structure_profiles import run_structure_profiles
from .interactions.hydrogen_bonds import run_hydrogen_bonds
from .mechanisms.local_environment import run_local_environment
from .mechanisms.potassium_bridge_audit import run_potassium_bridge_audit
from .observables.lipid_enrichment import FOOTPRINT_OUTPUT, run_lipid_enrichment
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

    if analysis_enabled(cfg, "potassium_bridge_audit"):
        yield AnalysisTask("potassium_bridge_audit", run_potassium_bridge_audit)

    for name,fn in (("lipid_enrichment",run_lipid_enrichment),
                    ("bilayer_normal",run_bilayer_normal),
                    ("lipid_hydrogen_bonds",run_lipid_hydrogen_bonds)):
        if cfg.get("analyses",{}).get(name,{}).get("enabled",False):
            output = FUSED_OUTPUTS[name] if name != "lipid_hydrogen_bonds" else None
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
    "lipid_enrichment": f"{FOOTPRINT_OUTPUT}/timeseries.parquet",
    "hydrogen_bonds": "interactions/hydrogen_bonds/events.parquet",
    "lipid_hydrogen_bonds": "interactions/lipid_hydrogen_bonds/events.parquet",
    "dssp": "profiles/dssp/assignments.parquet",
    "sasa": "profiles/sasa/residue_sasa.parquet",
    "local_environment": "extensions/local_environment/timeseries.parquet",
    "solvation": "extensions/solvation/counts.parquet",
    "hydration": "extensions/hydration/counts.parquet",
    "sdc4_mechanisms": "extensions/sdc4_mechanisms/timeseries.parquet",
    "potassium_bridge_audit": "extensions/potassium_bridge_audit/topology_audit.json",
}


# Precise schema details live with the artifact registry, not the GUI.
# A missing column/shape description means the corresponding producer has
# a dynamic or as-yet-uncataloged contract; it is *not* a zero-length table.
def artifact_schema(entry: "FusedOutput", path: Path) -> dict:
    """Known persisted fields, shapes and units for a registered artifact."""
    name = Path(path).name
    analysis = entry.name
    fields = None
    shape = None
    units = None
    dtype = None
    if analysis == "sasa":
        if name == "residue_sasa.parquet":
            fields = ["time_ps", "topology_residue_index", "residue", "sasa_nm2"]
            units = {"time_ps": "ps", "sasa_nm2": "nm^2"}
        elif name == "atom_sasa.parquet":
            fields = ["time_ps", "atom_index", "topology_residue_index",
                      "atom_name", "sasa_nm2"]
            units = {"time_ps": "ps", "sasa_nm2": "nm^2"}
        elif name == "residue_burial_summary.parquet":
            fields = ["topology_residue_index", "residue", "mean_buried_sasa_nm2",
                      "std_buried_sasa_nm2", "mean_buried_fraction",
                      "mean_environment_sasa_nm2", "mean_isolated_sasa_nm2"]
            units = {"mean_buried_sasa_nm2": "nm^2",
                     "std_buried_sasa_nm2": "nm^2",
                     "mean_buried_fraction": "fraction"}
    elif analysis == "residue_membrane_distance":
        if name == "distance_matrix_A.npy":
            shape, dtype, units = ["sampled_frames", "selected_residues"], "float32", "angstrom"
        elif name == "time_ps.npy":
            shape, units = ["sampled_frames"], "ps"
        elif name == "membrane_atom_indices.npy":
            shape, dtype = ["membrane_atoms"], "int64"
        elif name == "residues.parquet":
            fields = ["analysis_instance", "analysis_label", "residue_order",
                      "topology_residue_index", "segid", "resid", "resname",
                      "residue", "n_selected_atoms"]
        elif name == "median_profile.parquet":
            fields = ["analysis_instance", "analysis_label", "residue_order",
                      "topology_residue_index", "segid", "resid", "resname",
                      "residue", "n_selected_atoms", "n_frames",
                      "mean_min_distance_A", "sd_min_distance_A",
                      "median_min_distance_A", "q1_min_distance_A",
                      "q3_min_distance_A", "iqr_min_distance_A",
                      "min_distance_A", "max_distance_A"]
            units = {"distance_metrics_suffix_A": "angstrom"}
    elif analysis == "membrane_order":
        if name == "profile_index.parquet":
            fields = ["profile", "file"]
        elif name == "site_mapping.parquet":
            fields = ["resname", "tail", "carbon_position", "atom_name", "n_lipids"]
    elif analysis == "solvation":
        if name == "counts.parquet":
            from .observables.solvation import COUNT_COLUMNS
            fields = list(COUNT_COLUMNS)
            units = {"time_ps": "ps", "n_solvent": "entities"}
        elif name == "rdf.parquet":
            fields = ["solvent", "target_class", "site_class",
                      "r_inner_A", "r_outer_A", "pair_count", "g_bulk"]
            units = {"r_inner_A": "angstrom", "r_outer_A": "angstrom",
                     "g_bulk": "dimensionless; full-cell density reference"}
        elif name == "target_residues.parquet":
            fields = ["target_id", "topology_residue_index"]
    elif analysis == "hydration":
        if name in {"counts.parquet", "aggregate_counts.parquet", "target_metadata.parquet"}:
            from .mechanisms.hs_hydration import HydrationObserver
            fields = list({
                "counts.parquet": HydrationObserver.COUNT_COLUMNS,
                "aggregate_counts.parquet": HydrationObserver.AGGREGATE_COLUMNS,
                "target_metadata.parquet": HydrationObserver.TARGET_COLUMNS,
            }[name])
            units = {"time_ps": "ps", "water_cutoff_A": "angstrom",
                     "water_count": "unique water entities"}
        elif name == "hs_sugar_hydration.parquet":
            from .mechanisms.hs_hydration import HydrationObserver
            fields = list(HydrationObserver.SUGAR_COLUMNS)
            units = {"time_ps": "ps", "water_cutoff_A": "angstrom",
                     "water_count": "unique water entities"}
    return {
        "columns": fields,
        "columns_status": "verified_fixed" if fields is not None else "not_registered",
        "shape_axes": shape,
        "dtype": dtype,
        "units": units,
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
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        entries.append(fused_output(repout, "potassium_bridge_audit"))
    return tuple(entries)


def required_artifacts(entry: FusedOutput, cfg: dict) -> tuple[Path, ...]:
    """Required companions, without opening trajectories or inventing a manifest.

    Sparse event products that are legitimately absent after no events are not
    required. Large numerical tables are checked by the caller using metadata.
    """
    directory, name = entry.path.parent, entry.name
    companions = {
        'lipid_contacts': ('matrix.parquet', 'selection_audit.json'),
        'lipid_enrichment': ('counts.parquet', 'selection_audit.parquet', 'reference.json'),
        'membrane_order': ('site_mapping.parquet', 'profile_index.parquet'),
        'hydrogen_bonds': ('count_timeseries.parquet',),
        'lipid_hydrogen_bonds': ('count_timeseries.parquet', 'molecule_audit.parquet'),
        'residue_membrane_distance': ('time_ps.npy', 'residues.parquet',
                                    'median_profile.parquet', 'selection_audit.parquet',
                                    'membrane_atom_indices.npy'),
        'hydration': ('aggregate_counts.parquet', 'target_metadata.parquet'),
        'solvation': ('rdf.parquet', 'target_residues.parquet', 'metadata.json',
                      'species_timeseries.parquet', 'target_summary.parquet',
                      'shell_candidates.parquet'),
        'potassium_bridge_audit': ('chemical_sites.parquet', 'bridge_frames.parquet',
                                  'bridge_summary.parquet', 'rdf.parquet',
                                  'rdf_peak_candidates.parquet'),
        'sdc4_mechanisms': ('selection_audit.json', 'lipid_density_maps.npz',
                            'hs_conditioned_solvation.parquet'),
    }.get(name, ())
    paths = [entry.path, *(directory / filename for filename in companions)]
    spec = cfg['analyses'].get(name, {})
    if name == 'sasa':
        spec = next(run.spec for run in analysis_runs(cfg, name)
                    if run.run_id == (entry.run_id or 'default'))
        if spec.get('atom_sasa_selection'):
            paths.append(directory/'atom_sasa.parquet')
        if spec.get('reference_run') is not None:
            paths.extend(directory/f for f in ('residue_burial.parquet', 'residue_burial_summary.parquet'))
    if name == 'sdc4_mechanisms' and spec.get('residue_occlusion', False):
        from .mechanisms.residue_occlusion import OCCLUSION_FILE
        paths.append(directory / OCCLUSION_FILE)
    if name == 'sdc4_mechanisms' and spec.get('motif_lipid_residence', False):
        paths.extend(directory/f for f in ('motif_lipid_contact_events.parquet',
            'motif_lipid_residence.parquet', 'motif_lipid_residence_intermittent.parquet',
            'lipid_region_specificity.parquet', 'lipid_residue_contacts.parquet'))
    if name == 'hydration' and cfg['analyses'].get('sdc4_mechanisms', {}).get('hs_sugar_water_cutoff_A') is not None:
        from .mechanisms.hs_hydration import OUTPUT
        paths.append(Path(entry.path).parents[2] / OUTPUT)
    if name == 'solvation':
        options = spec.get('options') or {}
        bridges = options.get('bridges', spec.get('preset', 'standard') != 'standard')
        if options.get('events', spec.get('preset', 'standard') != 'standard') or bridges or options.get('water_hbond_bridges'):
            paths.extend(directory/f for f in ('episodes.parquet', 'survival.parquet'))
        if options.get('bridge_target_pairs'):
            paths.append(directory/'bridge_pair_summary.parquet')
        if options.get("bridge_spatial_ion") is not None:
            paths.extend(directory/name for name in ("bridge_spatial_reference.parquet",
                                                      "bridge_pair_timeseries.parquet"))
        if spec.get('aqueous_reference'):
            paths.append(directory/'aqueous_reference.parquet')
        if spec.get('grid_frames'):
            paths.append(directory/'water_density_3d.parquet')
    return tuple(paths)
