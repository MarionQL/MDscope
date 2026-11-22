"""Authoritative persisted-output contracts for MDscope analyses.

Planning and execution both consume these paths, companion artifacts and
schemas. Executable replica tasks live in pipeline.replica_tasks.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import analysis_enabled, analysis_runs
from .output_paths import FUSED_OUTPUTS, RUN_SCOPED_OUTPUTS


# Precise schema details live with the artifact registry, not the GUI.
# A missing column/shape description means the corresponding producer has
# a dynamic or as-yet-uncataloged contract; it is *not* a zero-length table.
def artifact_schema(entry: FusedOutput, path: Path) -> dict:
    """Known persisted fields, shapes and units for a registered artifact."""
    name = Path(path).name
    analysis = entry.name
    fields = None
    shape = None
    units = None
    dtype = None
    row_grain = None
    array_keys = None
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
    elif analysis in {"hydrogen_bonds", "lipid_hydrogen_bonds"} and name == "events.parquet":
        fields = ["time_ps", "donor", "hydrogen", "acceptor"]
        fields += (["donor_chain", "acceptor_chain", "scope", "interaction_class"]
                   if analysis == "hydrogen_bonds" else
                   ["donor_molecule_id", "acceptor_molecule_id", "donor_species",
                    "acceptor_species", "species_pair", "scope", "interaction_class"])
        fields += ["distance_A", "angle_deg"]
        units = {"time_ps": "ps", "distance_A": "angstrom", "angle_deg": "degrees"}
        row_grain = "one detected donor-hydrogen-acceptor event per analyzed frame per replicate"
    elif analysis == "dssp" and name == "assignments.parquet":
        fields = ["time_ps", "topology_residue_index", "residue", "dssp"]
        units = {"time_ps": "ps", "dssp": "secondary-structure code"}
        row_grain = "one residue per analyzed frame per replicate"
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
        elif name == "glycan_attachments.parquet":
            from .mechanisms.hs_hydration import HydrationObserver
            fields = list(HydrationObserver.ATTACHMENT_COLUMNS)
            row_grain = "one covalently attached glycan chain per replicate"
        elif name == "hs_sugar_hydration.parquet":
            from .mechanisms.hs_hydration import HydrationObserver
            fields = list(HydrationObserver.SUGAR_COLUMNS)
            units = {"time_ps": "ps", "water_cutoff_A": "angstrom",
                     "water_count": "unique water entities"}
    # Wide observer tables have topology/configuration-dependent metric columns.
    # Their sampled-time column and row population are still exact contracts.
    wide_series = {"bilayer_normal", "membrane_order", "lipid_enrichment", "lipid_contacts", "local_environment", "sdc4_mechanisms"}
    dynamic_columns = False
    if name == "timeseries.parquet" and analysis in wide_series:
        fields = ["time_ps"]
        units = {"time_ps": "ps"}
        row_grain = "one analyzed frame per replicate; metric columns depend on configured targets/topology"
        dynamic_columns = True
        if analysis == "bilayer_normal":
            fields += ["bilayer_normal_valid_fraction", "bilayer_normal_leaflet_changes",
                       "bilayer_normal_upper_mean_deg", "bilayer_normal_upper_sd_deg",
                       "bilayer_normal_lower_mean_deg", "bilayer_normal_lower_sd_deg"]
        elif analysis == "lipid_contacts":
            fields += ["lipid_contact_pair_count", "lipid_contact_cross_leaflet_pair_count", "lipid_contact_leaflet_changes"]
            dynamic_columns = False
    if analysis == "lipid_enrichment" and name == "selection_audit.parquet":
        fields = ["species", "resindex", "leaflet"]
        row_grain = "one selected lipid anchor with initial leaflet per replicate"
    if analysis == "lipid_contacts" and name == "matrix.parquet":
        fields = ["time_ps", "leaflet", "species_a", "species_b", "molecules_a", "molecules_b",
                  "pair_count_undirected", "possible_pairs", "total_contact_pairs", "neighbors_per_molecule", "enrichment"]
        row_grain = "one species pair per leaflet per sampled frame per replicate"
        units = {"time_ps": "ps", "enrichment": "finite-population reference ratio", "neighbors_per_molecule": "neighbors/molecule"}
    if analysis == "lipid_enrichment" and name == "counts.parquet":
        fields = ["time_ps", "leaflet", "species", "lower_A", "upper_A", "count_local", "n_local", "count_bulk", "n_bulk", "count_contact", "mean_min_distance_A", "bulk_fraction", "enrichment_fold", "area_density_ratio", "shell_accessible_area_A2", "leaflet_accessible_area_A2"]
        row_grain = "one species per leaflet per spatial shell per frame per replicate"
    if analysis == "solvation":
        if name == "species_timeseries.parquet":
            fields = ["solvent", "target_class", "site_class", "time_ps", "target_count", "mean_n_solvent", "total_target_counts", "hydrated_fraction"]
            row_grain = "one solvent/target/site class per analyzed frame per replicate"
        elif name == "target_summary.parquet":
            fields = ["solvent", "target_id", "target_class", "site_class", "n_frames", "mean_n_solvent", "sd_n_solvent", "hydrated_fraction"]
            row_grain = "one solvent/target identity per replicate"
        elif name == "episodes.parquet":
            from .observables.solvation_common import EPISODE_COLUMNS
            fields = list(EPISODE_COLUMNS)
            row_grain = "one sampled solvent/target contact episode; censor flags retained"
        elif name == "survival.parquet":
            fields = ["solvent", "target_id", "duration_ps", "at_risk", "dissociations", "right_censored", "survival_probability"]
            row_grain = "one target/solvent per observed episode duration; not independent replicas"
        elif name == "shell_candidates.parquet":
            fields = ["solvent", "target_class", "site_class", "suggested_cutoff_A", "peak_r_A", "minimum_g_bulk", "status"]
            row_grain = "one chemical target/site class per solvent; diagnostic candidate, not binding inference"
        elif name == "rdf.parquet":
            row_grain = "one chemical target/site class per solvent per radial shell per replicate"
    if analysis == "sdc4_mechanisms" and name == "lipid_density_maps.npz":
        array_keys = {key: ["species_leaflet", "x_bin", "y_bin"] for key in ("counts", "density_per_nm2", "enrichment_fold")}
        array_keys.update(bulk_density_per_nm2=["species_leaflet"], species=["species_leaflet"],
                          leaflet=["species_leaflet"], x_edges_A=["x_edge"], y_edges_A=["y_edge"], nframes=[])
        row_grain = "species/leaflet resolved spatial histogram; array dimensions depend on configured map grid"
        units = {"density_per_nm2": "lipids/nm^2", "enrichment_fold": "ratio", "x_edges_A": "angstrom", "y_edges_A": "angstrom"}
    if analysis == "sasa" and name == "atom_sasa.parquet":
        row_grain = "one selected atom per analyzed frame per replicate"
    if analysis in {"hydrogen_bonds", "lipid_hydrogen_bonds"} and name == "count_timeseries.parquet":
        prefix = "hydrogen_bond" if analysis == "hydrogen_bonds" else "lipid_hydrogen_bond"
        fields = ["time_ps", prefix + "_total_count"]
        units = {"time_ps": "ps", prefix + "_total_count": "unique events/frame"}
        row_grain = "one analyzed frame per replicate; additional interaction-class count columns depend on detected chemistry"
        dynamic_columns = True
    if analysis == "lipid_hydrogen_bonds" and name == "molecule_audit.parquet":
        fields = ["raw_hbond_events", "same_residue_events", "same_molecule_events", "cross_residue_same_molecule_events", "intermolecular_events", "selected_atoms", "selected_molecules"]
        row_grain = "one molecule-filtering audit per replicate"
    if analysis == "residue_membrane_distance" and name == "selection_audit.parquet":
        fields = ["analysis_instance", "analysis_label", "protein_selection", "membrane_selection", "protein_atom_count", "protein_residue_count", "membrane_atom_count", "membrane_residue_count", "protein_atom_index_hash", "membrane_atom_index_hash", "membrane_resnames"]
        row_grain = "one selection audit per configured run per replicate"
    if analysis == "sdc4_mechanisms" and name == "hs_conditioned_solvation.parquet":
        fields = ["hs_chain", "membrane_contact_state", "frames", "hydration_water_per_heavy_atom_mean"]
        row_grain = "one HS chain per observed membrane-contact state per replicate; optional ion-count columns depend on configured species"
        dynamic_columns = True
    if analysis == "potassium_bridge_audit":
        fields = {
            "chemical_sites.parquet": ["owner", "kind", "resindex", "resid", "resname", "center_name", "atom_index", "atom_name", "charge_e"],
            "bridge_frames.parquet": ["frame", "time_ps", "acid_patch", "hs_chain", "hs_kind", "shift_x", "shift_y", "contact_cutoff_A", "steric_cutoff_A", "min_target_separation_A", "eligible_ions", "bridging_ions"],
            "bridge_summary.parquet": ["acid_patch", "hs_chain", "hs_kind", "contact_cutoff_A", "steric_cutoff_A", "reference", "occupancy", "eligible_ions", "bridging_ions", "bridge_fraction_of_eligible_ion_samples"],
            "rdf.parquet": ["kind", "r_inner_A", "r_outer_A", "pair_count", "g_box", "shared_water_pair_count"],
            "rdf_peak_candidates.parquet": ["kind", "frames", "smoothing_A", "first_peak_A", "intervening_minimum_A", "second_peak_A"],
        }.get(name)
        row_grain = {"chemical_sites.parquet": "one anionic oxygen site", "bridge_frames.parquet": "one target pair/cutoff/reference shift per frame", "bridge_summary.parquet": "one target pair/cutoff/reference aggregated over frames", "rdf.parquet": "one chemical site kind per radial bin", "rdf_peak_candidates.parquet": "one chemical site kind; diagnostic peaks"}.get(name)
    if name.endswith(".json"):
        row_grain = "one metadata/audit object per analysis run per replicate; keys may depend on topology"
    if row_grain is None and fields is not None:
        if "time_ps" in fields and "topology_residue_index" in fields:
            row_grain = "one selected target residue per analyzed frame per replicate"
        elif name.endswith("summary.parquet"):
            row_grain = "one configured target identity per replicate; summary over sampled frames"
        else:
            row_grain = "identity/coordinate fields shown; exact population depends on configured target and artifact role"
    dynamic_columns = dynamic_columns or (analysis == "solvation" and name == "rdf.parquet") or (
        analysis == "residue_membrane_distance" and name == "median_profile.parquet"
    )
    return {
        "row_grain": row_grain,
        "contract_status": ("tabular_base_plus_dynamic" if dynamic_columns else "tabular_fixed" if fields is not None else "array_axes" if shape is not None or array_keys is not None else "metadata_or_dynamic_unresolved"),
        "unresolved_reason": (None if fields is not None or shape is not None or array_keys is not None else "Metadata keys or topology-dependent fields are not fixed; no numerical values or dimensions inferred"),
        "columns": fields,
        "columns_status": ("verified_base_plus_optional" if dynamic_columns
                           else "verified_fixed" if fields is not None
                           else "not_registered"),
        "array_keys": array_keys,
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
    """Return the expected artifact for a named analysis instance.

    A run named ``default`` maps to the unscoped product when one exists.
    """
    if run_id == "default" and name in FUSED_OUTPUTS:
        run_id = None
    if run_id is None:
        if name not in FUSED_OUTPUTS:
            raise ValueError(f"{name} requires a run_id")
        return FusedOutput(name, Path(repout) / FUSED_OUTPUTS[name])
    if name not in RUN_SCOPED_OUTPUTS:
        raise ValueError(f"{name} does not yet support independently named runs")
    from .config import _validate_path_component
    run_id = _validate_path_component(str(run_id), f"{name} run ID")
    directory, filename = RUN_SCOPED_OUTPUTS[name]
    return FusedOutput(name, Path(repout) / directory / run_id / filename, run_id)


def fused_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Plan enabled fused outputs without creating readers or opening trajectories."""
    # Share production scheduling predicates without importing heavy producers.
    from .pipeline.analysis_selection import (
        distance_enabled_names as distance_enabled,
    )
    from .pipeline.analysis_selection import (
        extensions_enabled_names as extensions_enabled,
    )
    from .pipeline.analysis_selection import (
        membrane_enabled_names as membrane_enabled,
    )
    from .pipeline.analysis_selection import (
        sdc4_enabled,
    )
    entries = [fused_output(repout, name) for name in membrane_enabled(cfg)]
    extensions = extensions_enabled(cfg)
    entries.extend(fused_output(repout, name) for name in extensions
                   if name != "sasa")
    if "sasa" in extensions:
        entries.extend(fused_output(repout, "sasa", run_id=run.run_id)
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
                       for run in analysis_runs(cfg, 'sasa'))
    if analysis_enabled(cfg, "potassium_bridge_audit"):
        entries.append(fused_output(repout, "potassium_bridge_audit"))
    return tuple(entries)


def saved_output_plan(cfg: dict, repout: Path) -> tuple[FusedOutput, ...]:
    """Readiness contracts, including products written by the core stream.

    Paths come from the same contracts used by authoring and statistics; this
    does not add work to the executable fused schedule.
    """
    entries = list(analysis_output_plan(cfg, repout))
    names = {entry.name for entry in entries}
    for contract in scalar_contracts(cfg) + structural_profile_contracts(cfg):
        name = contract["analysis"]
        if name not in names:
            entries.append(FusedOutput(name, Path(repout) / contract["path_within_replica"]))
            names.add(name)
    return tuple(entries)


def core_saved_contracts(cfg, name):
    """Core stream companions and columns, from the public output contracts."""
    if name in FUSED_OUTPUTS:
        return []
    return [item for item in scalar_contracts(cfg) + structural_profile_contracts(cfg)
            if item["analysis"] == name]


def required_artifacts(entry: FusedOutput, cfg: dict) -> tuple[Path, ...]:
    """Required companions, without opening trajectories or inventing a manifest.

    Sparse event products that are legitimately absent after no events are not
    required. Large numerical tables are checked by the caller using metadata.
    """
    directory, name = entry.path.parent, entry.name
    core = core_saved_contracts(cfg, name)
    if core:
        # Core primaries are one directory below the replica root.
        return tuple(dict.fromkeys(entry.path.parent.parent / item["path_within_replica"]
                                   for item in core))
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
    if name == 'hydration' and spec.get('enabled') and spec.get('glycans'):
        # Optional for HS-only runs, so pre-existing SDC4 outputs stay reusable.
        from .mechanisms.hs_hydration import ATTACHMENTS
        paths.append(directory / ATTACHMENTS)
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


STRUCTURAL_SCALAR_METRICS = {"rmsd": "rmsd_A",
                             "radius_of_gyration": "radius_of_gyration_A"}


def structural_scalar_contracts(cfg: dict) -> list[dict]:
    """Configured scalar identities from the shared streaming summary writer.

    Group RMSD measures the intersection with selections.features. Rg uses
    its own selection intersection. The global scalar remains written even
    when a groups list is supplied; lists filter grouped rows only.
    """
    from .config import analysis_groups
    from .reporting.plot_ids import describe_plot, safe_plot_stem

    groups = analysis_groups(cfg)
    result = []
    for analysis, metric in STRUCTURAL_SCALAR_METRICS.items():
        spec = cfg["analyses"].get(analysis, {})
        if not spec.get("enabled", False):
            continue
        selected = spec.get("groups")
        ids = [None] + [key for key in groups if selected is None or key in selected]
        for group in ids:
            stem = safe_plot_stem(metric + (f"__group_{group}" if group else ""))
            measurement = (cfg["selections"]["features"] if analysis == "rmsd"
                           else spec.get("selection") or cfg["selections"]["protein"])
            fit = (groups[group].get("alignment_selection") if group else None)
            fitted = analysis == "rmsd" and cfg["analyses"]["alignment"].get("enabled", True)
            fit = ((fit or cfg["analyses"]["alignment"].get("selection")
                    or cfg["selections"]["alignment"]) if fitted else None)
            result.append({
                "analysis": analysis, "metric": metric, "group": group,
                "label": groups[group]["label"] if group else "Global",
                "columns": ["time_ps", *(["group"] if group else []), metric],
                "columns_status": "verified_subset_of_shared_table",
                "units": {"time_ps": "ps", metric: "angstrom"},
                "row_grain": "one analyzed frame per replicate" + (" per group" if group else ""),
                "path_within_replica": "summary/" + (
                    "group_timeseries.parquet" if group else "timeseries.parquet"),
                "measurement_selection": measurement,
                "group_selection": groups[group]["selection"] if group else None,
                "selection_semantics": "intersection with group_selection" if group else "direct",
                "alignment_selection": fit,
                "coordinate_source": "fitted" if fitted else "source",
                "statistics_selector": {"metric": metric, "group": group},
                "figures": [dict(describe_plot(f"time_series:{stem}:{variant}"),
                                 requires_statistics=variant not in {"timeseries", "distribution"})
                            for variant in ("timeseries", "distribution", "autocorrelation",
                                            "bootstrap_table", "comparison_table")],
            })
    return result


def structural_profile_contracts(cfg):
    """Streaming RMSF products; separate from frame-wise scalar tables."""
    from .config import analysis_groups
    if not cfg["analyses"]["rmsf"].get("enabled", True):
        return []
    base = ["atom_index", "topology_residue_index", "residue", "atom_name", "rmsf_A"]
    result = []
    for grouped in ([False, True] if analysis_groups(cfg) else [False]):
        result.append({
            "analysis": "rmsf", "run_id": "default", "role": "conditional_summary",
            "format": "Parquet", "columns_status": "verified_fixed",
            "columns": ["group", *base] if grouped else [*base, "equivalent_bfactor_A2"],
            "path_within_replica": "raw/group_rmsf.parquet" if grouped else "raw/rmsf.parquet",
            "row_grain": "one selected feature atom per replicate" + (" per group" if grouped else ""),
            "units": {"rmsf_A": "angstrom", **({} if grouped else {"equivalent_bfactor_A2": "angstrom^2"})},
            "coordinate_source": "fitted" if cfg["analyses"]["alignment"].get("enabled", True) else "source",
            "statistics_requirement": "Temporal RMSF inference additionally requires saved aligned-coordinate cache; a static RMSF table alone is insufficient.",
        })
    return result


def scalar_contracts(cfg):
    """Deterministic pre-run scalar identities, including membrane and DSSP."""
    from .reporting.plot_ids import describe_plot
    result = structural_scalar_contracts(cfg)
    specs = []
    membrane = cfg["analyses"]["membrane"]
    if membrane.get("enabled"):
        if membrane.get("area_per_lipid", {}).get("enabled", True):
            specs.extend(("membrane", f"area_per_lipid{suffix}_A2", "angstrom^2")
                         for suffix in ("_upper", "_lower", ""))
        if membrane.get("thickness", {}).get("enabled", True):
            specs.append(("membrane", "membrane_thickness_A", "angstrom"))
    if cfg["analyses"]["dssp"].get("enabled"):
        specs.extend(("dssp", f"{kind}_fraction", "fraction")
                     for kind in ("helix", "sheet", "coil"))
    if cfg["analyses"].get("bilayer_normal", {}).get("enabled"):
        specs += [("bilayer_normal", "bilayer_normal_valid_fraction", "fraction"),
                  ("bilayer_normal", "bilayer_normal_leaflet_changes", "lipids")]
        for leaflet in ("upper", "lower"):
            specs.extend(("bilayer_normal", f"bilayer_normal_{leaflet}_{stat}_deg", "degrees")
                         for stat in ("mean", "sd"))
    if cfg["analyses"].get("lipid_contacts", {}).get("enabled"):
        specs.extend(("lipid_contacts", metric, "pairs" if "pair" in metric else "lipids") for metric in (
            "lipid_contact_pair_count", "lipid_contact_cross_leaflet_pair_count", "lipid_contact_leaflet_changes"))
    for analysis, metric, unit in specs:
        result.append({
            "analysis": analysis, "metric": metric, "group": None, "label": metric,
            "columns": ["time_ps", metric], "columns_status": "verified_subset_of_shared_table",
            "units": {"time_ps": "ps", metric: unit},
            "row_grain": "one analyzed frame per replicate", "path_within_replica": (FUSED_OUTPUTS[analysis] if analysis in {"lipid_contacts", "bilayer_normal"} else "summary/timeseries.parquet"),
            "measurement_selection": (membrane.get("headgroup_selection") if analysis == "membrane"
                                      else cfg["selections"]["protein"] if analysis == "dssp"
                                      else cfg["analyses"][analysis].get("selection") or cfg["analyses"][analysis].get("species_selections")),
            "group_selection": None, "selection_semantics": "producer-defined population; see scientific contracts",
            "alignment_selection": None, "coordinate_source": "source",
            "statistics_selector": {"metric": metric, "group": None},
            "figures": [dict(describe_plot(f"time_series:{metric}:{variant}"),
                             requires_statistics=variant not in {"timeseries", "distribution"})
                        for variant in ("timeseries", "distribution", "autocorrelation", "bootstrap_table", "comparison_table")],
        })
    return result


def configured_profile_variables(cfg, *, include_optional=False):
    """Pre-run identities consumed by temporal_rmsf/structure_profiles readers."""
    from .config import analysis_groups
    variables = []
    def add(analysis, profile, unit, prerequisite, run_id="default"):
        variables.append(dict(analysis=analysis, profile=profile, run_id=run_id,
                              units=unit, prerequisite=prerequisite,
                              statistics_selector={"profile": profile}))
    if cfg["analyses"]["rmsf"].get("enabled", True):
        add("rmsf", "rmsf", "angstrom", "Saved aligned-coordinate cache and compatible atom identities required")
        for group in analysis_groups(cfg):
            add("rmsf", f"rmsf__group_{group}", "angstrom", "Group feature intersection and saved coordinate cache required")
    for run in analysis_runs(cfg, "sasa"):
        suffix = "" if run.run_id == "default" else f"__{run.run_id}"
        add("sasa", "sasa_mean" + suffix, "nm^2", "Per-frame residue SASA with compatible residue identities", run.run_id)
        if run.run_id != "default" and run.spec.get("reference_run"):
            add("sasa", "sasa_buried_nm2" + suffix, "nm^2", "Compatible isolated/reference-run burial tables required", run.run_id)
    if cfg["analyses"]["dssp"].get("enabled", False):
        for kind in ("helix", "sheet", "coil"):
            add("dssp", f"dssp_{kind}_occupancy", "fraction", "Per-frame residue assignments required")
    for run in analysis_runs(cfg, "density_1d"):
        for key, raw in run.spec.get("selections", {}).items():
            target = {"selection": raw} if isinstance(raw, str) else raw
            unit = "amu/A^3" if target.get("mode", "number") == "mass" else "1/A^3"
            add("density_1d", f"density_1d__{run.run_id}__{key}", unit,
                "Compatible physical bin grids and per-frame density cache required", run.run_id)
            if target.get("groups"):
                for group in analysis_groups(cfg):
                    add("density_1d", f"density_1d__{run.run_id}__{key}__{group}", unit,
                        "Conditional: target/group intersection must be nonempty; compatible grids required", run.run_id)
    for run in analysis_runs(cfg, "residue_membrane_distance"):
        names = [("mean_min_distance_A", "angstrom"), ("median_min_distance_A", "angstrom")]
        if include_optional or run.spec.get("contact_cutoff_A") is not None:
            names.append(("fraction_frames_within_cutoff", "fraction"))
        prefixes = [f"residue_membrane_distance__{run.run_id}"]
        if run.spec.get("groups"):
            prefixes.extend(f"residue_membrane_distance__{run.run_id}__group_{group}" for group in analysis_groups(cfg))
        for prefix in prefixes:
            for metric, unit in names:
                add("residue_membrane_distance", f"{prefix}__{metric}", unit,
                    "Compatible residue identities and per-frame matrix; grouped outputs require nonempty intersections", run.run_id)
    return variables


def density_output_contracts(cfg):
    """Mean density tables saved by DensityProfileAccumulator.finalize."""
    from .config import analysis_groups, profile_dir
    outputs = []
    for run in analysis_runs(cfg, "density_1d"):
        for key, raw in run.spec.get("selections", {}).items():
            target = {"selection": raw} if isinstance(raw, str) else raw
            unit = "amu/A^3" if target.get("mode", "number") == "mass" else "1/A^3"
            outputs.append(dict(analysis="density_1d", run_id=run.run_id, role="conditional_summary",
                format="Parquet", columns_status="verified_fixed", row_grain="one spatial bin per target per replicate",
                columns=[f"{run.spec.get('axis', 'z')}_A", "density", "analysis_instance", "selection_key",
                         "selection_label", "mode", "density_unit", "n_frames"],
                units={"density": unit}, path_within_replica=str(profile_dir(Path(), "density_1d", run.run_id) / f"{key}.parquet")))
            if target.get("groups"):
                for group in analysis_groups(cfg):
                    output = dict(outputs[-1])
                    output["path_within_replica"] = str(profile_dir(Path(), "density_1d", run.run_id) / f"{key}__{group}.parquet")
                    output["row_grain"] = "conditional on nonempty target/group intersection; one spatial bin per replicate"
                    outputs.append(output)
    return outputs


def state_feature_catalog(cfg):
    """Available static identities plus explicit saved-only discovery contracts.

    Topology-dependent molecule/bridge identities are templates until the common
    scalar reader discovers saved rows; they are never advertised as measured.
    """
    from .output_paths import SAVED_SCALAR_SOURCE_FILES
    candidates = []
    for scalar in scalar_contracts(cfg):
        if scalar["path_within_replica"].startswith("summary/"):
            candidates.append({"source": "summary", "group": scalar["group"],
                "metric": scalar["metric"], "units": scalar["units"][scalar["metric"]],
                "analysis": scalar["analysis"], "availability": "conditional_on_saved_data"})
    solvation = cfg["analyses"].get("solvation", {})
    if solvation.get("enabled"):
        solvents = (["water"] if solvation.get("water_selection") else [])
        solvents += ["ions/" + name for name in (solvation.get("ion_selections") or {})]
        if solvation.get("ion_selection"):
            solvents.append("ions")
        for kind, targets in (solvation.get("targets") or {}).items():
            for label in targets:
                for solvent in solvents:
                    for metric, units in (("mean_n_solvent", "molecules"), ("hydrated_fraction", "fraction")):
                        candidates.append({"source": "solvation",
                            "group": f"solvation/{solvent}/{kind}/{label}", "metric": metric,
                            "units": units, "analysis": "solvation",
                            "availability": "conditional_on_saved_data"})
    hydration = cfg["analyses"].get("hydration", {})
    for name in (hydration.get("glycans") or []) if hydration.get("enabled") else []:
        for metric, units in (("unique_water_count", "molecules"),
                              ("unique_water_per_site", "molecules/site")):
            candidates.append({"source": "hydration", "group": f"hydration/{name}/total/all",
                               "metric": metric, "units": units, "analysis": "hydration",
                               "availability": "conditional_on_saved_data"})
    sources = {
        "summary": {"group": "null or shared group", "metrics": "numerical saved columns",
                    "prerequisite": "Saved frame-wise global/group summary with selected numerical columns"},
        "solvation": {"group": "solvation/<solvent>/<target_class>/<site_class>",
                      "metrics": ["mean_n_solvent", "hydrated_fraction"],
                      "prerequisite": "Saved per-frame solvent/site aggregates; requested motif targets must match topology"},
        "hydration": {"group": "hydration/<compartment>/<level>/<label>",
                      "metrics": ["unique_water_count", "unique_water_per_site"],
                      "prerequisite": "Saved aggregate hydration counts with topology-resolved group identities"},
        "solvation_bridges": {"group": "solvation_bridges/<JSON [solvent,target_a,target_b]>",
                      "metrics": ["bridge_entity_count", "bridge_present_fraction", "observed_minus_spatial_bridge_fraction"],
                      "prerequisite": "Saved bridge-pair traces; spatial contrast requires its optional reference"},
        "lipid_contacts": {"group": "lipid_contacts/<leaflet>/<species_a>/<species_b>",
                      "metrics": ["neighbors_per_molecule", "enrichment"],
                      "prerequisite": "Saved species/leaflet matrix; unavailable enrichment remains nonfinite"},
    }
    for name, source in sources.items():
        source["paths_within_replica"] = list(SAVED_SCALAR_SOURCE_FILES[name])
        source["identity_status"] = "template_until_saved_product_validation"
    return {"candidates": candidates, "sources": sources,
            "time_requirement": "Exact common time grid within each replicate; complete feature vectors only"}
