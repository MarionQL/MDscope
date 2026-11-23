"""Dependency-free persisted artifact descriptors (issue #260).

Producer implementations, execution planning and Studio consume these stable
paths and schemas. Config-dependent optional outputs are resolved by registry.
No mdscope package imports belong in this layer.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# One entry per output-producing analysis. Preserve historic path and ordering
# contracts; conditional extras are resolved by registry.required_artifacts.
@dataclass(frozen=True)
class AnalysisOutputDescriptor:
    """Static persisted-output contract for one analysis."""
    name: str
    primary_path: str | None = None
    named_run: tuple[str, str] | None = None
    companions: tuple[str, ...] = ()

    def schema(self, path: Path) -> dict:
        return artifact_schema(self, path)


OUTPUT_DESCRIPTORS = {
    "lipid_contacts": AnalysisOutputDescriptor(
        "lipid_contacts", "extensions/lipid_contacts/timeseries.parquet",
        companions=("matrix.parquet", "selection_audit.json")),
    "membrane_order": AnalysisOutputDescriptor(
        "membrane_order", "profiles/membrane_order/timeseries.parquet",
        companions=("site_mapping.parquet", "profile_index.parquet")),
    "bilayer_normal": AnalysisOutputDescriptor(
        "bilayer_normal", "extensions/bilayer_normal/timeseries.parquet"),
    "lipid_enrichment": AnalysisOutputDescriptor(
        "lipid_enrichment", "extensions/lipid_footprint/timeseries.parquet",
        companions=("counts.parquet", "selection_audit.parquet", "reference.json")),
    "hydrogen_bonds": AnalysisOutputDescriptor(
        "hydrogen_bonds", "interactions/hydrogen_bonds/events.parquet",
        companions=("count_timeseries.parquet",)),
    "lipid_hydrogen_bonds": AnalysisOutputDescriptor(
        "lipid_hydrogen_bonds", "interactions/lipid_hydrogen_bonds/events.parquet",
        companions=("count_timeseries.parquet", "molecule_audit.parquet")),
    "dssp": AnalysisOutputDescriptor(
        "dssp", "profiles/dssp/assignments.parquet"),
    "sasa": AnalysisOutputDescriptor(
        "sasa", "profiles/sasa/residue_sasa.parquet",
        named_run=("profiles/sasa", "residue_sasa.parquet")),
    "local_environment": AnalysisOutputDescriptor(
        "local_environment", "extensions/local_environment/timeseries.parquet"),
    "solvation": AnalysisOutputDescriptor(
        "solvation", "extensions/solvation/counts.parquet",
        companions=("rdf.parquet", "target_residues.parquet", "metadata.json",
                    "species_timeseries.parquet", "target_summary.parquet",
                    "shell_candidates.parquet")),
    "hydration": AnalysisOutputDescriptor(
        "hydration", "extensions/hydration/counts.parquet",
        companions=("aggregate_counts.parquet", "target_metadata.parquet")),
    "sdc4_mechanisms": AnalysisOutputDescriptor(
        "sdc4_mechanisms", "extensions/sdc4_mechanisms/timeseries.parquet",
        companions=("selection_audit.json", "lipid_density_maps.npz",
                    "hs_conditioned_solvation.parquet")),
    "potassium_bridge_audit": AnalysisOutputDescriptor(
        "potassium_bridge_audit", "extensions/potassium_bridge_audit/topology_audit.json",
        companions=("chemical_sites.parquet", "bridge_frames.parquet",
                    "bridge_summary.parquet", "rdf.parquet", "rdf_peak_candidates.parquet")),
    "residue_membrane_distance": AnalysisOutputDescriptor(
        "residue_membrane_distance",
        named_run=("profiles/residue_membrane_distance", "distance_matrix_A.npy"),
        companions=("time_ps.npy", "residues.parquet", "median_profile.parquet",
                    "selection_audit.parquet", "membrane_atom_indices.npy")),
}

# Legacy public tables are derived views: they are not separate registrations.
FUSED_OUTPUTS = {
    name: descriptor.primary_path
    for name, descriptor in OUTPUT_DESCRIPTORS.items()
    if descriptor.primary_path is not None
}
RUN_SCOPED_OUTPUTS = {
    name: descriptor.named_run
    for name, descriptor in OUTPUT_DESCRIPTORS.items()
    if descriptor.named_run is not None
}


SOLVATION_COUNT_COLUMNS = ["frame", "time_ps", "solvent", "target_id", "target_class",
                 "site_class", "n_solvent"]
SOLVATION_EPISODE_COLUMNS = ["solvent", "target_id", "solvent_id", "start_frame", "end_frame",
                   "duration_ps", "left_censored", "right_censored",
                   "sampled_duration_ps", "duration_convention"]
HYDRATION_COUNT_COLUMNS = [
        "frame", "time_ps", "compartment", "target_id", "water_count",
        "water_cutoff_A", "n_sites", "water_per_site",
    ]
HYDRATION_AGGREGATE_COLUMNS = [
        "frame", "time_ps", "compartment", "level", "group_id",
        "n_targets", "water_count", "water_cutoff_A", "n_sites", "water_per_site",
    ]
HYDRATION_TARGET_COLUMNS = [
        "compartment", "target_id", "chain", "sugar_position", "resname",
        "resid", "topology_residue_index", "n_sites", "water_cutoff_A",
    ]
HYDRATION_ATTACHMENT_COLUMNS = [
        "compartment", "chain", "protein_chain", "protein_resid", "protein_resname",
        "protein_atom", "protein_atom_index", "protein_topology_residue_index",
        "n_sugars", "sugar_types",
    ]
HYDRATION_SUGAR_COLUMNS = [
        "hs_chain", "sugar_position", "sugar_type", "sugar_resid",
        "oxygen_count", "time_ps", "water_cutoff_A", "water_count",
    ]
HS_SUGAR_HYDRATION_PATH = "extensions/sdc4_mechanisms/hs_sugar_hydration.parquet"
GLYCAN_ATTACHMENTS_FILE = "glycan_attachments.parquet"
RESIDUE_OCCLUSION_FILE = "residue_occlusion.parquet"

# Canonical persisted scalar source paths used by saved-data readers and schemas.
SAVED_SCALAR_SOURCE_FILES = {
    "summary": ("summary/timeseries.parquet", "summary/group_timeseries.parquet"),
    "solvation": ("extensions/solvation/species_timeseries.parquet",),
    "hydration": ("extensions/hydration/aggregate_counts.parquet",),
    "solvation_bridges": ("extensions/solvation/bridge_pair_timeseries.parquet",),
    "lipid_contacts": ("extensions/lipid_contacts/matrix.parquet",),
}

MEMBRANE_ANALYSES = ("membrane_order", "bilayer_normal", "lipid_enrichment", "lipid_contacts")
EXTENSION_ANALYSES = ("hydrogen_bonds", "lipid_hydrogen_bonds", "dssp", "sasa", "local_environment", "solvation")


def artifact_schema(entry, path: Path) -> dict:
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
            fields = list(SOLVATION_COUNT_COLUMNS)
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
            fields = list({
                "counts.parquet": HYDRATION_COUNT_COLUMNS,
                "aggregate_counts.parquet": HYDRATION_AGGREGATE_COLUMNS,
                "target_metadata.parquet": HYDRATION_TARGET_COLUMNS,
            }[name])
            units = {"time_ps": "ps", "water_cutoff_A": "angstrom",
                     "water_count": "unique water entities"}
        elif name == "glycan_attachments.parquet":
            fields = list(HYDRATION_ATTACHMENT_COLUMNS)
            row_grain = "one covalently attached glycan chain per replicate"
        elif name == "hs_sugar_hydration.parquet":
            fields = list(HYDRATION_SUGAR_COLUMNS)
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
            fields = list(SOLVATION_EPISODE_COLUMNS)
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


# Legacy products have intentionally smaller reuse requirements than the full
# writer schema. Keep this compatibility policy beside that schema.
MINIMUM_REUSE_COLUMNS = {
    "solvation": ("frame", "time_ps", "target_id", "n_solvent"),
    "dssp": ("time_ps", "topology_residue_index", "residue", "dssp"),
    "sasa": ("time_ps", "topology_residue_index", "residue", "sasa_nm2"),
    "hydrogen_bonds": ("time_ps", "donor", "acceptor"),
    "lipid_hydrogen_bonds": ("time_ps", "donor", "acceptor"),
}


def reuse_columns(entry, path):
    """Required columns for an existing product, preserving legacy reuse."""
    if Path(path) == entry.path:
        return MINIMUM_REUSE_COLUMNS.get(entry.name, ("time_ps",))
    return artifact_schema(entry, path)["columns"]
