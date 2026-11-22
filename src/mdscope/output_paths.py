"""Persisted analysis output paths shared by producers and the registry.

This dependency-free module prevents output producers from importing the
higher-level registry (which may inspect production enabled predicates).
Paths and on-disk formats are existing MDscope contracts, not new defaults.
"""

FUSED_OUTPUTS = {
    "lipid_contacts": "extensions/lipid_contacts/timeseries.parquet",
    "membrane_order": "profiles/membrane_order/timeseries.parquet",
    "bilayer_normal": "extensions/bilayer_normal/timeseries.parquet",
    "lipid_enrichment": "extensions/lipid_footprint/timeseries.parquet",
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

# Analyses whose named runs write separate products at <directory>/<run_id>/<file>.
# SASA's run named "default" keeps the unscoped FUSED_OUTPUTS path; residue
# distances have no unscoped product, so every instance needs a run ID.
RUN_SCOPED_OUTPUTS = {
    "sasa": ("profiles/sasa", "residue_sasa.parquet"),
    "residue_membrane_distance": ("profiles/residue_membrane_distance", "distance_matrix_A.npy"),
}

# Scalar readers, state-feature provenance and authoring use the same sources.
SAVED_SCALAR_SOURCE_FILES = {
    "summary": ("summary/timeseries.parquet", "summary/group_timeseries.parquet"),
    "solvation": ("extensions/solvation/species_timeseries.parquet",),
    "hydration": ("extensions/hydration/aggregate_counts.parquet",),
    "solvation_bridges": ("extensions/solvation/bridge_pair_timeseries.parquet",),
    "lipid_contacts": ("extensions/lipid_contacts/matrix.parquet",),
}
