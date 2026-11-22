"""Persisted analysis output paths shared by producers and the registry.

This dependency-free module prevents output producers from importing the
higher-level registry (which may inspect production enabled predicates).
Paths and on-disk formats are existing MDscope contracts, not new defaults.
"""

from .contracts.outputs import FUSED_OUTPUTS, RUN_SCOPED_OUTPUTS

# Scalar readers, state-feature provenance and authoring use the same sources.
SAVED_SCALAR_SOURCE_FILES = {
    "summary": ("summary/timeseries.parquet", "summary/group_timeseries.parquet"),
    "solvation": ("extensions/solvation/species_timeseries.parquet",),
    "hydration": ("extensions/hydration/aggregate_counts.parquet",),
    "solvation_bridges": ("extensions/solvation/bridge_pair_timeseries.parquet",),
    "lipid_contacts": ("extensions/lipid_contacts/matrix.parquet",),
}
