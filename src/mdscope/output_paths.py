"""Persisted analysis output paths shared by producers and the registry.

This dependency-free module prevents output producers from importing the
higher-level registry (which may inspect production enabled predicates).
Paths and on-disk formats are existing MDscope contracts, not new defaults.
"""

from .contracts.outputs import FUSED_OUTPUTS as FUSED_OUTPUTS
from .contracts.outputs import RUN_SCOPED_OUTPUTS as RUN_SCOPED_OUTPUTS

# Scalar readers, state-feature provenance and authoring use the same sources.
from .contracts.outputs import SAVED_SCALAR_SOURCE_FILES as SAVED_SCALAR_SOURCE_FILES
