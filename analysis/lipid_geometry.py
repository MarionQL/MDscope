"""Compatibility exports for lipid geometry and PBC primitives.

Canonical implementations live in :mod:`mdscope.analysis.core.lipid_geometry`.
"""
from .core.lipid_geometry import (
    AnchorTracker,
    classify_leaflets,
    aligned_leaflets,
    periodic_center,
    anchor_groups,
    centers,
    planar_offsets,
    lateral_distances,
)
