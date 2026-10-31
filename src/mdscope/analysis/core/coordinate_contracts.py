"""Coordinate provenance contracts for the streaming analysis pipeline.

Labels describe how coordinates are used, not whether a trajectory file is
actually unwrapped or centered. Those properties cannot reliably be inferred
from a frame of positions.
"""
from __future__ import annotations

from enum import Enum


class CoordinateSpace(str, Enum):
    SOURCE = "source"                  # untouched coordinates from the supplied trajectory
    MSD_SOURCE = "msd_source"          # independent reader when msd_trajectory is supplied
    FITTED = "fitted"                  # rigid-body fit, stored in separate arrays


RAW_FRAME_CONSUMERS = frozenset({
    "fused_membrane", "fused_extensions", "fused_distances",
    "sdc4_mechanisms", "density", "membrane",
})


def require_coordinate_space(consumer: str, actual: CoordinateSpace, expected: CoordinateSpace):
    """Fail explicitly if a caller attempts an incompatible coordinate handoff."""
    if actual != expected:
        raise ValueError(
            f"{consumer} requires {expected.value} coordinates, "
            f"received {actual.value}"
        )


def msd_source_space(has_separate_reader: bool) -> CoordinateSpace:
    """The membrane MSD source must never be a protein-fitted coordinate array."""
    return CoordinateSpace.MSD_SOURCE if has_separate_reader else CoordinateSpace.SOURCE
