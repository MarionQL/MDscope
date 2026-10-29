"""Compatibility exports for rigid-body alignment utilities.

Canonical implementations live in :mod:`mdscope.analysis.core.alignment`.
The old module remains importable for existing workflows.
"""
from .core.alignment import (
    _valid_box,
    _fit_transform,
    _apply_transform,
    resolved_selection_counts,
)
