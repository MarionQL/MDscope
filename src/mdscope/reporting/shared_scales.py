"""Shared visual scales for panels that measure the same quantity.

Call only on homogeneous y-axis panels. Matplotlib's auto limits include
uncertainty bands, reference lines, and normal padding; preserve those bounds.
"""
from __future__ import annotations

import numpy as np


def synchronize_y_axes(axes, *, symmetric: bool = False):
    """Give all panels a common y-range without treating missing panels as data.

    Returns shared limits or None when there is no measured y-extent.
    Excludes unused and text-only axes; a no-data panel still receives the
    shared range so a visual comparison cannot imply different scales.
    """
    panels = list(np.asarray(axes, dtype=object).flat)
    populated = []
    for ax in panels:
        if not ax.has_data():
            continue
        bounds = np.asarray(ax.get_ylim(), dtype=float)
        if np.isfinite(bounds).all() and bounds[0] < bounds[1]:
            populated.append(bounds)
    if not populated:
        return None
    low = min(float(pair[0]) for pair in populated)
    high = max(float(pair[1]) for pair in populated)
    if symmetric:
        bound = max(abs(low), abs(high))
        low, high = -bound, bound
    for ax in panels:
        ax.set_ylim(low, high)
    return low, high
