"""Small, explicit leaflet partition shared by legacy median-based analyses."""
from __future__ import annotations
import numpy as np


def median_leaflet_mask(coordinates):
    """Upper when coordinate >= median; no PBC reconstruction or relabeling.

    Intended for initial, already centered bilayers only. For PBC-crossing
    or drifting membranes use the circular classifier in lipid_geometry.
    """
    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        raise ValueError("Leaflet coordinates must be finite 1D data with >=2 entries")
    return values >= np.median(values)
