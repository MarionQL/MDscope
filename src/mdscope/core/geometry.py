"""Shared minimum-image geometry primitives.

The minimum-image operation is shared; callers retain their own axis,
shape, dtype and box-validation semantics. This module does not guess a
membrane normal or silently turn an invalid box into non-periodic geometry.
"""
from __future__ import annotations
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors


def minimum_image_displacements(delta, box, *, dtype=None):
    """Apply triclinic-compatible minimum image to displacement vectors (Å).

    Parameters
    ----------
    delta : array-like, (..., 3)
        Displacement vectors, not absolute atom coordinates.
    box : array-like, (6,)
        MDAnalysis periodic dimensions: lengths in Å, angles in degrees.
    dtype : numpy dtype, optional
        Preserve caller's intentional precision conversion. None leaves
        NumPy's input dtype unchanged where supported.
    """
    vectors = np.asarray(delta, dtype=dtype)
    dimensions = np.asarray(box, dtype=float)
    if vectors.ndim < 1 or vectors.shape[-1] != 3:
        raise ValueError("Displacements require a final coordinate dimension of 3")
    if dimensions.shape != (6,) or not np.all(np.isfinite(dimensions)) or np.any(dimensions[:3] <= 0):
        raise ValueError("Minimum-image distances require a finite valid (6,) periodic box")
    return minimize_vectors(vectors, dimensions)
