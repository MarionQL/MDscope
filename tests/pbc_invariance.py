"""PBC wrap-invariance oracle for local bilayer normals, robust to hard-cutoff ties.

``fit_local_normals`` selects neighbours with a hard 3D cutoff. Shifting an
anchor by whole box vectors is a physical no-op, but it perturbs coordinates by
float32 rounding (~1e-5 A for ~300 A boxes). A neighbour lying within that
rounding of the cutoff can then enter or leave a neighbourhood and legitimately
change that anchor's fitted tilt (observed: 1.03 deg, carb test data frame 25).
Such *tie* anchors cannot be compared exactly; every other anchor must be
invariant. Ties are counted and reported, never silently dropped.
"""
from __future__ import annotations

import numpy as np
from MDAnalysis.lib.distances import capped_distance

from mdscope.observables.bilayer_normal import fit_local_normals

# Measured on the test-data-v1 trajectories, a whole-box shift changes pair
# distances by <= 2.9e-5 A (float32 rounding at ~185-300 A coordinates). The band
# is ~10x that, yet far below lipid spacing, so genuine PBC errors cannot hide in
# it. Over all 401 real frames it flags at most 12 of 1632 anchors per frame
# (122/201 carb and 117/200 nocarb frames have any tie); a 1e-3 A band flagged ~2x more.
CUTOFF_TIE_TOLERANCE_A = 2.5e-4
# Non-tie anchors must agree to this after a whole-box shift (float rounding
# alone gives ~1e-7-1e-5 deg; one neighbour change gives ~1 deg).
WRAP_TOLERANCE_DEG = 1e-3


def cutoff_tie_anchors(xyz, upper, box, radius, tol=CUTOFF_TIE_TOLERANCE_A):
    """Boolean mask of anchors with a same-leaflet neighbour within ``tol`` of ``radius``."""
    xyz = np.asarray(xyz, float)
    upper = np.asarray(upper, bool)
    ties = np.zeros(len(xyz), bool)
    pairs = capped_distance(xyz, xyz, max_cutoff=radius + tol, min_cutoff=radius - tol,
                            box=box, return_distances=False)
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        keep = (i != j) & (upper[i] == upper[j])
        ties[i[keep]] = True
        ties[j[keep]] = True
    return ties


def wrap_normal_difference(xyz, upper, shifted_xyz, shifted_upper, box, *, radius,
                           min_neighbors, fit=fit_local_normals, tol=CUTOFF_TIE_TOLERANCE_A):
    """Compare local normals before/after a whole-box shift, separating cutoff ties.

    Returns the maximum absolute tilt difference over non-tie anchors, the number
    of non-tie anchors whose validity (finite vs NaN) changed, and tie statistics.
    """
    before, _ = fit(xyz, upper, box, radius=radius, min_neighbors=min_neighbors, normal=2)
    after, _ = fit(shifted_xyz, shifted_upper, box, radius=radius,
                   min_neighbors=min_neighbors, normal=2)
    ties = cutoff_tie_anchors(xyz, upper, box, radius, tol)
    both = np.isfinite(before) & np.isfinite(after)
    diff = np.zeros(len(before))
    diff[both] = np.abs(before[both] - after[both])
    validity_changed = np.isfinite(before) != np.isfinite(after)
    return {
        "max_non_tie_diff_deg": float(diff[~ties].max(initial=0.0)),
        "non_tie_validity_changes": int((validity_changed & ~ties).sum()),
        "tie_anchors": int(ties.sum()),
        "max_tie_diff_deg": float(diff[ties].max(initial=0.0)),
        "n_anchors": len(before),
    }
