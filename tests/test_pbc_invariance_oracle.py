"""The wrap-invariance oracle tolerates hard-cutoff ties but still catches PBC errors (#250)."""
import numpy as np
from pbc_invariance import (
    WRAP_TOLERANCE_DEG,
    cutoff_tie_anchors,
    wrap_normal_difference,
)

from mdscope.observables.bilayer_normal import fit_local_normals

BOX = np.array([61.3, 59.7, 100.0, 90.0, 90.0, 90.0])
RADIUS = 15.0


def _membrane():
    """Gently tilted upper leaflet on an irregular lattice; no accidental cutoff ties."""
    rng = np.random.default_rng(250)
    xs, ys = np.meshgrid(np.arange(0.2, 61.3, 4.37), np.arange(0.3, 59.7, 4.61))
    xyz = np.column_stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)])
    xyz[:, :2] += rng.uniform(-0.4, 0.4, size=(len(xyz), 2))
    xyz[:, 2] = 40.0 + 2.0 * np.sin(2 * np.pi * xyz[:, 0] / BOX[0]) + rng.normal(0, 0.3, len(xyz))
    return xyz, np.ones(len(xyz), bool)


def _box_shift(xyz):
    shifted = xyz.copy()
    shifted[::4, 0] += BOX[0]
    shifted[1::4, 1] -= BOX[1]
    shifted[2::4, 2] += BOX[2]
    return shifted


def _place_tie(xyz):
    """Move anchor 1 to exactly RADIUS from anchor 0 (3D), above the local plane."""
    tied = xyz.copy()
    direction = np.array([0.6, 0.0, 0.8])
    tied[1] = tied[0] + RADIUS * direction / np.linalg.norm(direction)
    return tied


def test_lattice_has_no_ties_and_exact_invariance():
    xyz, upper = _membrane()
    assert not cutoff_tie_anchors(xyz, upper, BOX, RADIUS).any()
    result = wrap_normal_difference(xyz, upper, _box_shift(xyz), upper, BOX,
                                    radius=RADIUS, min_neighbors=4)
    assert result["tie_anchors"] == 0
    assert result["max_non_tie_diff_deg"] < WRAP_TOLERANCE_DEG
    assert result["non_tie_validity_changes"] == 0


def test_cutoff_tie_is_detected_and_isolated_from_the_invariance_check():
    xyz, upper = _membrane()
    tied = _place_tie(xyz)
    ties = cutoff_tie_anchors(tied, upper, BOX, RADIUS)
    assert ties[0] and ties[1]
    shifted = _box_shift(tied)
    # Emulate the float32-rounding flip: the tied neighbour lands just outside.
    shifted[1] += 2e-5 * (tied[1] - tied[0]) / RADIUS
    flipped = wrap_normal_difference(tied, upper, shifted, upper, BOX,
                                     radius=RADIUS, min_neighbors=4)
    assert flipped["max_tie_diff_deg"] > WRAP_TOLERANCE_DEG, "construction must change a tie anchor"
    assert flipped["max_non_tie_diff_deg"] < WRAP_TOLERANCE_DEG
    assert flipped["non_tie_validity_changes"] == 0
    assert 2 <= flipped["tie_anchors"] <= 4


def test_genuine_pbc_error_still_fails_the_oracle():
    xyz, upper = _membrane()
    # Suppress minimum-image wrapping without allocating a million-Angstrom
    # neighbor-search grid in the deliberately broken implementation.
    huge = np.array([1000., 1000., 1000., 90.0, 90.0, 90.0])
    assert np.all(np.ptp(np.vstack([xyz, _box_shift(xyz)]), axis=0) < huge[:3] / 2)

    def no_minimum_image(points, leaflet, box, **kwargs):
        return fit_local_normals(points, leaflet, huge, **kwargs)

    result = wrap_normal_difference(xyz, upper, _box_shift(xyz), upper, BOX, radius=RADIUS,
                                    min_neighbors=4, fit=no_minimum_image)
    assert result["tie_anchors"] == 0
    assert (result["max_non_tie_diff_deg"] > WRAP_TOLERANCE_DEG
            or result["non_tie_validity_changes"] > 0), result
