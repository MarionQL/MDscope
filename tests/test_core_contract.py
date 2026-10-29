"""Compatibility and scientific regressions for the core-module extraction."""
import numpy as np

from mdscope.analysis import alignment, lipid_geometry
from mdscope.analysis.core import alignment as core_alignment
from mdscope.analysis.core import lipid_geometry as core_lipid


def test_legacy_alignment_imports_are_identical():
    assert alignment._fit_transform is core_alignment._fit_transform
    assert alignment._apply_transform is core_alignment._apply_transform
    assert alignment._valid_box is core_alignment._valid_box
    assert alignment.resolved_selection_counts is core_alignment.resolved_selection_counts


def test_legacy_lipid_geometry_imports_are_identical():
    for name in ("AnchorTracker", "classify_leaflets", "aligned_leaflets",
                 "periodic_center", "anchor_groups", "centers",
                 "planar_offsets", "lateral_distances"):
        assert getattr(lipid_geometry, name) is getattr(core_lipid, name)


def test_alignment_rigid_rotation_and_translation():
    rng = np.random.default_rng(123)
    reference = rng.normal(size=(20, 3))
    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                         [np.sin(angle), np.cos(angle), 0.],
                         [0., 0., 1.]])
    mobile = reference @ rotation + np.array([15., -3., 4.])
    transform = core_alignment._fit_transform(mobile, reference)
    np.testing.assert_allclose(
        core_alignment._apply_transform(mobile, transform), reference,
        atol=1e-12,
    )


def test_periodic_center_stable_under_independent_wrapping():
    box = np.array([100., 100., 100., 90., 90., 90.])
    xyz = np.array([[99., 20., 30.], [1., 21., 30.], [98., 19., 31.]])
    shifted = xyz.copy()
    shifted[0, 0] -= 100.
    shifted[1, 0] += 100.
    np.testing.assert_allclose(
        core_lipid.periodic_center(xyz, box),
        core_lipid.periodic_center(shifted, box),
        atol=1e-10,
    )
