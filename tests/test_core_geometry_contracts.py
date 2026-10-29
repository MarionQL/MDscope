"""Geometry deduplication contracts: retain differences that are scientific."""
import numpy as np
import pytest
from MDAnalysis.lib.distances import minimize_vectors

from mdscope.analysis.core.geometry import minimum_image_displacements
from mdscope.analysis.extended.common import minimum_image, lateral_distance as extended_lateral
from mdscope.analysis.mechanistic_common import planar_distances
from mdscope.analysis.biophysical_common import lateral_distance as pairwise_lateral

BOX = np.array([40., 50., 60., 90., 90., 90.])


def test_minimum_image_matches_independent_mdanalysis_reference():
    offsets = np.array([[39., -49., 2.], [-21., 26., 31.]], dtype=np.float64)
    np.testing.assert_allclose(
        minimum_image_displacements(offsets, BOX),
        minimize_vectors(offsets, BOX),
    )


def test_extended_lateral_uses_fixed_xy_and_float32():
    xyz = np.array([[39., 3., 9.], [2., 48., 20.]])
    anchor = np.array([1., 1., 5.])
    expected = minimize_vectors((xyz-anchor).astype(np.float32), BOX)
    np.testing.assert_allclose(extended_lateral(xyz, anchor, BOX),
                               np.linalg.norm(expected[:, :2], axis=1))
    assert minimum_image(xyz-anchor, BOX).dtype == np.float32


def test_mechanistic_lateral_normal_is_configurable():
    xyz = np.array([[39., 3., 9.], [2., 48., 20.]])
    anchor = np.array([1., 1., 5.])
    offsets = minimize_vectors(xyz-anchor, BOX)
    for normal in (0, 1, 2):
        expected = np.linalg.norm(np.delete(offsets, normal, axis=1), axis=1)
        np.testing.assert_allclose(planar_distances(xyz, anchor, BOX, normal), expected)


def test_pairwise_lateral_is_not_equivalent_to_pointwise_lateral():
    a = np.array([[1., 2., 3.], [5., 7., 11.]])
    b = np.array([[4., 5., 6.], [8., 9., 10.], [2., 3., 4.]])
    assert pairwise_lateral(a, b, BOX).shape == (2, 3)
    assert extended_lateral(a, b[0], BOX).shape == (2,)


def test_invalid_geometry_fails_explicitly():
    with pytest.raises(ValueError, match="final coordinate"):
        minimum_image_displacements(np.ones((4, 2)), BOX)
    with pytest.raises(ValueError, match="periodic box"):
        minimum_image_displacements(np.ones((4, 3)), [0, 50, 60, 90, 90, 90])
