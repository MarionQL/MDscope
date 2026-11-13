"""Geometry deduplication contracts: retain differences that are scientific."""
import numpy as np
import pytest
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors

from mdscope.core import alignment as core_alignment
from mdscope.core import lipid_geometry as core_lipid
from mdscope.core.geometry import minimum_image_displacements
from mdscope.core.lipid_geometry import planar_offsets, lateral_distances, periodic_center
from mdscope.mechanisms.mechanistic_common import planar_distances
from mdscope.biophysical_common import lateral_distance

BOX = np.array([40., 50., 60., 90., 90., 90.])


@pytest.mark.parametrize("box", [BOX, np.array([30., 35., 40., 75., 85., 105.])])
def test_minimum_image_matches_independent_mdanalysis_reference(box):
    offsets = np.array([[39., -49., 2.], [-21., 26., 31.]], dtype=np.float64)
    np.testing.assert_allclose(
        minimum_image_displacements(offsets, box),
        minimize_vectors(offsets, box),
    )


def test_minimum_image_preserves_explicit_float32_precision():
    xyz = np.array([[39., 3., 9.], [2., 48., 20.]])
    anchor = np.array([1., 1., 5.])
    expected = minimize_vectors((xyz-anchor).astype(np.float32), BOX)
    actual = minimum_image_displacements(xyz-anchor, BOX, dtype=np.float32)
    np.testing.assert_allclose(np.linalg.norm(actual[:, :2], axis=1),
                               np.linalg.norm(expected[:, :2], axis=1))
    assert actual.dtype == np.float32


def test_mechanistic_lateral_normal_is_configurable():
    xyz = np.array([[39., 3., 9.], [2., 48., 20.]])
    anchor = np.array([1., 1., 5.])
    offsets = minimize_vectors(xyz-anchor, BOX)
    for normal in (0, 1, 2):
        expected = np.linalg.norm(np.delete(offsets, normal, axis=1), axis=1)
        np.testing.assert_allclose(planar_distances(xyz, anchor, BOX, normal), expected)


def test_invalid_geometry_fails_explicitly():
    with pytest.raises(ValueError, match="final coordinate"):
        minimum_image_displacements(np.ones((4, 2)), BOX)
    with pytest.raises(ValueError, match="periodic box"):
        minimum_image_displacements(np.ones((4, 3)), [0, 50, 60, 90, 90, 90])



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


def test_pairwise_and_pointwise_planar_distances_are_distinct_contracts():
    box = np.array([10., 20., 30., 90., 90., 90.])
    points = np.array([[9., 1., 2.], [1., 2., 3.]])
    anchors = np.array([[1., 1., 2.], [9., 2., 3.]])
    all_pairs = lateral_distance(points, anchors, box)
    assert all_pairs.shape == (2, 2)
    np.testing.assert_allclose(all_pairs.diagonal(), [2., 2.])
    one_reference = lateral_distances(points, anchors[0], box)
    np.testing.assert_allclose(one_reference, all_pairs[:, 0])
    np.testing.assert_allclose(
        planar_distances(points, anchors[0], box),
        one_reference,
    )


def test_triclinic_planar_methods_preserve_distinct_policies():
    box = np.array([15., 17., 19., 75., 85., 105.])
    points = np.array([[14., 4., 9.], [1., 2., 3.]])
    anchor = np.array([1., 2., 3.])
    with pytest.raises(ValueError, match="orthorhombic"):
        lateral_distance(points, anchor[None, :], box)
    np.testing.assert_allclose(
        planar_distances(points, anchor, box),
        lateral_distances(points, anchor, box),
    )


@pytest.mark.parametrize("box", [
    [30., 35., 40., 90., 90., 90.],
    [30., 35., 40., 75., 85., 105.],
])
def test_planar_offsets_match_independent_mdanalysis_reference(box):
    box = np.asarray(box, float)
    points = np.array([[29., 2., 5.], [1., 33., 38.], [12., 20., 9.]])
    reference = np.array([2., 3., 4.])
    expected = minimize_vectors(points - reference, box)
    for normal in range(3):
        np.testing.assert_allclose(
            planar_offsets(points, reference, box, normal),
            np.delete(expected, normal, axis=-1), rtol=0, atol=1e-5,
        )
        np.testing.assert_allclose(
            lateral_distances(points, reference, box, normal),
            np.linalg.norm(np.delete(expected, normal, axis=-1), axis=1),
            rtol=0, atol=1e-5,
        )


def test_triclinic_periodic_center_at_crossing():
    box = np.array([30., 35., 40., 75., 85., 105.])
    lattice = np.asarray(triclinic_vectors(box), float)
    fractional = np.array([[.98, .99, .97], [.02, .01, .03]])
    center = periodic_center(fractional @ lattice, box)
    frac_center = (center @ np.linalg.inv(lattice)) % 1.
    np.testing.assert_allclose(np.minimum(frac_center, 1 - frac_center), [0., 0., 0.], atol=1e-5)


def test_invalid_boxes_and_empty_periodic_center():
    for dimensions in ([0, 30, 30, 90, 90, 90], [30, 30, 30], [np.nan, 30, 30, 90, 90, 90]):
        with pytest.raises(ValueError):
            minimum_image_displacements([[1., 2., 3.]], dimensions)
    with pytest.raises(ValueError, match="nonempty"):
        periodic_center(np.empty((0, 3)), [30, 30, 30, 90, 90, 90])
