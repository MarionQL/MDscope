"""Cross-module numerical contracts which must *not* be accidentally unified."""
import numpy as np
import pytest

from mdscope.analysis.biophysical_common import lateral_distance
from mdscope.analysis.mechanistic_common import planar_distances
from mdscope.analysis.core.lipid_geometry import lateral_distances


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
