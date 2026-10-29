"""Independent geometry and topology contracts for issue #31."""
import numpy as np
import pytest
from MDAnalysis.lib.distances import minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors

from mdscope.analysis.core.geometry import minimum_image_displacements
from mdscope.analysis.core.lipid_geometry import planar_offsets, lateral_distances, periodic_center


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
