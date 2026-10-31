"""Legacy median leaflet partition and PBC-specific behavior."""
import numpy as np
import pytest

from mdscope.core.leaflet_partition import median_leaflet_mask
from mdscope.core.lipid_geometry import classify_leaflets


@pytest.mark.parametrize("values", [
    [-8., -6., 5., 7.],
    [0., 1., 2., 2., 5.],
    [10., 10., 10., 10.],
])
def test_median_is_exactly_legacy_rule(values):
    z = np.array(values, float)
    np.testing.assert_array_equal(median_leaflet_mask(z), z >= np.median(z))


def test_median_rule_is_not_periodic_safe():
    # Intentionally document different semantics, not a replacement for
    # circular clustering. Bilayer spans the periodic boundary.
    z = np.array([98., 99., 1., 2., 48., 49., 51., 52.])
    raw = median_leaflet_mask(z)
    periodic = classify_leaflets(z, 100.)
    assert not np.array_equal(raw, periodic)
    assert np.count_nonzero(raw) == 4


def test_bad_coordinates_fail():
    with pytest.raises(ValueError, match="finite"):
        median_leaflet_mask([1., np.nan])
