"""Coordinate-space contract tests; no changes to numerical estimators."""
import pytest

from mdscope.core.coordinate_contracts import (
    CoordinateSpace, require_coordinate_space, msd_source_space,
)


def test_source_roles_are_distinct_from_fitted_arrays():
    require_coordinate_space("membrane", CoordinateSpace.SOURCE, CoordinateSpace.SOURCE)
    with pytest.raises(ValueError, match="membrane requires source"):
        require_coordinate_space("membrane", CoordinateSpace.FITTED, CoordinateSpace.SOURCE)


def test_msd_reader_provenance():
    assert msd_source_space(False) is CoordinateSpace.SOURCE
    assert msd_source_space(True) is CoordinateSpace.MSD_SOURCE
    with pytest.raises(ValueError, match="MSD requires msd_source"):
        require_coordinate_space("MSD", CoordinateSpace.FITTED, CoordinateSpace.MSD_SOURCE)


def test_space_labels_are_stable_strings():
    assert CoordinateSpace.FITTED.value == "fitted"
