"""Independent geometry and observable checks for the leaflet footprint."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from mdscope.observables.lipid_enrichment import footprint_bands, summarize_lipid_composition
from mdscope_test_support import (
    membrane_config as _cfg, membrane_universe as _universe, execute_membrane as _execute,
)


def test_periodic_footprint_uses_the_same_geometry_for_lipids_and_area():
    cell = np.array([100., 100., 100.])
    protein = np.array([[98., 15., 25.], [6., 15., 25.], [40., 30., 65.]])
    anchors = np.array([[0., 15., 25.], [2., 15., 25.], [6., 20., 25.], [40., 30., 25.]])
    fraction = (np.arange(200) + .5) / 200
    grid = np.stack(np.meshgrid(fraction, fraction), axis=-1).reshape(-1, 2)
    band, areas, full = footprint_bands(anchors, protein, cell, 2, [0, 2, 5], 3., 6., grid)
    np.testing.assert_array_equal(band, [-1, 0, 1, 2])
    # Independent pointwise minimum-image disk-union oracle for this planar fixture.
    xy = grid * 100 + protein[0, :2]
    difference = xy[:, None] - protein[None, :2, :2]
    difference -= 100 * np.rint(difference / 100)
    distance = np.sqrt((difference**2).sum(axis=2)).min(axis=1) - 3
    np.testing.assert_array_equal(areas, [((distance >= lo) & (distance < hi)).sum() / 4
                                         for lo, hi in [(0, 2), (2, 5)]])
    assert full == (distance >= 0).sum() / 4
    assert full == pytest.approx(10000 - 2*np.pi*3**2, abs=2.)


def test_footprint_resume_absent_species_and_direct_contacts(tmp_path):
    cfg = deepcopy(_cfg())
    spec = cfg["analyses"]["lipid_enrichment"]
    spec.update(allow_absent_species=True, slab_half_width_A=25., contact_cutoff_A=24.)
    spec["species_selections"] = {**spec["species_selections"], "absent": "resname NONE"}
    for name, resume in (("clean", False), ("resume", True)):
        _execute(_universe(), cfg, tmp_path / name, resume=resume)
    relative = "extensions/lipid_footprint"
    left, right = [tmp_path / name / relative for name in ("clean", "resume")]
    counts = pd.read_parquet(left / "counts.parquet")
    pd.testing.assert_frame_equal(counts, pd.read_parquet(right / "counts.parquet"))
    assert counts[counts.species == "absent"].area_density_ratio.isna().all()
    series = pd.read_parquet(left / "timeseries.parquet")
    assert series.lipid_footprint_POPC_upper_contact_fraction.eq(.25).all()
    assert series.lipid_footprint_absent_upper_mean_min_distance_A.isna().all()
    distances = np.sqrt(np.array([9., 24., 9., 24.])**2 + np.array([5., 5., 10., 10.])**2 + 20.**2)
    assert series.lipid_footprint_POPC_upper_mean_min_distance_A.iloc[0] == pytest.approx(distances.mean())
    summarize_lipid_composition(tmp_path / "clean")
    inventory = pd.read_parquet(left / "composition_summary.parquet")
    assert inventory[inventory.species == "absent"].fraction_frames_available.eq(0).all()
    assert inventory[inventory.species == "POPC"].mean_molecules.eq(4).all()
