"""Scientific MSD regression tests: PBC, drift, batched COM and source frames."""
from types import SimpleNamespace

import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader

from mdscope.analysis.observables.msd import _residue_centers, initialize_msd, msd_metrics


def _toy_universe():
    universe = mda.Universe.empty(
        4, n_residues=4, atom_resindex=np.arange(4), trajectory=True,
    )
    universe.add_TopologyAttr("names", ["P"] * 4)
    universe.add_TopologyAttr("resnames", ["POPC"] * 4)
    universe.add_TopologyAttr("resids", [1, 2, 3, 4])
    universe.add_TopologyAttr("masses", [31.] * 4)
    # All four lipids translate by +2 Å in x; two cross the 10 Å PBC.
    frames = np.asarray([
        [[9, 1, 8], [3, 1, 8], [9, 3, 2], [3, 3, 2]],
        [[1, 1, 8], [5, 1, 8], [1, 3, 2], [5, 3, 2]],
    ], np.float32)
    boxes = np.tile([10., 10., 10., 90., 90., 90.], (len(frames), 1))
    universe.load_new(frames, format=MemoryReader, dimensions=boxes)
    return universe


def _state():
    return SimpleNamespace(ai=2, upper_resix=np.asarray([0, 1]),
                           lower_resix=np.asarray([2, 3]))


def test_msd_uses_raw_minimum_image_lipid_displacements():
    universe = _toy_universe()
    state = _state()
    initialize_msd(universe, state, {
        "headgroup_selection": "all",
        "msd": {"enabled": True, "include_3d": True},
    })
    raw_before = universe.atoms.positions.copy()
    assert msd_metrics(universe, state)["membrane_msd_lateral_A2"] == 0.
    universe.trajectory[1]
    result = msd_metrics(universe, state)
    np.testing.assert_allclose(result["membrane_msd_lateral_A2"], 4., atol=1e-6)
    np.testing.assert_allclose(result["membrane_msd_lateral_upper_A2"], 4., atol=1e-6)
    np.testing.assert_allclose(result["membrane_msd_lateral_lower_A2"], 4., atol=1e-6)
    np.testing.assert_allclose(result["membrane_msd_3d_A2"], 4., atol=1e-6)
    # MSD must not overwrite the trajectory with aligned or unwrapped coordinates.
    np.testing.assert_array_equal(universe.atoms.positions, universe.trajectory.ts.positions)
    assert raw_before[0, 0] == 9.


def test_drift_correction_removes_common_translation_without_modifying_source():
    universe = _toy_universe()
    state = _state()
    initialize_msd(universe, state, {
        "headgroup_selection": "all",
        "msd": {"enabled": True, "remove_drift": True},
    })
    universe.trajectory[1]
    before = universe.atoms.positions.copy()
    result = msd_metrics(universe, state)
    np.testing.assert_allclose(result["membrane_msd_lateral_A2"], 0., atol=1e-10)
    np.testing.assert_array_equal(universe.atoms.positions, before)


def test_batched_center_of_mass_equals_per_residue_reference_and_falls_back():
    universe = _toy_universe()
    reference = np.stack([res.atoms.center_of_mass() for res in universe.residues])
    np.testing.assert_allclose(_residue_centers(universe.atoms), reference)
    # Missing masses must still allow center-of-geometry fallback.
    universe = mda.Universe.empty(
        4, n_residues=2, atom_resindex=[0, 0, 1, 1], trajectory=True,
    )
    universe.atoms.positions = np.asarray(
        [[0., 0., 0.], [2., 0., 0.], [1., 0., 0.], [5., 0., 0.]],
        np.float32,
    )
    np.testing.assert_allclose(
        _residue_centers(universe.atoms),
        [[1., 0., 0.], [3., 0., 0.]],
    )


def test_drift_removal_does_not_eliminate_relative_diffusion():
    universe = _toy_universe()
    state = _state()
    initialize_msd(universe, state, {
        "headgroup_selection": "all",
        "msd": {"enabled": True, "remove_drift": True},
    })
    universe.trajectory[1]
    # One lipid has one extra Å of relative lateral movement.
    universe.atoms[0].position = np.asarray([2., 1., 8.], np.float32)
    result = msd_metrics(universe, state)
    # 4 lipids: drifts 2.25 Å; corrected steps (0.75,-0.25,-0.25,-0.25).
    np.testing.assert_allclose(result["membrane_msd_lateral_A2"], 0.1875, atol=1e-8)
