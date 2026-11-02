"""Independent, interpretable numerical oracles for MDscope physical contracts.

These tests use hand-countable configurations and formulas that do not call
the implementation to compute their expected answers. Passing them does not
establish converged simulation estimates or transferable binding affinities.
"""
import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda

from mdscope.mechanisms.motif_residence import (
    summarize_lipid_residence, summarize_intermittent_residence,
)
from mdscope.mechanisms.sdc4_mechanisms import (
    _composition_enrichment, _contacts, _accessible_area_grid,
)


def _events(times):
    return pd.DataFrame({
        "time_ps": np.asarray(times, float),
        "motif": ["NxIP"] * len(times),
        "molecule_id": [7] * len(times),
        "species": ["GM1a"] * len(times),
    })


def test_composition_reference_two_species_and_missing_denominator():
    # One GM1a versus three POPC molecules in the sampled leaflet.
    # Both have one contacting molecule. Contact shares are 1/2 each;
    # abundance shares are 1/4 and 3/4 respectively.
    table = pd.DataFrame({
        "time_ps": [0.] * 5,
        "region": ["TM_A"] * 5,
        "leaflet": ["upper"] * 5,
        "species": ["GM1a", "POPC", "GM1a", "BGAL", "BGAL"],
        "moiety": ["any", "any", "glycan", "any", "ceramide_tail"],
        "contacted_molecules": [1, 1, 1, 0, 0],
        "available_molecules": [1, 3, 1, 0, 0],
        "locally_available_molecules": [1, 3, 1, 0, 0],
    })
    for column in ("available_molecules", "locally_available_molecules"):
        result = _composition_enrichment(table, column)
        assert result[0] == pytest.approx(2.0)
        assert result[1] == pytest.approx(2 / 3)
        assert np.isnan(result[2])  # no double-counting moiety subclasses
        assert np.isnan(result[3])  # BGAL unavailable, not zero-enrichment
        assert np.isnan(result[4])


def test_local_population_mismatch_reference_is_detectable():
    # Contacting GM1a is beyond the radius used for local anchors.
    # This is a deliberately diagnostic counterexample: 1 contact / 0
    # locally available is not a physically valid local occupancy.
    table = pd.DataFrame({
        "time_ps": [0., 0.], "region": ["TM_A"] * 2,
        "leaflet": ["upper"] * 2, "species": ["GM1a", "POPC"],
        "moiety": ["any", "any"], "contacted_molecules": [1, 1],
        "available_molecules": [1, 3],
        "locally_available_molecules": [0, 3],
    })
    result = _composition_enrichment(table, "locally_available_molecules")
    assert np.isnan(result[0])
    assert result[1] == pytest.approx(.5)  # contact share 1/2, local share 3/3
    assert table.loc[0, "contacted_molecules"] > table.loc[0, "locally_available_molecules"]
    # The oracle exposes a population mismatch; it does not silently
    # claim local enrichment for the unavailable GM1a.


def test_orthorhombic_minimum_image_contact_oracle():
    # Across a 10 Å periodic boundary: 9.7 -> 0.2 is 0.5 Å,
    # while 9.7 -> 4.2 is 4.5 Å. Test strict cutoff behavior too.
    class Atoms:
        def __init__(self, coords):
            self.positions = np.asarray(coords, dtype=np.float32)
        def __len__(self):
            return len(self.positions)
    box = np.array([10., 10., 10., 90., 90., 90.])
    left = Atoms([[9.7, 5., 5.]])
    right = Atoms([[.2, 5., 5.], [4.2, 5., 5.]])
    actual = _contacts(left, right, 1., box)
    assert actual.tolist() == [[0, 0]]
    assert _contacts(left, right, .4, box).size == 0
    # Rigid translation by one box length must preserve contacts.
    translated = Atoms([[19.7, 5., 5.]])
    np.testing.assert_array_equal(_contacts(translated, right, 1., box), actual)


def test_projected_area_grid_periodicity_and_slab_oracle():
    grid = np.array([[x, y] for x in (-1., 0., 1.)
                     for y in (-1., 0., 1.)], float)
    box = np.array([10., 10., 10., 90., 90., 90.])
    # One protein atom at periodic x=-0.2 blocks only (0,0) on a
    # 1 Å grid when exclusion radius is 0.4 Å.
    wrapped = _accessible_area_grid([0., 0.], 0., [[9.8, 0., 0.]],
                                    box, grid, 1., .4, 1.)
    direct = _accessible_area_grid([0., 0.], 0., [[-.2, 0., 0.]],
                                   box, grid, 1., .4, 1.)
    outside_slab = _accessible_area_grid([0., 0.], 0., [[9.8, 0., 2.]],
                                         box, grid, 1., .4, 1.)
    assert wrapped == pytest.approx(8.)
    assert direct == pytest.approx(8.)
    assert outside_slab == pytest.approx(9.)


def test_residence_stride_and_gap_tolerance_are_observation_dependent():
    # Same hypothetical event: contact at 0, 2 and 4 ps, absent
    # at 1 and 3. A 2-ps stride cannot detect the two missing frames.
    dense_times = np.arange(5., dtype=float)
    sparse_times = np.array([0., 2., 4.])
    events = _events([0., 2., 4.])
    dense = summarize_lipid_residence(events, dense_times)
    sparse = summarize_lipid_residence(events, sparse_times)
    assert dense.n_frames.tolist() == [1, 1, 1]
    assert sparse.n_frames.tolist() == [3]
    assert dense.occupancy_fraction.tolist() == pytest.approx([3 / 5] * 3)
    assert sparse.occupancy_fraction.tolist() == pytest.approx([1.])
    bridged = summarize_intermittent_residence(events, dense_times, gap_frames=1)
    assert bridged.n_frames.tolist() == [3]
    assert bridged.duration_ps.tolist() == pytest.approx([5.])


def test_residence_censoring_and_irregular_last_interval_oracle():
    times = np.array([0., 3., 7., 12.])
    events = _events([0., 3., 12.])
    result = summarize_lipid_residence(events, times)
    assert len(result) == 2
    assert result.left_censored.tolist() == [True, False]
    assert result.right_censored.tolist() == [False, True]
    assert result.n_frames.tolist() == [2, 1]
    # The first dwell includes its next sampled interval (3->7).
    assert result.duration_ps.iloc[0] == pytest.approx(7.)
    # The terminal episode has no observed next interval. Current
    # convention extrapolates the median spacing (3,4,5)->4.
    assert result.duration_ps.iloc[1] == pytest.approx(4.)
    assert result.occupancy_fraction.tolist() == pytest.approx([.75, .75])


def test_ceramide_head_tail_glycan_chemistry_selection_oracle():
    # A synthetic topology has no force-field chemistry inference:
    # label identities and atom names are explicit, checkable inputs.
    u = mda.Universe.empty(
        7, n_residues=3,
        atom_resindex=[0, 0, 0, 0, 1, 1, 2],
        trajectory=True,
    )
    u.add_TopologyAttr("names", ["NF", "O3", "C4S", "C2F",
                                   "C1", "O5", "O1"])
    u.add_TopologyAttr("resnames", ["CER160", "BGLCA", "CER160"])
    u.add_TopologyAttr("chainIDs", ["GM1a"] * 6 + ["CER160"])
    head = u.select_atoms(
        "chainID GM1a CER160 and resname CER160 and name NF O3 O1"
    )
    tail = u.select_atoms(
        "chainID GM1a CER160 and resname CER160 and name C4S C2F"
    )
    glycan = u.select_atoms(
        "chainID GM1a and not resname CER160 and not name H*"
    )
    assert set(head.indices.tolist()) == {0, 1, 6}
    assert set(tail.indices.tolist()) == {2, 3}
    assert set(glycan.indices.tolist()) == {4, 5}
    assert not (set(head.indices) & set(tail.indices))
    assert not (set(glycan.indices) & (set(head.indices) | set(tail.indices)))
    assert "O1" in set(head.names)  # free ceramide hydroxyl in control
