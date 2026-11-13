"""Per-sugar hydration uses bond order, periodic distances and unique waters."""
from types import SimpleNamespace

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from mdscope.mechanisms import hs_hydration


def test_sugar_hydration_with_reversed_residue_ids_and_shared_water(tmp_path, monkeypatch):
    u = mda.Universe.empty(7, n_residues=5, atom_resindex=[0, 1, 1, 2, 3, 4, 4],
                           trajectory=True)
    u.add_TopologyAttr("names", ["OG", "O1", "O2", "O3", "O4", "OH2", "H1"])
    u.add_TopologyAttr("resnames", ["SER", "BXYL", "BGAL", "BGAL", "TIP3"])
    u.add_TopologyAttr("resids", [39, 99, 40, 15, 102])
    u.add_TopologyAttr("chainIDs", ["PROA"] * 7)
    u.add_bonds([(0, 1), (1, 2), (2, 3), (3, 4), (5, 6)])
    u.atoms.positions = [[1, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0],
                         [8, 0, 0], [19, 0, 0], [18, 0, 0]]
    u.dimensions = [20, 20, 20, 90, 90, 90]
    cfg = {
        "trajectory": {"stride": 1},
        "time": {"dt_ps": 1.},
        "analyses": {"sdc4_mechanisms": {
            "hs_sugar_water_cutoff_A": 4.5, "water_oxygen_selection": "name OH2",
            "hs_auto": {"glycan_selection": "resname BXYL BGAL",
                        "attachment_selection": "name OG"},
        }},
    }
    monkeypatch.setattr(hs_hydration.mda, "Universe", lambda *_: u)
    hs_hydration.run_hs_hydration(SimpleNamespace(topology="top", trajectory="traj"),
                                  cfg, tmp_path, np.array([0.]))
    data = pd.read_parquet(tmp_path / hs_hydration.OUTPUT)
    assert list(data.sugar_resid) == [99, 40, 15]
    assert list(data.water_count) == [1, 1, 0]
    profiles = hs_hydration.temporal_hs_hydration_profiles(tmp_path, cfg, np.array([0.]))
    by_position = profiles["hs_water_by_position__PROA_S39"]
    by_type = profiles["hs_water_by_sugar_type__PROA_S39"]
    assert list(by_position["labels"]) == ["1:BXYL", "2:BGAL", "3:BGAL"]
    assert list(by_type["labels"]) == ["BGAL (2,3)", "BXYL (1)"]
    np.testing.assert_array_equal(by_type["values"], [[.5, 1.]])


@pytest.mark.parametrize("cfg", [
    {},
    {"statistics": {"bootstrap_draws": 32}},
    {"analyses": {"sdc4_mechanisms": {"hs_sugar_water_cutoff_A": None}}},
])
def test_optional_hs_hydration_is_ignored_without_a_cutoff(tmp_path, cfg):
    # Other profile statistics must not require the optional hydration analysis.
    assert hs_hydration.temporal_hs_hydration_profiles(tmp_path, cfg) == {}


def test_configured_hs_hydration_requires_saved_data(tmp_path):
    cfg = {"analyses": {"sdc4_mechanisms": {"hs_sugar_water_cutoff_A": 4.5}}}
    with pytest.raises(ValueError, match="Missing HS sugar hydration"):
        hs_hydration.temporal_hs_hydration_profiles(tmp_path, cfg)
