"""Regressions for provenance, merged outputs, and lipid PBC geometry."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import MDAnalysis as mda

from mdscope.pipeline import runner
from mdscope.core.lipid_geometry import AnchorTracker


def test_merge_rejects_conflicting_observable_names():
    first = pd.DataFrame({"time_ps": [0., 1.], "msd": [0., 1.]})
    second = pd.DataFrame({"time_ps": [0., 1.], "msd": [0., 2.]})
    with pytest.raises(ValueError, match="Conflicting observable"):
        runner._merge_timeseries([first, second])


def test_merge_rejects_duplicate_times():
    table = pd.DataFrame({"time_ps": [1., 1.], "count": [0, 1]})
    with pytest.raises(ValueError, match="Duplicate time-series keys"):
        runner._merge_timeseries([table])


def test_merge_preserves_missing_counts():
    a = pd.DataFrame({"time_ps": [0., 1.], "contact_count": [2., 3.]})
    b = pd.DataFrame({"time_ps": [1., 2.], "other_count": [4., 5.]})
    merged = runner._merge_timeseries([a, b])
    assert np.isnan(merged.loc[merged.time_ps == 0., "other_count"].iloc[0])
    assert np.isnan(merged.loc[merged.time_ps == 2., "contact_count"].iloc[0])


def test_multiatom_anchor_straddling_pbc():
    u = mda.Universe.empty(
        4, n_residues=2, atom_resindex=[0, 0, 1, 1], trajectory=True
    )
    u.dimensions = [100., 100., 100., 90., 90., 90.]
    u.atoms.positions = np.array([
        [99., 4., 5.], [1., 4., 5.],
        [30., 4., 5.], [32., 4., 5.],
    ], dtype=np.float32)
    tracker = AnchorTracker(u.atoms, np.array([0, 0, 1, 1]), np.array([2, 2]))
    pos = tracker.positions()
    assert min(abs(pos[0, 0]), abs(pos[0, 0] - 100.)) < 1e-5
    np.testing.assert_allclose(pos[1, 0], 31.)
