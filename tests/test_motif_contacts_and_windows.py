
import numpy as np
import pandas as pd

from mdscope.mechanisms.sdc4_mechanisms import _leaflet_plane_height
from mdscope.mechanisms.shared_states import _write_matched_populations
from mdscope.statistics import analyze_project


def test_hs_episodes_and_dual_scalar_windows(tmp_path):
    prefix = "sdc4_motif_proline_hydrophobic_A_"
    cfg = {"statistics": {"bootstrap_draws": 32, "block": {"length": 1}, "matched_window_ps": "auto"}}
    repouts = {}
    assignments = []
    for system, start in (("first", 0), ("second", 2)):
        repouts[system] = []
        for name in ("rep1", "rep2"):
            rep = tmp_path / system / name
            (rep / "summary").mkdir(parents=True)
            repouts[system].append(rep)
            times = np.arange(start, start + 10, dtype=float)
            contact = (times.astype(int) % 4 <= 1)
            pd.DataFrame({
                "time_ps": times, prefix + "glycan_contacts": np.where(contact, 3, 0),
                prefix + "hs_1_contact_state": contact.astype(int), "scalar": times,
            }).to_parquet(rep / "summary" / "timeseries.parquet")
            assignments.append(pd.DataFrame({
                "system": system, "replicate": name, "time_ps": times,
                "state": np.where(contact, 0, 1),
            }))

    analyze_project(cfg, tmp_path, repouts)
    stats = tmp_path / "project_summary" / "statistics"
    full = pd.read_parquet(stats / "time_series" / "replicate_statistics.parquet")
    matched = pd.read_parquet(stats / "matched_time_series" / "replicate_statistics.parquet")
    np.testing.assert_allclose(full.loc[full.metric.eq("scalar"), "mean"], [4.5, 4.5, 6.5, 6.5])
    np.testing.assert_allclose(matched.loc[matched.metric.eq("scalar"), "mean"], 5.5)
    np.testing.assert_allclose(matched.loc[matched.metric.eq(prefix + "hs_contacts_during_contact"), "mean"], 3)
    episodes = pd.read_parquet(stats / "matched_time_series" / "hs_contact_episodes.parquet")
    assert set(episodes.observable) == {prefix + "hs_any_contact_state", prefix + "hs_1_contact_state"}
    assert set(episodes.start_ps) == {4, 8}
    assert episodes.loc[episodes.start_ps.eq(8), "right_censored"].all()

    out = tmp_path / "shared_states"
    _write_matched_populations(cfg, out, pd.concat(assignments), repouts, 2)
    full = pd.read_parquet(out / "full_population_comparisons.parquet")
    matched = pd.read_parquet(out / "matched_population_comparisons.parquet")
    assert abs(full.difference_system2_minus_system1).max() > 0
    np.testing.assert_allclose(matched.difference_system2_minus_system1, 0)


def test_local_phosphate_plane_signed_height_and_periodic_xy():
    position = np.array([0.5, 0.5, 8.0])
    x = np.array([-2., -2., 2., 2.])
    y = np.array([-2., 2., -2., 2.])
    anchors = np.column_stack(((position[0] + x) % 20, (position[1] + y) % 20, 5 + .2 * x - .1 * y))
    box = np.array([20., 20., 30., 90., 90., 90.])
    np.testing.assert_allclose(_leaflet_plane_height(position, anchors, box, 1, 4, 4), 3)
    assert np.isnan(_leaflet_plane_height(position, anchors, box, 1, 2, 4))
