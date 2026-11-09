"""Independent precision and end-to-end scalar hierarchy references."""

import numpy as np
import pandas as pd
import pytest

from mdscope.statistics.common import _adaptive_intervals
from mdscope.statistics.time_series import analyze_time_series


OPTIONS = dict(min_draws=800, max_draws=2400, batch_size=400, stable_checks=2,
               absolute_tolerance=0.01)


def test_adaptive_endpoints_require_repeated_joint_precision():
    draws, meta = _adaptive_intervals(lambda n: np.ones((n, 2)), .95, OPTIONS)
    assert draws.shape == (2, 1200)
    assert meta["bootstrap_mc_status"] == "converged"
    assert meta["mean_ci_low_mc_error_bound"] == 0
    assert meta["median_ci_high_mc_error_bound"] == 0

    # One constant statistic cannot stop a row with a rare second-statistic tail.
    rng = np.random.default_rng(813)
    draws, meta = _adaptive_intervals(
        lambda n: np.column_stack((np.ones(n), 20 * (rng.random(n) < .025))),
        .95, OPTIONS)
    assert draws.shape == (2, 2400)
    assert meta["bootstrap_mc_status"] == "not_converged"
    assert meta["median_ci_high_mc_error_bound"] == 20


@pytest.mark.parametrize("design", ["paired", "independent"])
def test_adaptive_time_series_output_and_unavailable_sampling(tmp_path, design):
    systems = {}
    for name, offset in (("a", 0), ("b", 5)):
        systems[name] = []
        for rep in ("r1", "r2"):
            directory = tmp_path / name / rep
            (directory / "summary").mkdir(parents=True)
            pd.DataFrame({"time_ps": np.arange(8.), "value": offset,
                          "gapped": [offset] * 7 + [np.nan]}).to_parquet(
                              directory / "summary/timeseries.parquet")
            systems[name].append(directory)

    cfg = {"statistics": {"seed": 13, "bootstrap_draws": 100,
                          "hypothesis_tests": {"design": design},
                          "families": {"time_series": {"adaptive_hierarchy": OPTIONS}}}}
    analyze_time_series(cfg, tmp_path, systems)
    path = tmp_path / "project_summary/statistics/time_series/system_comparisons.parquet"
    rows = pd.read_parquet(path).set_index("metric")
    assert rows.loc["value", "mean_difference_ci_low"] == 5
    assert rows.loc["value", "median_difference_ci_high"] == 5
    assert rows.loc["value", "bootstrap_draws"] == 1200
    assert rows.loc["value", "bootstrap_mc_status"] == "converged"
    assert rows.loc["gapped", "bootstrap_draws"] == 0
    assert rows.loc["gapped", "bootstrap_mc_status"] == "unavailable"
    analyze_time_series(cfg, tmp_path, systems)
    pd.testing.assert_frame_equal(rows, pd.read_parquet(path).set_index("metric"))


def test_skewed_draws_match_analytic_endpoint_precision():
    rng = np.random.default_rng(941)
    options = {**OPTIONS, "max_draws": 16000, "batch_size": 800, "absolute_tolerance": .4}
    draws, meta = _adaptive_intervals(
        lambda n: rng.exponential(size=(n, 2)), .95, options)
    assert meta["bootstrap_mc_status"] == "converged"
    assert 1600 <= draws.shape[1] < 16000
    expected = -np.log1p(-np.array([.025, .975]))
    for stat, values in zip(("mean", "median"), draws):
        for side, q, target in zip(("low", "high"), (.025, .975), expected):
            assert abs(np.quantile(values, q) - target) <= meta[f"{stat}_ci_{side}_mc_error_bound"]


def test_temporary_precision_does_not_stop_after_one_look():
    transient = np.zeros((400, 2))
    transient[:30] = 10.
    batches = iter((np.zeros((800, 2)), transient, np.zeros((400, 2))))
    draws, meta = _adaptive_intervals(
        lambda n: next(batches), .95, {**OPTIONS, "max_draws": 1600})
    assert draws.shape == (2, 1600)
    assert meta["bootstrap_mc_status"] == "not_converged"
