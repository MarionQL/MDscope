"""Time-based display smoothing must leave source observations and statistics intact."""
import numpy as np
import pandas as pd
import pytest

from copy import deepcopy
from mdscope.config import DEFAULTS
from mdscope.reporting.plotting import rolling_time_mean


def test_trailing_window_uses_timestamps_not_number_of_frames():
    # 100 ns window: last point includes only 150 and 200 ns, not 0 ns.
    time_ps = [0., 20_000., 150_000., 200_000.]
    values = [10., 20., 30., 50.]
    result = rolling_time_mean(time_ps, values, 100)
    np.testing.assert_allclose(result, [10., 15., 30., 40.])


def test_missing_samples_excluded_but_measured_zero_is_retained():
    times = [0., 10_000., 20_000., 200_000., 210_000.]
    values = [np.nan, 0., 20., np.nan, np.nan]
    original = np.asarray(values).copy()
    result = rolling_time_mean(times, values, 100)
    assert np.isnan(result[0])
    np.testing.assert_allclose(result[1:3], [0., 10.])
    assert np.isnan(result[3:]).all()
    np.testing.assert_array_equal(np.isnan(values), np.isnan(original))


def test_rolling_window_requires_valid_time_grid():
    with pytest.raises(ValueError, match="strictly increasing"):
        rolling_time_mean([0., 0.], [1., 2.], 100)
    with pytest.raises(ValueError, match="one-dimensional"):
        rolling_time_mean([0., 10.], [1.], 100)


def _valid_cfg():
    cfg = deepcopy(DEFAULTS)
    cfg["systems"] = {"demo": {
        "topology": "topology.pdb",
        "replicates": [{"label": "rep1", "trajectory": "trajectory.xtc"}],
    }}
    return cfg


@pytest.mark.parametrize("window", [0, -1, float("inf"), float("nan"), True, "100"])
def test_plot_rolling_window_rejects_invalid_config(window):
    cfg = _valid_cfg()
    cfg["plots"]["time_series_rolling_window_ns"] = window
    # Config validation should reject the display option before opening files.
    from mdscope.config import validate_config
    with pytest.raises(ValueError, match="plots.time_series_rolling_window_ns"):
        validate_config(cfg, check_paths=False)


def test_window_config_accepts_null_and_100_ns():
    from mdscope.config import validate_config
    for window in (None, 100, 12.5):
        cfg = _valid_cfg()
        cfg["plots"]["time_series_rolling_window_ns"] = window
        validate_config(cfg, check_paths=False)


def test_renderer_uses_rolling_mean_only_for_trace(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt
    import mdscope.reporting.plotting as plotting

    source = pd.DataFrame({
        "time_ps": [0., 20_000., 150_000., 200_000.],
        "value": [10., 20., 30., 50.],
    })
    monkeypatch.setattr(plotting, "collect_time_series", lambda *a, **k: {
        (None, "test_metric"): source,
    })
    observed = []

    def inspect(fig, _path, _cfg, _id):
        observed.extend(line.get_ydata().copy()
                        for ax in fig.axes for line in ax.lines)
        plt.close(fig)

    monkeypatch.setattr(plotting, "_save_selected", inspect)
    cfg = {
        "plots": {
            "enabled": True, "time_series": True,
            "time_series_rolling_window_ns": 100,
            "distributions": False, "autocorrelation": False,
            "bootstrap_tables": False,
        },
        "systems": {
            "control": {"label": "Control",
                        "replicates": [{"label": "rep1", "trajectory": "unused"}]},
        },
    }
    plotting.plot_time_series(cfg, tmp_path, {"control": [tmp_path / "control/rep1"]})
    assert len(observed) == 1
    np.testing.assert_allclose(observed[0], [10., 15., 30., 40.])
    np.testing.assert_array_equal(source["value"].to_numpy(), [10., 20., 30., 50.])
