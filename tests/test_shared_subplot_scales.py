"""Production subplot scaling: consistent units mean comparable visual limits."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from mdscope.reporting.shared_scales import synchronize_y_axes


def test_shared_axis_scales_two_unequal_panels_and_ignores_empty():
    fig, axes = plt.subplots(1, 3)
    axes[0].plot([0, 1], [1, 2])
    axes[1].plot([0, 1], [100, 200])
    axes[2].text(.5, .5, "No finite data")
    bounds = synchronize_y_axes(axes)
    assert bounds[0] < 1 and bounds[1] > 200
    assert all(ax.get_ylim() == pytest.approx(bounds) for ax in axes)
    plt.close(fig)


def test_shared_axis_symmetric_difference_and_zero_constant():
    fig, axes = plt.subplots(1, 2)
    axes[0].plot([0, 1], [-4, -2])
    axes[1].plot([0, 1], [0, 1])
    low, high = synchronize_y_axes(axes, symmetric=True)
    assert -low == pytest.approx(high)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2)
    axes[0].plot([0, 1], [0, 0])
    axes[1].plot([0, 1], [0, 0])
    low, high = synchronize_y_axes(axes)
    assert low < 0 < high
    plt.close(fig)


def test_no_measured_data_does_not_create_limits():
    fig, axes = plt.subplots(1, 2)
    axes[0].text(.5, .5, "not measured")
    assert synchronize_y_axes(axes) is None
    plt.close(fig)


def _plot_config():
    return {
        "plots": {"enabled": True, "distributions": False,
                  "autocorrelation": False, "bootstrap_tables": False},
        "systems": {
            "control": {"replicates": [{"label": "rep1"}]},
            "treated": {"replicates": [{"label": "rep1"}]},
        },
    }


def test_time_series_production_panels_same_limits_after_smoothing(tmp_path, monkeypatch):
    from pathlib import Path

    import mdscope.reporting.plotting as plotting

    input_values = {
        "control": [1., 2., 3., 4.],
        "treated": [100., 120., 140., 160.],
    }
    def collect(path, **kwargs):
        values = input_values[Path(path).parent.name]
        return {(None, "test_metric"): pd.DataFrame({
            "time_ps": [0., 10_000., 20_000., 30_000.],
            "value": values,
        })}
    monkeypatch.setattr(plotting, "collect_time_series", collect)
    recorded = []
    def check(fig, path, cfg, plot_id):
        if plot_id.endswith(":timeseries"):
            recorded.append([ax.get_ylim() for ax in fig.axes])
        plt.close(fig)
    monkeypatch.setattr(plotting, "_save_selected", check)
    cfg = _plot_config()
    cfg["plots"]["time_series_rolling_window_ns"] = 100
    reps = {name: [tmp_path / name / "rep1"] for name in input_values}
    plotting.plot_time_series(cfg, tmp_path, reps)
    assert len(recorded) == 1
    assert recorded[0][0] == pytest.approx(recorded[0][1])
    assert recorded[0][0][0] < 1
    assert recorded[0][0][1] > 100


def test_profile_production_panels_include_uncertainty_extent(tmp_path, monkeypatch):
    import mdscope.reporting.plotting as plotting

    folder = tmp_path / "project_summary/statistics/profiles"
    folder.mkdir(parents=True)
    pd.DataFrame([
        {"profile": "p", "system": name, "replicate": "rep1",
         "position": j, "value": value}
        for name, values in (("control", [1., 2.]),
                             ("treated", [10., 20.]))
        for j, value in enumerate(values)
    ]).to_parquet(folder / "replicate_profiles.parquet")
    pd.DataFrame([
        {"profile": "p", "system": name, "position": j,
         "mean": value, "mean_ci_low": value - 1,
         "mean_ci_high": value + 3}
        for name, values in (("control", [1., 2.]),
                             ("treated", [10., 20.]))
        for j, value in enumerate(values)
    ]).to_parquet(folder / "system_profiles.parquet")
    recorded = []
    def check(fig, path, cfg, plot_id):
        if plot_id.endswith(":profile"):
            recorded.append([ax.get_ylim() for ax in fig.axes])
        plt.close(fig)
    monkeypatch.setattr(plotting, "_save_selected", check)
    cfg = _plot_config()
    plotting.plot_profiles(cfg, tmp_path, {
        "control": [tmp_path / "control/rep1"],
        "treated": [tmp_path / "treated/rep1"],
    })
    assert len(recorded) == 1
    assert recorded[0][0] == pytest.approx(recorded[0][1])
    assert recorded[0][0][1] > 23


def test_solvation_rdf_production_panels_share_y_limits(tmp_path, monkeypatch):
    import mdscope.reporting.solvation_figures as solvation

    folder = tmp_path / "project_summary/statistics/solvation"
    folder.mkdir(parents=True)
    pd.DataFrame([
        dict(solvent="water", target_class=target, site_class="protein",
             system="control", r_inner_A=1., r_outer_A=2.,
             mean_g_bulk=y, sd_g_bulk=0.2, n_replicates=3)
        for target, y in (("acidic", 1.0), ("neutral", 8.0))
    ]).to_parquet(folder / "system_rdf.parquet")
    actual = plt.Figure.savefig
    ranges = []
    def capture(self, *args, **kwargs):
        if len(self.axes) >= 2:
            ranges.append([ax.get_ylim() for ax in self.axes])
        return actual(self, *args, **kwargs)
    monkeypatch.setattr(plt.Figure, "savefig", capture)
    result = solvation.solvation_rdf_figure(tmp_path, {})
    assert result is not None
    assert len(ranges) == 1
    assert ranges[0][0] == pytest.approx(ranges[0][1])


def test_distance_heatmaps_share_time_axis_and_distance_colorbar(tmp_path, monkeypatch):
    from types import SimpleNamespace

    import mdscope.config as config
    import mdscope.reporting.plotting as plotting

    run = SimpleNamespace(run_id="default", label="Test", spec={})
    monkeypatch.setattr(config, "analysis_runs", lambda cfg, analysis: [run])
    monkeypatch.setattr(
        config, "profile_dir",
        lambda repout, analysis, run_id: repout / "profiles" / analysis / run_id,
    )
    repouts = [tmp_path / "control" / f"rep{i}" for i in (1, 2)]
    for repout, last_time, peak in zip(repouts, (10_000., 20_000.), (2., 20.)):
        folder = repout / "profiles/residue_membrane_distance/default"
        folder.mkdir(parents=True)
        np.save(folder / "distance_matrix_A.npy",
                np.array([[1., peak], [1., peak]], dtype=np.float32))
        np.save(folder / "time_ps.npy", np.array([0., last_time]))
        pd.DataFrame({"residue": ["ALA1", "GLY2"]}).to_parquet(
            folder / "residues.parquet")

    observed = []
    def capture(fig, path, cfg):
        observed.append([(ax.get_ylim(), ax.images[0].get_clim()) for ax in fig.axes
                         if ax.images])
        plt.close(fig)
    monkeypatch.setattr(plotting, "_save", capture)
    cfg = {"plots": {"enabled": True},
           "systems": {"control": {"replicates": [
               {"label": "rep1"}, {"label": "rep2"}]}}}
    plotting._plot_residue_distance_heatmaps(
        cfg, tmp_path, {"control": repouts})
    assert len(observed) == 1 and len(observed[0]) == 2
    assert observed[0][0][0] == pytest.approx(observed[0][1][0])
    assert observed[0][0][1] == pytest.approx(observed[0][1][1])
    assert observed[0][0][1][1] == pytest.approx(20.)
