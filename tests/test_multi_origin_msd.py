"""Independent direct-sum oracle for multi-time-origin MSD."""
import numpy as np
import pytest

from mdscope.observables.multi_origin_msd import (
    multi_origin_msd, multi_origin_diffusion,
)


def brute_force(xyz, normal=2):
    lateral = [k for k in range(3) if k != normal]
    return np.array([
        np.mean(np.sum(
            (xyz[lag:][:, :, lateral] - xyz[:len(xyz)-lag][:, :, lateral]) ** 2,
            axis=-1,
        ))
        for lag in range(len(xyz))
    ])


@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("frames,lipids", [(7, 1), (17, 4), (32, 9)])
def test_fft_matches_independent_direct_sum(axis, frames, lipids):
    rng = np.random.default_rng(frames * 10 + lipids)
    xyz = np.cumsum(rng.normal(size=(frames, lipids, 3)), axis=0)
    xyz += 1e3
    lags, got = multi_origin_msd(xyz, normal_axis=axis, particle_batch=3)
    np.testing.assert_array_equal(lags, np.arange(frames))
    np.testing.assert_allclose(got, brute_force(xyz, axis), rtol=1e-9, atol=1e-7)


def test_deterministic_translation_exact_lag_square():
    t = np.arange(20, dtype=float)
    xyz = np.zeros((20, 3, 3))
    xyz[:, :, 0] = t[:, None] * 2.
    lags, result = multi_origin_msd(xyz)
    np.testing.assert_allclose(result, 4. * lags ** 2, atol=1e-9)


def test_independent_lag_fit_for_linear_msd():
    # The fitter should use lag times; synthetic Brownian trajectories only
    # converge statistically, so compare FFT and direct sum rather than
    # assuming a finite sample returns the true diffusion constant exactly.
    rng = np.random.default_rng(8)
    xyz = np.cumsum(rng.normal(size=(120, 12, 3)), axis=0)
    times = np.arange(120, dtype=float) * 5.
    result = multi_origin_diffusion(times, xyz, max_lag=35)
    oracle = brute_force(xyz)[:36]
    np.testing.assert_allclose(result["msd_lateral_A2"], oracle, atol=1e-9)
    mask = (result["lag_ps"] >= result["fit_start_ps"]) & (
        result["lag_ps"] <= result["fit_end_ps"])
    mask[0] = False
    expected_slope = np.polyfit(result["lag_ps"][mask], oracle[mask], 1)[0]
    np.testing.assert_allclose(result["diffusion_nm2_ns"], expected_slope * 2.5)


def test_invalid_inputs_fail_explicitly():
    xyz = np.zeros((10, 2, 3))
    with pytest.raises(ValueError, match="regularly spaced"):
        multi_origin_diffusion([0, 1, 2, 3, 4, 5, 6, 7, 8, 11], xyz)
    with pytest.raises(ValueError, match="Non-finite"):
        invalid = xyz.copy()
        invalid[0, 0, 0] = np.nan
        multi_origin_msd(invalid)
    with pytest.raises(ValueError, match="max_lag"):
        multi_origin_msd(xyz, max_lag=10)


@pytest.mark.parametrize("exponent, expected", [(1., "within_screen"), (.55, "window_or_scaling_sensitive")])
def test_diffusion_fit_screen_distinguishes_scaling(exponent, expected):
    from mdscope.observables.multi_origin_msd import diffusion_fit_sensitivity

    lag = np.arange(101, dtype=float)
    rows = diffusion_fit_sensitivity(lag, 3. * lag**exponent, 101 - lag)
    assert {row["screen"] for row in rows} == {expected}
    assert [row["min_origins_observed"] for row in rows] == [21, 41, 21]
    np.testing.assert_allclose([row["log_log_slope"] for row in rows], exponent)
    if exponent == 1.:
        np.testing.assert_allclose([row["diffusion_nm2_ns"] for row in rows], 7.5)


def test_diffusion_fit_screen_flags_poor_support_and_negative_fit():
    from mdscope.observables.multi_origin_msd import diffusion_fit_sensitivity

    lag = np.arange(8, dtype=float)
    rows = diffusion_fit_sensitivity(lag, lag, 8 - lag)
    assert {row["fit_status"] for row in rows} == {"insufficient_support"}
    assert all(np.isnan(row["diffusion_nm2_ns"]) for row in rows)
    lag = np.arange(101, dtype=float)
    rows = diffusion_fit_sensitivity(lag, np.where(lag > 0, 101 - lag, 0), 101 - lag)
    assert {row["fit_status"] for row in rows} == {"nonpositive_fit"}
    assert {row["screen"] for row in rows} == {"unavailable"}


def test_project_msd_diagnostics_uses_saved_curves_and_preserves_default(tmp_path, monkeypatch):
    import MDAnalysis as mda
    import pandas as pd
    from mdscope.observables.msd import summarize_project_multi_origin

    def forbidden(*args, **kwargs):
        raise AssertionError("diagnostics reopened a trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)
    repout = tmp_path / "condition" / "rep1"
    summary = repout / "summary"
    summary.mkdir(parents=True)
    pd.DataFrame({"system": ["condition"], "replicate": ["rep1"], "metric": ["lateral"],
                  "diffusion_nm2_ns": [7.5]}).to_parquet(summary / "multi_origin_diffusion.parquet")
    lag = np.arange(101, dtype=float)
    pd.DataFrame({"metric": "lateral", "lag_ps": lag, "msd_A2": 3. * lag,
                  "n_time_origins": 101 - lag}).to_parquet(summary / "multi_origin_msd.parquet")
    reps = {"condition": [repout, tmp_path / "condition" / "rep2"]}
    summarize_project_multi_origin(tmp_path, reps)
    dest = tmp_path / "project_summary" / "diffusion"
    output = dest / "multi_origin_fit_diagnostics.parquet"
    assert not output.exists()
    previous = pd.read_parquet(dest / "multi_origin_replicates.parquet")
    cfg = {"analyses": {"membrane": {"msd": {"multi_origin_fit_diagnostics": True}}}}
    summarize_project_multi_origin(tmp_path, reps, cfg)
    rows = pd.read_parquet(output)
    assert rows.loc[rows.replicate == "rep1", "screen"].eq("within_screen").all()
    assert rows.loc[rows.replicate == "rep2", "fit_status"].item() == "missing_curve_or_fit"
    pd.testing.assert_frame_equal(previous, pd.read_parquet(dest / "multi_origin_replicates.parquet"))
    summarize_project_multi_origin(tmp_path, reps)
    assert not output.exists()
