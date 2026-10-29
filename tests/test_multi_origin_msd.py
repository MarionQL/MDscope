"""Independent direct-sum oracle for multi-time-origin MSD."""
import numpy as np
import pytest

from mdscope.analysis.observables.multi_origin_msd import (
    multi_origin_msd, multi_origin_diffusion,
)


def brute_force(xyz, normal=2):
    lateral = [k for k in range(3) if k != normal]
    return np.array([
        np.mean(np.sum(
            (xyz[lag:, :, lateral] - xyz[:len(xyz)-lag, :, lateral]) ** 2,
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
