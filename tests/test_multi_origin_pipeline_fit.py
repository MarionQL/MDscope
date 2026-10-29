"""Independent numerical checks for the opt-in pipeline's FFT and lag fits."""
import numpy as np
import pytest
from mdscope.analysis.observables.multi_origin_msd import (
    multi_origin_msd, fit_multi_origin_curve,
)


def reference_multi_origin(xyz):
    out = []
    for lag in range(len(xyz)):
        shift = xyz[lag:, :, :2] - xyz[:len(xyz)-lag, :, :2]
        out.append(np.mean(np.sum(shift**2, axis=-1)))
    return np.array(out)


def test_lag_curve_and_batching_against_direct_oracle():
    rng = np.random.default_rng(123)
    unwrapped = np.cumsum(rng.normal(size=(72, 9, 3)), axis=0)
    expected = reference_multi_origin(unwrapped)
    lag, got = multi_origin_msd(unwrapped, particle_batch=3)
    np.testing.assert_array_equal(lag, np.arange(72))
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-8)
    fit = fit_multi_origin_curve(lag * 10., got)
    window = ((lag * 10.) >= fit["fit_start_ps"]) & (
        (lag * 10.) <= fit["fit_end_ps"]) & (lag > 0)
    slope = np.polyfit((lag * 10.)[window], expected[window], 1)[0]
    np.testing.assert_allclose(fit["diffusion_nm2_ns"], slope * 2.5)


def test_invalid_fit_window_fails():
    with pytest.raises(ValueError, match="at least four"):
        fit_multi_origin_curve(np.arange(8.), np.arange(8.), fit_start_fraction=.9)
