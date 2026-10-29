"""Multi-time-origin lateral MSD for *unwrapped* lipid-center trajectories.

This numerical kernel deliberately does not unwrap coordinates or read trajectories.
Callers must reconstruct displacement with appropriate PBC, preserve consistent
lipid identity and supply regularly sampled times. Independent, quadratic-time
reference calculations are provided in the tests.
"""
from __future__ import annotations

import numpy as np


def multi_origin_msd(
    positions_A: np.ndarray,
    *,
    normal_axis: int = 2,
    particle_batch: int = 128,
    max_lag: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (lag_frames, lateral MSD in Å²), averaging all time origins/lipids.

    positions_A: shape (n_frames, n_lipids, 3), unwrapped coordinates in Å.
    FFT autocorrelation gives O(n_lipids*n_frames*log(n_frames)) rather than
    an O(n_frames²) lag/origin loop. Batches bound temporary memory.
    """
    xyz = np.asarray(positions_A, dtype=np.float64)
    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise ValueError("positions_A must have shape (frames, lipids, 3)")
    n, particles, _ = xyz.shape
    if n < 3 or particles < 1:
        raise ValueError("At least 3 frames and 1 lipid are required")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("Non-finite coordinates are not supported")
    if normal_axis not in (0, 1, 2):
        raise ValueError("normal_axis must be 0, 1 or 2")
    if particle_batch < 1:
        raise ValueError("particle_batch must be positive")
    maximum = n - 1 if max_lag is None else int(max_lag)
    if not 1 <= maximum < n:
        raise ValueError("max_lag must be between 1 and n_frames - 1")
    lateral = [axis for axis in range(3) if axis != normal_axis]
    sums = np.zeros(maximum + 1, dtype=np.float64)
    fft_length = 1 << ((2 * n - 1).bit_length())
    # Shifting the origin reduces catastrophic cancellation for coordinates
    # located far away from zero.
    for first in range(0, particles, particle_batch):
        p = xyz[:, first:first + particle_batch, lateral].copy()
        p -= p[0:1]
        q = np.sum(p * p, axis=2)
        prefix = np.vstack([np.zeros((1, q.shape[1])), np.cumsum(q, axis=0)])
        freq = np.fft.rfft(p, n=fft_length, axis=0)
        autocorr = np.fft.irfft(freq.conj() * freq, n=fft_length, axis=0)
        autocorr = np.sum(autocorr[:maximum + 1], axis=2)
        lags = np.arange(maximum + 1)
        squared = prefix[n - lags] + (prefix[n] - prefix[lags])
        sums += np.sum(squared - 2 * autocorr, axis=1)
    lags = np.arange(maximum + 1, dtype=int)
    msd = sums / ((n - lags) * particles)
    msd[0] = 0.0
    # Clamp only machine-precision round-off, never significant negative values.
    tolerance = 1e-10 * max(1.0, float(np.max(np.abs(msd))))
    if np.any(msd < -tolerance):
        raise ArithmeticError("FFT MSD has unphysical negative values")
    msd = np.maximum(msd, 0.0)
    return lags, msd


def multi_origin_diffusion(
    times_ps: np.ndarray,
    positions_A: np.ndarray,
    *,
    normal_axis: int = 2,
    max_lag: int | None = None,
    fit_start_fraction: float = 0.2,
    fit_end_fraction: float = 0.8,
) -> dict:
    """Fit Einstein lateral D from multi-origin MSD at regular sampling.

    Fits lag times, NOT absolute trajectory timestamps. Returns Å²/ps slope
    and nm²/ns D; estimates are descriptive and do not imply uncertainty.
    """
    t = np.asarray(times_ps, dtype=float)
    if t.ndim != 1 or t.shape[0] != len(positions_A) or not np.all(np.isfinite(t)):
        raise ValueError("times_ps must be a finite time per frame")
    delta = np.diff(t)
    if np.any(delta <= 0) or not np.allclose(delta, delta[0], rtol=1e-5, atol=1e-8):
        raise ValueError("Multi-origin MSD requires regularly spaced, increasing times")
    if not 0 <= fit_start_fraction < fit_end_fraction <= 1:
        raise ValueError("Require 0 <= fit_start_fraction < fit_end_fraction <= 1")
    lags, msd = multi_origin_msd(positions_A, normal_axis=normal_axis, max_lag=max_lag)
    lag_ps = lags * float(delta[0])
    low = lag_ps[-1] * fit_start_fraction
    high = lag_ps[-1] * fit_end_fraction
    mask = (lag_ps >= low) & (lag_ps <= high) & (lags > 0)
    if mask.sum() < 4:
        raise ValueError("Fit window must contain at least four nonzero lag points")
    slope, intercept = np.polyfit(lag_ps[mask], msd[mask], 1)
    prediction = slope * lag_ps[mask] + intercept
    residual = float(np.sum((msd[mask] - prediction) ** 2))
    total = float(np.sum((msd[mask] - msd[mask].mean()) ** 2))
    return {
        "method": "multi_origin_fft",
        "diffusion_nm2_ns": float(slope * 2.5),
        "slope_A2_ps": float(slope),
        "intercept_A2": float(intercept),
        "r_squared": 1 - residual / total if total > 0 else np.nan,
        "n_fit": int(mask.sum()),
        "fit_start_ps": float(low),
        "fit_end_ps": float(high),
        "lag_ps": lag_ps,
        "msd_lateral_A2": msd,
        "origins_per_lag": (len(t) - lags).astype(int),
    }
