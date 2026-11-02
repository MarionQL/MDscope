"""Sampling conventions shared by lipid, water, and ion contact episodes."""
import numpy as np


def sampled_episode_durations(times, first_index, last_index):
    """Return observed spans and interval-imputed durations, in ps.

    Indices may be scalars or equally shaped integer arrays. Evaluate all
    episodes together so time validation and the terminal median run once.
    Imputation adds the next interval, or the median interval at the right
    edge. A single observation has no estimable interval (NaN duration).
    Neither return value is a continuously observed microscopic dwell time.
    """
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or not len(times) or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("Sampled times must be finite and increase strictly")
    first, last = np.asarray(first_index), np.asarray(last_index)
    if (first.shape != last.shape or first.dtype.kind not in "iu"
            or last.dtype.kind not in "iu" or np.any(first < 0)
            or np.any(first > last) or np.any(last >= len(times))):
        raise ValueError("Invalid sampled episode time indices")
    span = times[last] - times[first]
    intervals = np.diff(times)
    following = np.r_[intervals, np.median(intervals) if len(intervals) else np.nan]
    return span, span + following[last]
