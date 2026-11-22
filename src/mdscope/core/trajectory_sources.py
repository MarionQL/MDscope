"""Frame-synchronized iteration of a primary and an optional raw MSD trajectory."""
from __future__ import annotations

import numpy as np


def paired_frames(primary_frames, raw_msd_frames, *, tolerance_ps=0.1):
    """Advance both streams together and reject mismatched sample timestamps.

    The returned coordinates come from the primary (possibly whole/centered)
    analysis reader; the secondary MDAnalysis universe independently advances
    to the raw frame used by msd_metrics.
    """
    if raw_msd_frames is None:
        yield from primary_frames
        return
    if len(primary_frames) != len(raw_msd_frames):
        raise ValueError("Main and raw MSD trajectories have different sampled frame counts")
    for analysis_ts, msd_ts in zip(primary_frames, raw_msd_frames):
        t_analysis = float(analysis_ts.time)
        t_msd = float(msd_ts.time)
        if not np.isfinite(t_analysis) or not np.isfinite(t_msd) or abs(t_analysis - t_msd) > tolerance_ps:
            raise ValueError(
                "Main and raw MSD frames are not time synchronized: "
                f"analysis={t_analysis:g} ps, raw MSD={t_msd:g} ps"
            )
        yield analysis_ts
