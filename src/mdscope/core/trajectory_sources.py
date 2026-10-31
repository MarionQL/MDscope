"""Coordinate provenance for raw MSD alongside preprocessed analysis trajectories."""
from __future__ import annotations

import numpy as np
import MDAnalysis as mda

from ..observables.msd import initialize_msd


def open_msd_universe(rep, cfg, membrane, main_universe, first_frame):
    """Open a second trajectory only when explicitly requested for active MSD.

    Both readers must use the identical topology and frame indexing. The
    separate MSD reader never receives protein-fit transforms.
    """
    source = getattr(rep, "msd_trajectory", None)
    if source is None or not (membrane.enabled and membrane.msd):
        return None
    raw = mda.Universe(str(rep.topology), str(source))
    if len(raw.atoms) != len(main_universe.atoms):
        raise ValueError(
            "msd_trajectory and trajectory have different atom counts; "
            "both sources must use the identical topology and atom ordering"
        )
    if len(raw.trajectory) != len(main_universe.trajectory):
        raise ValueError(
            "msd_trajectory and trajectory have different frame counts; "
            "provide frame-synchronized trajectories"
        )
    raw.trajectory[int(first_frame)]
    initialize_msd(raw, membrane, cfg["analyses"]["membrane"])
    return raw


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
