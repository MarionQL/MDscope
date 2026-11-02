"""Continuous lipid-molecule contact episodes from sparse motif contact events.

Events are sampled frames, not physical binding free energies. Lifetimes are
time-resolution-limited, and episodes touching either edge are censored.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


COLUMNS = [
    "motif", "species", "molecule_id", "start_ps", "end_ps", "duration_ps",
    "n_frames", "left_censored", "right_censored", "occupancy_fraction",
    "episode_contact_fraction", "observed_span_ps", "duration_convention",
]


def summarize_lipid_residence(events: pd.DataFrame, frame_times, gap_frames=0):
    """Sampled contact episodes with at most gap_frames absent frames.

    Bridged gaps do not imply continuous microscopic contact. Durations
    include the gap and are distinct from occupancy, which counts only
    observed contact frames.
    """
    if not isinstance(gap_frames, int) or gap_frames < 0:
        raise ValueError("gap_frames must be a nonnegative integer")
    times = np.asarray(frame_times, dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or (
        len(times) > 1 and np.any(np.diff(times) <= 0)
    ):
        raise ValueError("Frame times must increase strictly")
    if events.empty:
        return pd.DataFrame(columns=COLUMNS)
    required = {"time_ps", "motif", "molecule_id", "species"}
    if not required.issubset(events.columns):
        raise ValueError("Missing residence event fields")
    if events.duplicated(["time_ps", "motif", "molecule_id"]).any():
        raise ValueError("Duplicate molecule contact events")
    indices = pd.Index(times).get_indexer(events["time_ps"].to_numpy(float))
    if np.any(indices < 0):
        raise ValueError("Contact events do not align with trajectory times")
    tagged = events.assign(frame_idx=indices)
    dt = float(np.median(np.diff(times))) if len(times) > 1 else np.nan
    rows = []
    for (motif, molecule), group in tagged.groupby(["motif", "molecule_id"], sort=True):
        if group.species.nunique() != 1:
            raise ValueError("Molecule changes species")
        frames = np.sort(group.frame_idx.to_numpy(int))
        cuts = np.r_[0, np.flatnonzero(np.diff(frames) > gap_frames + 1) + 1, len(frames)]
        for begin, end in zip(cuts[:-1], cuts[1:]):
            left, right = frames[begin], frames[end - 1]
            step = times[right+1] - times[right] if right+1 < len(times) else dt
            rows.append({
                "motif": motif, "species": group.species.iloc[0],
                "molecule_id": int(molecule),
                "start_ps": float(times[left]), "end_ps": float(times[right]),
                "duration_ps": float(times[right] - times[left] + step),
                "n_frames": int(end - begin),
                "left_censored": bool(left == 0),
                "right_censored": bool(right == len(times)-1),
                "occupancy_fraction": float(len(frames) / len(times)),
                "episode_contact_fraction": float((end - begin) / len(times)),
                "observed_span_ps": float(times[right] - times[left]),
                "duration_convention": "next_sample_interval_or_terminal_median",
            })
    return pd.DataFrame(rows, columns=COLUMNS)


def summarize_intermittent_residence(events: pd.DataFrame, frame_times, gap_frames=1):
    """Gap-tolerant residence episodes, reported separately from continuous ones."""
    return summarize_lipid_residence(events, frame_times, gap_frames=gap_frames)
