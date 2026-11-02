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
]


def summarize_lipid_residence(events: pd.DataFrame, frame_times) -> pd.DataFrame:
    times = np.asarray(frame_times, dtype=float)
    if times.ndim != 1 or len(times) == 0 or not np.isfinite(times).all():
        raise ValueError("Expected finite, nonempty time grid")
    if len(np.unique(times)) != len(times) or np.any(np.diff(times) <= 0):
        raise ValueError("Frame times must increase strictly")
    if events.empty:
        return pd.DataFrame(columns=COLUMNS)
    required = {"time_ps", "motif", "molecule_id", "species"}
    if not required.issubset(events.columns):
        raise ValueError(f"Missing contact event columns: {sorted(required - set(events.columns))}")
    if events.duplicated(subset=["time_ps", "motif", "molecule_id"]).any():
        raise ValueError("Duplicate contact events for the same molecule and frame")
    lookup = pd.Index(times)
    indices = lookup.get_indexer(events["time_ps"].to_numpy(float))
    if np.any(indices < 0):
        raise ValueError("Contact events do not align with observed trajectory times")
    tagged = events.copy()
    tagged["frame_idx"] = indices
    dt = float(np.median(np.diff(times))) if len(times) > 1 else np.nan
    rows = []
    for (motif, molecule), group in tagged.groupby(["motif", "molecule_id"], sort=True):
        if group["species"].nunique() != 1:
            raise ValueError("One lipid molecule cannot change species")
        frames = np.sort(group["frame_idx"].to_numpy(int))
        boundaries = np.r_[0, np.flatnonzero(np.diff(frames) > 1) + 1, len(frames)]
        for begin, end in zip(boundaries[:-1], boundaries[1:]):
            start, stop = int(frames[begin]), int(frames[end-1])
            # The last sampled contact represents at most one sampling interval.
            # For unequal spacing, include the next interval if available.
            step = float(times[stop+1] - times[stop]) if stop+1 < len(times) else dt
            rows.append({
                "motif": motif, "species": group["species"].iloc[0],
                "molecule_id": int(molecule),
                "start_ps": float(times[start]), "end_ps": float(times[stop]),
                "duration_ps": float(times[stop] - times[start] + step),
                "n_frames": int(end - begin),
                "left_censored": bool(start == 0),
                "right_censored": bool(stop == len(times) - 1),
                "occupancy_fraction": float(len(frames) / len(times)),
            })
    return pd.DataFrame(rows, columns=COLUMNS)


def summarize_intermittent_residence(events: pd.DataFrame, frame_times, gap_frames=1):
    """Episodes allowing at most gap_frames *absent* sampled frames.

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
                "occupancy_fraction": float((end - begin) / len(times)),
            })
    return pd.DataFrame(rows, columns=COLUMNS)
