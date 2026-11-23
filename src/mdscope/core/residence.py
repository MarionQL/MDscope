"""Sampling conventions shared by lipid, water, and ion contact episodes."""
import numpy as np
import pandas as pd


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


EPISODE_FIELDS = [ "start_ps", "end_ps", "duration_ps",
    "n_frames", "left_censored", "right_censored", "occupancy_fraction",
    "episode_contact_fraction", "observed_span_ps", "duration_convention",
]


def summarize_contact_residence(events, frame_times, *, identity_columns,
                                species_column, molecule_column, gap_frames=0):
    """Sampled contact episodes with at most gap_frames absent frames.

    Bridged gaps do not imply continuous microscopic contact. Durations
    include the gap and are distinct from occupancy, which counts only
    observed contact frames.
    """
    columns = [*identity_columns, species_column, molecule_column, *EPISODE_FIELDS]
    if not isinstance(gap_frames, int) or gap_frames < 0:
        raise ValueError("gap_frames must be a nonnegative integer")
    times = np.asarray(frame_times, dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all() or (
        len(times) > 1 and np.any(np.diff(times) <= 0)
    ):
        raise ValueError("Frame times must increase strictly")
    if events.empty:
        return pd.DataFrame(columns=columns)
    required = {"time_ps", *identity_columns, molecule_column, species_column}
    if not required.issubset(events.columns):
        raise ValueError("Missing residence event fields")
    if events.duplicated(["time_ps", *identity_columns, molecule_column]).any():
        raise ValueError("Duplicate molecule contact events")
    indices = pd.Index(times).get_indexer(events["time_ps"].to_numpy(float))
    if np.any(indices < 0):
        raise ValueError("Contact events do not align with trajectory times")
    tagged = events.assign(frame_idx=indices)
    rows, bounds = [], []
    for key, group in tagged.groupby([*identity_columns, molecule_column], sort=True):
        if group[species_column].nunique() != 1:
            raise ValueError("Molecule changes species")
        frames = np.sort(group.frame_idx.to_numpy(int))
        cuts = np.r_[0, np.flatnonzero(np.diff(frames) > gap_frames + 1) + 1, len(frames)]
        for begin, end in zip(cuts[:-1], cuts[1:]):
            left, right = frames[begin], frames[end - 1]
            bounds.append((left, right))
            rows.append({
                **dict(zip(identity_columns, key[:-1])),
                species_column: group[species_column].iloc[0],
                molecule_column: int(key[-1]),
                "start_ps": float(times[left]), "end_ps": float(times[right]),
                "duration_ps": np.nan,
                "n_frames": int(end - begin),
                "left_censored": bool(left == 0),
                "right_censored": bool(right == len(times)-1),
                "occupancy_fraction": float(len(frames) / len(times)),
                "episode_contact_fraction": float((end - begin) / len(times)),
                "observed_span_ps": np.nan,
                "duration_convention": "next_sample_interval_or_terminal_median",
            })
    result = pd.DataFrame(rows, columns=columns)
    first, last = np.asarray(bounds, dtype=int).T
    result["observed_span_ps"], result["duration_ps"] = sampled_episode_durations(
        times, first, last
    )
    return result
