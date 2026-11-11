"""Shared sampled timestamps for optional matched-window comparisons."""
import numpy as np
import pandas as pd

from .common import _regular_times


def matched_grid(cfg, system_repouts):
    requested = cfg.get("statistics", {}).get("matched_window_ps")
    if requested is None:
        return None
    grids = []
    for reps in system_repouts.values():
        for rep in reps:
            path = rep / "summary" / "timeseries.parquet"
            grid = pd.read_parquet(path, columns=["time_ps"]).time_ps.to_numpy(float)
            if not _regular_times(grid):
                raise ValueError(f"Matched windows require a regular frame grid in {path}")
            grids.append(grid)
    if not grids:
        raise ValueError("Matched windows require at least one replica")
    lower = max(grid[0] for grid in grids)
    upper = min(grid[-1] for grid in grids)
    if requested != "auto":
        if requested[0] < lower or requested[1] > upper:
            raise ValueError("Matched window extends outside a replica's saved frame range")
        lower, upper = requested
    trimmed = [grid[(grid >= lower) & (grid <= upper)] for grid in grids]
    if not _regular_times(trimmed[0]) or any(
            not np.array_equal(trimmed[0], grid) for grid in trimmed[1:]):
        raise ValueError("Matched windows require identical sampled timestamps across replicas")
    return trimmed[0]
