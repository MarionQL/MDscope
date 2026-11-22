"""Load explicitly identified saved scalar features onto one exact time grid."""
import numpy as np
import pandas as pd

from .time_series import collect_time_series


def load_state_features(repout, features, sources):
    frame = None
    for name in features:
        spec = sources[name]
        identity = (spec.get("group"), spec["metric"])
        available = collect_time_series(repout, sources={spec["source"]}, selected={identity})
        if identity not in available:
            raise ValueError(f"Shared-state feature {name!r} unavailable in {repout}: {spec}")
        series = available[identity].sort_values("time_ps").reset_index(drop=True)
        times = series.time_ps.to_numpy(float)
        if not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError(f"Shared-state feature {name!r} has invalid or duplicate time grid in {repout}")
        if frame is None:
            frame = pd.DataFrame({"time_ps": times})
        elif not np.array_equal(times, frame.time_ps.to_numpy()):
            raise ValueError(f"Shared-state feature {name!r} time grid differs in {repout}; exact alignment required")
        frame[name] = series.value.to_numpy(float)
    return frame
