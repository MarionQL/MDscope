"""Between-replica intervals for a shared-lag multi-origin MSD slope."""
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import atomic_parquet, statistics_family_config
from ..observables.multi_origin_msd import fit_multi_origin_curve
from .replicate_contrasts import replicate_mean_contrasts


def analyze_multi_origin_diffusion(cfg, root, system_repouts):
    cfg = statistics_family_config(cfg, "diffusion")
    if cfg is None:
        return
    settings = cfg.get("analyses", {}).get("membrane", {}).get("msd", {})
    if not settings.get("enabled") or not settings.get("multi_origin"):
        return
    curves = {}
    requested = settings.get("multi_origin_window_ps")
    time_columns = ["time_window_start_ps", "time_window_end_ps"]
    for system, reps in system_repouts.items():
        for rep in reps:
            path = Path(rep) / "summary" / "multi_origin_msd.parquet"
            if path.is_file():
                data = pd.read_parquet(path)
                if requested is not None and not set(time_columns) <= set(data):
                    raise ValueError(f"{path} lacks trajectory-window provenance. Regenerate its multi-origin MSD.")
                for metric, curve in data.groupby("metric", sort=True):
                    curve = curve.set_index("lag_ps").sort_index()
                    lag = curve.index.to_numpy(float)
                    origins = curve.n_time_origins.to_numpy(float)
                    msd = curve.msd_A2.to_numpy(float)
                    if (not np.isfinite(lag).all() or lag[0] != 0 or not curve.index.is_unique
                            or not np.isfinite(origins).all() or np.any(origins < 1)
                            or not np.isfinite(msd).all() or np.any(msd < 0)):
                        raise ValueError(f"Invalid multi-origin curve for {metric} in {path}")
                    if requested is not None:
                        observed = curve[time_columns].drop_duplicates().to_numpy(float)
                        step = float(np.diff(lag)[0]) if len(lag) > 1 else 0.0
                        if (observed.shape != (1, 2) or step <= 0
                                or not np.isfinite(observed).all()
                                or np.any(np.abs(observed[0] - requested) >= step)):
                            raise ValueError(f"{path} does not match the requested trajectory window. Regenerate its multi-origin MSD.")
                    curves.setdefault(metric, {})[system, rep.name] = curve
    start = float(settings.get("slope_start_fraction", .2))
    end = float(settings.get("slope_end_fraction", .8))
    if not 0 <= start < end <= 1:
        raise ValueError("Diffusion slope fractions require 0 <= start < end <= 1")
    rows = []
    for metric, by_replica in sorted(curves.items()):
        if requested is not None:
            observed = {tuple(curve[time_columns].iloc[0]) for curve in by_replica.values()}
            if len(observed) != 1:
                raise ValueError(f"{metric}: replica MSD curves have different trajectory windows")
        duration = min(float(curve.index[-1]) for curve in by_replica.values())
        shared = None
        for curve in by_replica.values():
            sufficient_origins = curve.n_time_origins.ge(max(10, np.ceil(.1 * curve.n_time_origins.iloc[0])))
            lag = curve.index.to_numpy(float)
            supported = lag[(lag > 0) & (lag >= start * duration) & (lag <= end * duration) & sufficient_origins]
            shared = supported if shared is None else np.intersect1d(shared, supported)
        for (system, replicate), curve in sorted(by_replica.items()):
            row = dict(metric=metric, system=system, replicate=replicate, value=np.nan,
                       fit_start_ps=np.nan, fit_end_ps=np.nan, n_fit=len(shared),
                       n_frames=int(curve.n_time_origins.iloc[0]), log_log_slope=np.nan,
                       fit_status="insufficient_common_lags")
            for column in time_columns:
                row[column] = float(curve[column].iloc[0]) if column in curve else np.nan
            if len(shared) >= 5:
                lag = np.r_[0.0, shared]
                msd = curve.loc[lag, "msd_A2"].to_numpy(float)
                fit = fit_multi_origin_curve(lag, msd, fit_start_fraction=0, fit_end_fraction=1)
                row.update(fit_start_ps=shared[0], fit_end_ps=shared[-1],
                           fit_status="nonpositive_fit")
                if fit["diffusion_nm2_ns"] > 0 and np.all(msd[1:] > 0):
                    scaling = float(np.polyfit(np.log(shared), np.log(msd[1:]), 1)[0])
                    row.update(value=fit["diffusion_nm2_ns"], log_log_slope=scaling,
                               fit_status="fit" if abs(scaling - 1) <= .2 else "scaling_sensitive")
            rows.append(row)
    if not rows:
        return
    replicas = pd.DataFrame(rows)
    comparisons = replicate_mean_contrasts(replicas, ["metric"], system_repouts, cfg["statistics"])
    dest = Path(root) / "project_summary" / "statistics" / "diffusion"
    atomic_parquet(dest / "replicate_common_lag_slopes.parquet", replicas)
    if not comparisons.empty:
        atomic_parquet(dest / "common_lag_comparisons.parquet", comparisons)
