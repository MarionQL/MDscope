"""Opt-in bounded-RAM adapter for multi-origin lipid diffusion.

This intentionally runs in a separate pass, leaving streaming/checkpoint
state untouched. Positions are staged in a temporary on-disk memmap; FFT
processing only loads a limited number of lipid tracks at a time.
"""
from __future__ import annotations
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors

from ..config import atomic_parquet, time_ps
from .membrane import initialize_membrane
from .msd import _residue_centers


def run_multi_origin_diffusion(rep, cfg, repout):
    settings = cfg["analyses"]["membrane"]["msd"]
    if not settings.get("multi_origin", False):
        return {}
    source = rep.msd_trajectory or rep.trajectory
    u = mda.Universe(str(rep.topology), str(source))
    tr = cfg["trajectory"]
    frames = u.trajectory[slice(tr.get("start"), tr.get("stop"), int(tr["stride"]))]
    n_frames = len(frames)
    if n_frames < 6:
        raise ValueError("Multi-origin MSD requires at least six sampled frames")
    u.trajectory[0 if tr.get("start") is None else int(tr["start"])]
    membrane = initialize_membrane(u, cfg)
    if not membrane.msd:
        raise ValueError("Multi-origin MSD requires analyses.membrane.msd.enabled")
    n_lipids = len(membrane.msd_resindices)
    axis = membrane.ai
    max_lag = settings.get("multi_origin_max_lag_frames")
    max_lag = None if max_lag is None else int(max_lag)
    batch = int(settings.get("multi_origin_particle_batch", 64))
    if batch < 1 or (max_lag is not None and not 4 <= max_lag < n_frames):
        raise ValueError("Invalid multi-origin particle batch or maximum lag")
    output = Path(repout) / "summary"
    output.mkdir(parents=True, exist_ok=True)
    times = np.empty(n_frames, dtype=np.float64)
    with TemporaryDirectory(prefix="mdscope-multi-origin-", dir=Path(repout)) as tmp:
        track = np.lib.format.open_memmap(
            Path(tmp) / "unwrapped.npy", mode="w+", dtype=np.float64,
            shape=(n_frames, n_lipids, 3),
        )
        last = None
        displacement = np.zeros((n_lipids, 3), dtype=np.float64)
        for i, ts in enumerate(frames):
            times[i] = time_ps(ts, i, cfg)
            current = _residue_centers(membrane.msd_atoms)
            if last is not None:
                step = np.asarray(minimize_vectors(current - last, u.dimensions), float)
                if membrane.msd_remove_drift:
                    step -= step.mean(axis=0, keepdims=True)
                displacement += step
            track[i] = displacement
            last = current
        if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0) or not np.allclose(np.diff(times), np.diff(times)[0], rtol=1e-5, atol=1e-8):
            raise ValueError("Multi-origin MSD requires regularly spaced finite time samples")
        track.flush()
        rows, curves = [], []
        from .multi_origin_msd import multi_origin_msd, fit_multi_origin_curve
        for suffix, mask in (
            ("", np.ones(n_lipids, dtype=bool)),
            ("_upper", membrane.msd_upper_mask),
            ("_lower", membrane.msd_lower_mask),
        ):
            if not np.any(mask):
                continue
            # Boolean fancy indexing copies the full selection. Preserve
            # bounded RAM by processing contiguous tracks in small batches.
            selected = np.flatnonzero(mask)
            lag = None
            weighted = None
            for start in range(0, len(selected), batch):
                indices = selected[start:start + batch]
                lags, curve = multi_origin_msd(
                    track[:, indices, :], normal_axis=axis, max_lag=max_lag,
                    particle_batch=batch,
                )
                if lag is None:
                    lag = lags
                    weighted = curve * len(indices)
                else:
                    weighted += curve * len(indices)
            # Fit against the independent lag MSD using the shared estimator.
            curve = weighted / len(selected)
            lag_times = lag * (times[1] - times[0])
            fit = fit_multi_origin_curve(
                lag_times, curve,
                fit_start_fraction=float(settings.get("slope_start_fraction", .2)),
                fit_end_fraction=float(settings.get("slope_end_fraction", .8)),
            )
            metric = "membrane_msd_lateral" + suffix + "_A2"
            rows.append({
                "system": rep.system, "replicate": rep.label, "metric": metric,
                "method": "multi_origin_fft", "n_lipids": int(len(selected)),
                "n_frames": int(n_frames), **fit,
            })
            curves.extend({
                "metric": metric, "lag_ps": float(t), "msd_A2": float(v),
                "n_time_origins": int(n_frames - j),
            } for j, (t, v) in enumerate(zip(lag_times, curve)))
        atomic_parquet(output / "multi_origin_diffusion.parquet", pd.DataFrame(rows))
        atomic_parquet(output / "multi_origin_msd.parquet", pd.DataFrame(curves))
    return {"frames": n_frames, "lipids": n_lipids, "method": "multi_origin_fft"}
