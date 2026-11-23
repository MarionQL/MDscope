"""Exploratory raw, per-replica correlations between matched residue time series.

No frame pooling across replicas, bootstrap, p-values or condition-level inference.
Time and topology-residue joins are exact; unavailable variables remain missing.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ..core.persistence import atomic_parquet

DEFAULT_PAIRS = (
    ("water_count", "protein_fraction"),
    ("water_count", "hs_fraction"),
    ("water_count", "lipid_fraction"),
    ("water_count", "sasa_isolated"),
    ("water_count", "sasa_environment"),
    ("water_count", "helix_state"),
    ("water_count", "interprotomer_fraction"),
    ("interprotomer_fraction", "helix_state"),
    ("interprotomer_fraction", "sasa_isolated"),
    ("helix_state", "protein_fraction"),
    ("helix_state", "hs_fraction"),
    ("helix_state", "sasa_isolated"),
    ("protein_fraction", "sasa_isolated"),
    ("protein_fraction", "sasa_environment"),
)
METRIC_SOURCES = {
    "water_count": ("extensions/sdc4_mechanisms/residue_occlusion.parquet", "water_count"),
    "protein_fraction": ("extensions/sdc4_mechanisms/residue_occlusion.parquet", "protein_fraction"),
    "hs_fraction": ("extensions/sdc4_mechanisms/residue_occlusion.parquet", "hs_fraction"),
    "lipid_fraction": ("extensions/sdc4_mechanisms/residue_occlusion.parquet", "lipid_fraction"),
    "sasa_isolated": ("profiles/sasa/isolated/residue_sasa.parquet", "sasa_nm2"),
    "sasa_environment": ("profiles/sasa/environment/residue_sasa.parquet", "sasa_nm2"),
    "helix_state": ("profiles/dssp/assignments.parquet", "dssp"),
    "interprotomer_fraction": (
        "extensions/sdc4_mechanisms/interprotomer_residue_frames.parquet",
        "contacting_heavy_atom_fraction",
    ),
}


def _series(repout: Path, metric: str, *, sasa_runs=None):
    relative, field = METRIC_SOURCES[metric]
    if metric.startswith("sasa_") and sasa_runs:
        run = sasa_runs.get(metric)
        if run is None:
            return None, {}
        relative = ("profiles/sasa/" + (run + "/" if run != "default" else "")
                    + "residue_sasa.parquet")
    path = repout / relative
    if not path.is_file():
        return None, {}
    frame = pd.read_parquet(path)
    needed = {"time_ps", "topology_residue_index", "residue", field}
    if not needed.issubset(frame):
        return None, {}
    if frame.duplicated(["time_ps", "topology_residue_index"]).any():
        raise ValueError(f"{path}: duplicate time/residue identity")
    labels = frame[["topology_residue_index", "residue"]].drop_duplicates()
    if labels.duplicated("topology_residue_index").any():
        raise ValueError(f"{path}: changing residue labels within a topology index")
    names = dict(zip(labels.topology_residue_index.astype(int), labels.residue.astype(str)))
    if metric == "helix_state":
        values = frame[field].map({"H": 1.0, "E": 0.0, "C": 0.0})
    else:
        values = pd.to_numeric(frame[field], errors="coerce")
    series = pd.Series(values.to_numpy(float), index=pd.MultiIndex.from_frame(
        frame[["time_ps", "topology_residue_index"]]
    ), name=metric)
    return series, names


COLUMNS = (
    "system", "replicate", "topology_residue_index", "residue",
    "metric_x", "metric_y", "n_frames_x", "n_frames_y", "n_paired_frames",
    "pearson_r", "spearman_rho", "status", "first_time_ps", "last_time_ps",
    "interpretation",
)


def correlate_replica(repout: Path, pairs=DEFAULT_PAIRS, *, min_frames=30,
                      sasa_runs=None):
    """Return one descriptive row per residue and configured metric pair."""
    metrics = sorted({m for pair in pairs for m in pair})
    sources = {}
    labels = {}
    for metric in metrics:
        series, mapping = _series(Path(repout), metric, sasa_runs=sasa_runs)
        if series is None:
            continue
        for index, name in mapping.items():
            if index in labels and labels[index] != name:
                raise ValueError(f"{repout}: residue identity differs between sources at index {index}")
            labels[index] = name
        sources[metric] = series
    if not sources:
        return pd.DataFrame(columns=COLUMNS)
    panel = pd.concat(sources.values(), axis=1, join="outer").sort_index()
    rows = []
    for residx, group in panel.groupby(level="topology_residue_index", sort=True):
        times = group.index.get_level_values("time_ps").to_numpy(float)
        for x, y in pairs:
            if x not in sources or y not in sources:
                rows.append(_row(int(residx), labels[int(residx)], x, y,
                                 0, 0, 0, np.nan, np.nan, "missing_metric",
                                 np.nan, np.nan))
                continue
            a = group[x].to_numpy(float)
            b = group[y].to_numpy(float)
            nx, ny = int(np.isfinite(a).sum()), int(np.isfinite(b).sum())
            mask = np.isfinite(a) & np.isfinite(b)
            n = int(mask.sum())
            status = "ok"
            pearson = spearman = np.nan
            if n < min_frames:
                status = "insufficient_paired_frames"
            elif np.ptp(a[mask]) == 0 or np.ptp(b[mask]) == 0:
                status = "constant_series"
            else:
                pearson = float(np.corrcoef(a[mask], b[mask])[0, 1])
                ra = pd.Series(a[mask]).rank(method="average").to_numpy(float)
                rb = pd.Series(b[mask]).rank(method="average").to_numpy(float)
                spearman = float(np.corrcoef(ra, rb)[0, 1])
            rows.append(_row(int(residx), labels[int(residx)], x, y, nx, ny, n,
                             pearson, spearman, status,
                             float(times[mask].min()) if n else np.nan,
                             float(times[mask].max()) if n else np.nan))
    return pd.DataFrame(rows, columns=COLUMNS)


def _row(index, label, x, y, nx, ny, n, pearson, spearman, status, first, last):
    return {
        "topology_residue_index": index, "residue": label,
        "metric_x": x, "metric_y": y,
        "n_frames_x": nx, "n_frames_y": ny, "n_paired_frames": n,
        "pearson_r": pearson, "spearman_rho": spearman, "status": status,
        "first_time_ps": first, "last_time_ps": last,
        "interpretation": "exploratory_raw_within_replica_not_independent_frames",
    }


def write_raw_residue_correlations(cfg, root, system_repouts):
    """Consume saved replica products; never reopen an XTC or infer uncertainty."""
    from ..config.validators.statistics import residue_correlations_settings
    settings = residue_correlations_settings(cfg)
    pairs = settings["pairs"] if settings["pairs"] is not None else DEFAULT_PAIRS
    mechanism = cfg.get("analyses", {}).get("sdc4_mechanisms", {})
    sasa_runs = {
        "sasa_isolated": mechanism.get("sasa_run", "default"),
        "sasa_environment": mechanism.get("environment_sasa_run"),
    } if mechanism.get("enabled") else None
    rows = []
    for system, reps in system_repouts.items():
        for repout in reps:
            table = correlate_replica(
                Path(repout), pairs, min_frames=settings["min_frames"],
                sasa_runs=sasa_runs,
            )
            if not table.empty:
                table.insert(0, "replicate", Path(repout).name)
                table.insert(0, "system", system)
                rows.append(table)
    outdir = Path(root) / "project_summary/statistics/residue_correlations"
    outdir.mkdir(parents=True, exist_ok=True)
    atomic_parquet(outdir / "replicate_correlations.parquet",
                   pd.concat(rows, ignore_index=True) if rows
                   else pd.DataFrame(columns=COLUMNS))
    return sum(len(table) for table in rows)
