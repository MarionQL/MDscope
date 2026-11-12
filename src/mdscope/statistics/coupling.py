"""Matched-window coupling estimates and between-condition intervals."""
from pathlib import Path

import numpy as np
import pandas as pd
from pyarrow.parquet import read_schema

from ..config import atomic_parquet, statistics_family_config
from .common import _regular_times
from .replicate_contrasts import replicate_mean_contrasts


PAIRS=[
 ("lipid_proximity","local_hs_lipid_contact_pairs_count","local_protein_lipid_contact_residues_count"),
 ("ion_screening","local_hs_ion_contacts_count","local_hs_lipid_contact_pairs_count"),
 ("lipid_contact_rigidity","local_protein_lipid_contact_residues_count","radius_of_gyration_A"),
 ("membrane_coupling","local_protein_lipid_contact_residues_count","membrane_thickness_A")
]

def block_means(x,block):
    n=len(x)//block
    return x[:n*block].reshape(n,block).mean(axis=1) if n else np.empty(0)

def correlation(x,y,min_blocks=8):
    a=np.asarray(x,float);b=np.asarray(y,float)
    good=np.isfinite(a)&np.isfinite(b)
    a,b=a[good],b[good]
    if len(a)<min_blocks or np.std(a)==0 or np.std(b)==0:return np.nan
    return float(np.corrcoef(a,b)[0,1])

def collect_couplings(cfg, system_repouts):
    spec = cfg.get("analyses", {}).get("coupling_screen", {})
    pairs = list(spec.get("pairs", []))
    if spec.get("include_defaults", True):
        pairs = [dict(name=name, x=x, y=y) for name, x, y in PAIRS] + pairs
    pairs = [(str(pair["name"]), str(pair["x"]), str(pair["y"])) for pair in pairs]
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate coupling pair")
    block = int(spec.get("block_frames", 50))
    if block < 1:
        raise ValueError("block_frames must be positive")

    requested = {"time_ps", *(column for _, x, y in pairs for column in (x, y))}
    sources = {}
    for system, reps in system_repouts.items():
        for rep in reps:
            path = Path(rep) / "summary" / "timeseries.parquet"
            if path.is_file():
                columns = [name for name in read_schema(path).names if name in requested]
                frame = pd.read_parquet(path, columns=columns)
                if "time_ps" in frame and _regular_times(frame.time_ps):
                    sources[system, rep.name] = frame.set_index("time_ps")
    common = None
    for frame in sources.values():
        common = frame.index if common is None else common.intersection(frame.index)
    shared_grid = common is not None and _regular_times(common)
    if shared_grid:
        step = float(np.diff(common)[0])
        shared_grid = all(np.isclose(np.diff(frame.index)[0], step)
                          for frame in sources.values())
    rows = []
    for system, reps in system_repouts.items():
        for rep in reps:
            frame = sources.get((system, rep.name))
            means = {
                column: block_means(values.to_numpy(float), block)
                for column, values in frame.loc[common].items()
            } if shared_grid and frame is not None else {}
            for name, x, y in pairs:
                r, n = np.nan, 0
                if x in means and y in means:
                    a, b = means[x], means[y]
                    n = int((np.isfinite(a) & np.isfinite(b)).sum())
                    r = correlation(a, b)
                rows.append(dict(system=system, replicate=rep.name, mechanism=name, x=x, y=y,
                                 pearson_block_r=r, block_frames=block, n_blocks=n,
                                 window_start_ps=float(common[0]) if shared_grid else np.nan,
                                 window_end_ps=float(common[-1]) if shared_grid else np.nan,
                                 interpretation="exploratory_association_not_causation"))
    return pd.DataFrame(rows)


def analyze_coupling(cfg, root, system_repouts):
    cfg = statistics_family_config(cfg, "coupling")
    if cfg is None or not cfg.get("analyses", {}).get("coupling_screen", {}).get("enabled", False):
        return
    frame = collect_couplings(cfg, system_repouts)
    if frame.empty:
        return
    dest = Path(root) / "project_summary" / "statistics" / "coupling"
    values = frame.rename(columns={"pearson_block_r": "value"})
    comparisons = replicate_mean_contrasts(
        values, ["mechanism", "x", "y"], system_repouts, cfg["statistics"]
    )
    atomic_parquet(dest / "replicate_correlations.parquet", frame)
    if not comparisons.empty:
        atomic_parquet(dest / "system_comparisons.parquet", comparisons)
