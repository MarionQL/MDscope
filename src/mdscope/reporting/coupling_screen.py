"""Report the replicate-level coupling screen using the statistics estimator."""
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import ensure_dir, atomic_parquet
from ..statistics.coupling import PAIRS, block_means, collect_couplings, correlation

__all__ = ["PAIRS", "block_means", "correlation", "screen_couplings"]


def screen_couplings(cfg, root, system_repouts):
    frame = collect_couplings(cfg, system_repouts)
    if frame.empty:
        return
    out=ensure_dir(Path(root)/"project_summary"/"mechanisms")
    atomic_parquet(out/"replicate_correlations.parquet",frame)
    summaries=[]
    for keys,d in frame.groupby(["system","mechanism","x","y"],dropna=False):
        rs=d["pearson_block_r"].to_numpy(float);rs=rs[np.isfinite(rs)]
        if not len(rs):continue
        summaries.append(dict(zip(("system","mechanism","x","y"),keys),
            n_replicates=len(rs),median_r=float(np.median(rs)),
            min_r=float(rs.min()),max_r=float(rs.max()),
            fraction_same_sign=float(max(np.mean(rs>0),np.mean(rs<0))),
            note="associational; no causal or statistical significance claim"))
    if summaries:atomic_parquet(out/"mechanism_consistency.parquet",pd.DataFrame(summaries))
    return out
