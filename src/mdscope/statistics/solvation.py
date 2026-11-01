"""Replicate-level descriptive RDF outputs, without treating RDF bins as independent."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import atomic_parquet, ensure_dir


RDF_KEYS = ["solvent", "target_class", "site_class", "r_inner_A", "r_outer_A"]


def analyze_solvation_rdfs(cfg, root, system_repouts):
    """Report RDF values per replica and mean±SD across independent replicas.

    RDF is a nonlinear ratio, so pointwise averaging g(r) equally across
    replicates is intentional. This is descriptive SD, *not* a CI or
    independent-bin p-value. Solvation scalar block bootstrapping is handled
    by the pre-existing analyze_time_series workflow.
    """
    if not cfg.get("statistics", {}).get("enabled", True):
        return
    frames = []
    for system, reps in system_repouts.items():
        for repout in reps:
            source = Path(repout) / "extensions" / "solvation" / "rdf.parquet"
            if not source.exists():
                continue
            rdf = pd.read_parquet(source)
            required = set(RDF_KEYS + ["g_bulk"])
            if not required.issubset(rdf):
                raise ValueError(
                    f"incomplete solvation RDF in {source}: {sorted(required-set(rdf))}")
            if rdf.empty:
                continue
            if rdf.duplicated(RDF_KEYS).any():
                raise ValueError(f"duplicate RDF radial bins in {source}")
            frames.append(rdf[RDF_KEYS+["g_bulk"]].assign(
                system=system, replicate=Path(repout).name))
    if not frames:
        return
    reps = pd.concat(frames, ignore_index=True)
    keys = ["system", *RDF_KEYS]
    system = reps.groupby(keys, as_index=False).agg(
        n_replicates=("g_bulk", "count"),
        mean_g_bulk=("g_bulk", "mean"),
        sd_g_bulk=("g_bulk", "std"))
    # A missing replicate/bin must not be presented as a complete
    # uncertainty estimate for that chemical class.
    expected = reps.groupby(
        ["system", "solvent", "target_class", "site_class"])["replicate"].nunique()
    count = system.set_index(["system", "solvent", "target_class", "site_class"]).index
    system["complete_replica_coverage"] = (
        system.n_replicates.to_numpy() == expected.reindex(count).to_numpy())
    dest = ensure_dir(Path(root)/"project_summary"/"statistics"/"solvation")
    atomic_parquet(dest/"replicate_rdf.parquet", reps)
    atomic_parquet(dest/"system_rdf.parquet", system)
