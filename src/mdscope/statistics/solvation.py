"""Replicate RDF curves and predefined first-shell contrasts."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import statistics_family_config, atomic_parquet, ensure_dir
from .replicate_contrasts import replicate_mean_contrasts


RDF_KEYS = ["solvent", "target_class", "site_class", "r_inner_A", "r_outer_A"]


def analyze_solvation_rdfs(cfg, root, system_repouts):
    """Report RDF values per replica and mean±SD across independent replicas.

    RDF is a nonlinear ratio, so pointwise averaging g(r) equally across
    replicates is intentional. This is descriptive SD, *not* a CI or
    independent-bin p-value. Solvation scalar block bootstrapping is handled
    by the pre-existing analyze_time_series workflow.
    """
    cfg = statistics_family_config(cfg, "solvation_rdfs")
    if cfg is None:
        return
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
            frames.append(rdf[RDF_KEYS+["g_bulk"] +
                (["g_aqueous_reference"] if "g_aqueous_reference" in rdf else [])].assign(
                system=system, replicate=Path(repout).name))
    if not frames:
        return
    reps = pd.concat(frames, ignore_index=True)
    keys = ["system", *RDF_KEYS]
    system = reps.groupby(keys, as_index=False).agg(
        n_replicates=("g_bulk", "count"),
        mean_g_bulk=("g_bulk", "mean"),
        sd_g_bulk=("g_bulk", "std"))
    if "g_aqueous_reference" in reps:
        reference = reps.groupby(keys, as_index=False).agg(
            n_aqueous_replicates=("g_aqueous_reference", "count"),
            mean_g_aqueous_reference=("g_aqueous_reference", "mean"),
            sd_g_aqueous_reference=("g_aqueous_reference", "std"))
        system = system.merge(reference, on=keys, validate="one_to_one")
    # A missing replicate/bin must not be presented as a complete
    # uncertainty estimate for that chemical class.
    # Compare against requested replicas, including those whose entire RDF
    # file/class is absent. Inferring the denominator from present data can
    # incorrectly label an incomplete analysis as complete.
    expected = {name: len(paths) for name, paths in system_repouts.items()}
    system["expected_replicates"] = system.system.map(expected).astype(int)
    system["complete_replica_coverage"] = (
        system.n_replicates == system.expected_replicates)
    if "n_aqueous_replicates" in system:
        system["complete_aqueous_replica_coverage"] = (
            system.n_aqueous_replicates == system.expected_replicates)
    dest = ensure_dir(Path(root)/"project_summary"/"statistics"/"solvation")
    atomic_parquet(dest/"replicate_rdf.parquet", reps)
    atomic_parquet(dest/"system_rdf.parquet", system)
    _write_shell_contrasts(reps, dest, cfg, system_repouts)


def _write_shell_contrasts(reps, dest, cfg, system_repouts):
    """Compare shell-volume-weighted mean g(r) at the configured cutoffs."""
    options = cfg.get("analyses", {}).get("solvation", {})
    cutoffs = {"water": 3.5, "ions": 4.0} | (options.get("shell_cutoffs_A") or {})
    keys = ["solvent", "target_class", "site_class"]
    rows = []
    for identity, rdf in reps.groupby(["system", "replicate", *keys], sort=True):
        cutoff = cutoffs.get(identity[2].split("/", 1)[0])
        if cutoff is None:
            continue
        cutoff = float(cutoff)
        shell = rdf.loc[rdf.r_inner_A.lt(cutoff)].sort_values("r_inner_A")
        left = shell.r_inner_A.to_numpy(float)
        right = shell.r_outer_A.to_numpy(float)
        complete = (len(left) > 0 and np.isclose(left[0], 0)
                    and (right[-1] >= cutoff or np.isclose(right[-1], cutoff))
                    and np.all(right > left) and np.allclose(left[1:], right[:-1]))
        for normalization in ("g_bulk", "g_aqueous_reference"):
            if normalization not in reps:
                continue
            density = shell[normalization].to_numpy(float)
            value = np.nan
            if complete and np.isfinite(density).all():
                weights = np.minimum(right, cutoff) ** 3 - left ** 3
                value = float(np.average(density, weights=weights))
            rows.append(dict(zip(["system", "replicate", *keys], identity),
                             normalization=normalization, shell_cutoff_A=cutoff, value=value))
    if not rows:
        return
    shell_values = pd.DataFrame(rows)
    comparisons = replicate_mean_contrasts(
        shell_values, [*keys, "normalization", "shell_cutoff_A"], system_repouts, cfg["statistics"]
    )
    atomic_parquet(dest / "replicate_rdf_shells.parquet", shell_values)
    if not comparisons.empty:
        atomic_parquet(dest / "rdf_shell_comparisons.parquet", comparisons)
