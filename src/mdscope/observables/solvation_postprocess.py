"""Offline hydration statistics: no new trajectory or SASA calculations."""
import numpy as np
import pandas as pd
from ..config import atomic_parquet


def suggest_shell_boundaries(rdf, *, min_peak_r_A=1.8,
                             min_prominence=0.15, smoothing_bins=3):
    """Suggest, but never automatically impose, a prominent first RDF minimum."""
    cols = ["solvent", "target_class", "site_class", "suggested_cutoff_A",
            "peak_r_A", "minimum_g_bulk", "status"]
    if smoothing_bins < 1 or smoothing_bins % 2 != 1:
        raise ValueError("smoothing_bins must be positive and odd")
    rows = []
    for (solvent, kind, label), frame in rdf.groupby(
            ["solvent", "target_class", "site_class"], sort=True):
        frame = frame.sort_values("r_inner_A")
        r = ((frame.r_inner_A + frame.r_outer_A) / 2).to_numpy(float)
        g = frame.g_bulk.to_numpy(float)
        cut = peak = minimum = np.nan
        status = "insufficient_or_invalid"
        if len(r) >= 7 and np.all(np.isfinite(g)):
            k = smoothing_bins // 2
            smoothed = np.convolve(np.pad(g, (k, k), mode="edge"),
                                   np.ones(smoothing_bins)/smoothing_bins,
                                   mode="valid")
            peaks = [i for i in range(1, len(r)-2) if r[i] >= min_peak_r_A
                     and smoothed[i-1] < smoothed[i] >= smoothed[i+1]]
            if not peaks:
                status = "no_first_peak"
            else:
                p = peaks[0]
                peak = r[p]
                minima = [i for i in range(p+1, len(r)-1)
                          if smoothed[i-1] > smoothed[i] <= smoothed[i+1]]
                if not minima:
                    status = "no_following_minimum"
                else:
                    m = minima[0]
                    minimum = smoothed[m]
                    if smoothed[p] - smoothed[m] < min_prominence:
                        status = "weak_peak_to_minimum_contrast"
                    else:
                        status, cut = "candidate_only", float(r[m])
        rows.append((solvent, kind, label, cut, peak, minimum, status))
    return pd.DataFrame(rows, columns=cols)


def kaplan_meier(episodes, *, group_columns=("solvent", "target_id")):
    """Censor-aware sampled-contact survival; exclude unknown left origins.

    Events from the same trajectory/solvent are correlated: do not treat
    individual episodes as independent condition-level replicates.
    """
    cols = [*group_columns, "duration_ps", "at_risk", "dissociations",
            "right_censored", "survival_probability"]
    rows = []
    if episodes.empty:
        return pd.DataFrame(columns=cols)
    for key, group in episodes.groupby(list(group_columns), sort=True):
        group = group.loc[~group.left_censored.astype(bool)]
        if group.empty:
            continue
        key = key if isinstance(key, tuple) else (key,)
        survival = 1.
        for duration, at_time in group.groupby("duration_ps", sort=True):
            risk = int((group.duration_ps >= duration).sum())
            dissociation = int((~at_time.right_censored.astype(bool)).sum())
            censored = int(at_time.right_censored.astype(bool).sum())
            survival *= 1. - dissociation/risk
            rows.append((*key, float(duration), risk, dissociation, censored, survival))
    return pd.DataFrame(rows, columns=cols)


def join_existing_sasa(counts, target_residues, sasa):
    """Match every constituent residue by exact frame time and topology index.

    Missing values are NaN, never zero. No interpolation, no rerun of FreeSASA.
    """
    cols = ["target_id", "topology_residue_index"]
    mapping = target_residues[cols].drop_duplicates()
    if mapping.empty or counts.empty or sasa.empty:
        return counts.assign(sasa_nm2=np.nan, sasa_coverage=0.)
    keys = ["time_ps", "topology_residue_index"]
    if sasa.duplicated(keys).any():
        raise ValueError("SASA has duplicate time/residue keys: select a single SASA run")
    unique_frames = counts[["target_id", "time_ps"]].drop_duplicates()
    required = unique_frames.merge(mapping, on="target_id", how="left",
                                    validate="many_to_many")
    paired = required.merge(sasa[keys+["sasa_nm2"]], on=keys, how="left",
                            validate="many_to_one")
    reduced = paired.groupby(["target_id","time_ps"], as_index=False).agg(
        sasa_nm2=("sasa_nm2", lambda x: x.sum(min_count=len(x))),
        sasa_coverage=("sasa_nm2", lambda x: float(x.notna().mean())))
    return counts.merge(reduced, on=["target_id", "time_ps"], how="left",
                        validate="many_to_one")



def calibration_stability(chunk_histograms, edges, *, min_chunks=4, tolerance_A=0.5):
    """Check RDF minimum consistency across chunk-level persisted histograms.

    Chunks are correlated trajectory segments, not independent replicas.
    This is a convergence diagnostic, not a confidence interval.
    """
    columns = ["solvent", "target_class", "site_class", "chunks_total",
               "chunks_resolved", "median_cutoff_A", "range_A", "status"]
    if chunk_histograms.empty:
        return pd.DataFrame(columns=columns)
    bin_cols = sorted((c for c in chunk_histograms if c.startswith("bin_")),
                      key=lambda c: int(c[4:]))
    from .solvation import rdf_from_counts
    output = []
    for key, group in chunk_histograms.groupby(
            ["solvent", "target_class", "site_class"], sort=True):
        estimates = []
        for _, part in group.iterrows():
            g = rdf_from_counts(part[bin_cols].to_numpy(float), edges,
                                float(part["exposure"]))
            rdf = pd.DataFrame({"solvent": key[0], "target_class": key[1],
                "site_class": key[2], "r_inner_A": edges[:-1],
                "r_outer_A": edges[1:], "g_bulk": g})
            candidate = suggest_shell_boundaries(rdf).suggested_cutoff_A.iloc[0]
            if np.isfinite(candidate):
                estimates.append(float(candidate))
        n, resolved = len(group), len(estimates)
        spread = float(np.ptp(estimates)) if estimates else np.nan
        status = ("insufficient_chunks" if n < min_chunks else
                  "unresolved_chunks" if resolved != n else
                  "unstable_minimum" if spread > tolerance_A else "stable")
        output.append((*key, n, resolved,
                       float(np.median(estimates)) if estimates else np.nan,
                       spread, status))
    return pd.DataFrame(output, columns=columns)


def aggregate_coordination(counts):
    """Derive class/time and target summaries; never sum per-site RDF counts.

    Each target contributes equally to per-class means. A target's
    multi-site waters are already de-duplicated within its molecule.
    This does not de-duplicate waters *across different molecules*.
    """
    class_keys = ["solvent", "target_class", "site_class"]
    class_columns = class_keys + ["time_ps", "target_count",
                                   "mean_n_solvent", "total_target_counts"]
    target_columns = ["solvent", "target_id", "target_class", "site_class",
                      "n_frames", "mean_n_solvent", "sd_n_solvent",
                      "hydrated_fraction"]
    if counts.empty:
        return (pd.DataFrame(columns=class_columns),
                pd.DataFrame(columns=target_columns))
    frame = counts.groupby(class_keys + ["time_ps"], as_index=False).agg(
        target_count=("target_id", "nunique"),
        mean_n_solvent=("n_solvent", "mean"),
        total_target_counts=("n_solvent", "sum"))
    by_target = counts.groupby(class_keys + ["target_id"], as_index=False).agg(
        n_frames=("n_solvent", "size"),
        mean_n_solvent=("n_solvent", "mean"),
        sd_n_solvent=("n_solvent", "std"),
        hydrated_fraction=("n_solvent", lambda x: float((x > 0).mean())))
    return frame[class_columns], by_target[target_columns]


def aggregate_coordination_parts(paths):
    """Bounded-memory aggregation from atomic count chunks.

    All frames within one saved part are complete, so per-class time
    means can be reduced independently. Per-target sufficient statistics
    are combined across chunks without ever loading all frames at once.
    """
    class_keys = ["solvent", "target_class", "site_class"]
    target_keys = class_keys + ["target_id"]
    class_tables = []
    partials = []
    for path in paths:
        counts = pd.read_parquet(path)
        if counts.empty:
            continue
        frame = counts.groupby(class_keys + ["time_ps"], as_index=False).agg(
            target_count=("target_id", "nunique"),
            mean_n_solvent=("n_solvent", "mean"),
            total_target_counts=("n_solvent", "sum"))
        class_tables.append(frame)
        numeric = counts.assign(
            positive=(counts["n_solvent"] > 0).astype(np.int64),
            sum_square=np.square(counts["n_solvent"].to_numpy(dtype=np.int64)))
        group = numeric.groupby(target_keys, as_index=False).agg(
            n_frames=("n_solvent", "size"),
            total=("n_solvent", "sum"),
            positive=("positive", "sum"),
            sum_square=("sum_square", "sum"))
        partials.append(group)
    empty = aggregate_coordination(pd.DataFrame())[0:2]
    if not partials:
        return empty
    class_frame = pd.concat(class_tables, ignore_index=True)
    summary = pd.concat(partials, ignore_index=True).groupby(
        target_keys, as_index=False).agg(
            n_frames=("n_frames", "sum"), total=("total", "sum"),
            positive=("positive", "sum"), sum_square=("sum_square", "sum"))
    n = summary.n_frames.to_numpy(float)
    total = summary.total.to_numpy(float)
    variance = (summary.sum_square.to_numpy(float)-total**2/n)/np.maximum(1., n-1.)
    summary["mean_n_solvent"] = total/n
    summary["sd_n_solvent"] = np.where(n > 1., np.sqrt(np.maximum(variance, 0.)), np.nan)
    summary["hydrated_fraction"] = summary.positive/n
    return class_frame, summary[
        ["solvent", "target_id", "target_class", "site_class",
         "n_frames", "mean_n_solvent", "sd_n_solvent", "hydrated_fraction"]]


def write_solvation_derivatives(directory, repout, *, sasa_run="default"):
    """Persist optional derived outputs from previous measurements only."""
    from pathlib import Path
    directory, repout = Path(directory), Path(repout)
    rdf = directory / "rdf.parquet"
    if rdf.exists():
        atomic_parquet(directory / "shell_candidates.parquet",
                       suggest_shell_boundaries(pd.read_parquet(rdf)))
    cache = repout / "cache" / "solvation"
    histogram_parts = sorted(cache.glob("hist-*.parquet"))
    if histogram_parts and rdf.exists():
        hist = pd.concat([pd.read_parquet(p) for p in histogram_parts],
                         ignore_index=True)
        full = pd.read_parquet(rdf)
        axis = full[["r_inner_A", "r_outer_A"]].drop_duplicates().sort_values("r_inner_A")
        if not axis.empty:
            edges = np.r_[axis.r_inner_A.to_numpy(float), float(axis.r_outer_A.iloc[-1])]
            atomic_parquet(directory / "shell_stability.parquet",
                           calibration_stability(hist, edges))
    counts = directory / "counts.parquet"
    count_parts = sorted(cache.glob("counts-*.parquet"))
    if counts.exists():
        by_class, by_target = aggregate_coordination_parts(count_parts)
        atomic_parquet(directory / "species_timeseries.parquet", by_class)
        atomic_parquet(directory / "target_summary.parquet", by_target)
    episodes = directory / "episodes.parquet"
    if episodes.exists():
        atomic_parquet(directory / "survival.parquet",
                       kaplan_meier(pd.read_parquet(episodes)))
    sasa = repout / "profiles/sasa"
    if sasa_run != "default":
        sasa = sasa / sasa_run
    sasa = sasa / "residue_sasa.parquet"
    if all(p.exists() for p in (directory/"counts.parquet",
                                directory/"target_residues.parquet", sasa)):
        atomic_parquet(directory / "hydration_sasa.parquet",
            join_existing_sasa(
                pd.read_parquet(directory/"counts.parquet"),
                pd.read_parquet(directory/"target_residues.parquet"),
                pd.read_parquet(sasa)))
