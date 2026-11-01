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
