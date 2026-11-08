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
        durations = group.assign(right_censored=group.right_censored.astype(bool)).groupby(
            "duration_ps", sort=True).right_censored.agg(['size', 'sum'])
        risk = durations['size'].iloc[::-1].cumsum().iloc[::-1].to_numpy(int)
        censored = durations['sum'].to_numpy(int)
        dissociation = durations['size'].to_numpy(int)-censored
        survival = np.cumprod(1.-dissociation/risk)
        rows.extend((*key, float(duration), int(n), int(d), int(c), float(s))
            for duration, n, d, c, s in zip(durations.index, risk, dissociation, censored, survival))
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
                                   "mean_n_solvent", "total_target_counts",
                                   "hydrated_fraction"]
    target_columns = ["solvent", "target_id", "target_class", "site_class",
                      "n_frames", "mean_n_solvent", "sd_n_solvent",
                      "hydrated_fraction"]
    if counts.empty:
        return (pd.DataFrame(columns=class_columns),
                pd.DataFrame(columns=target_columns))
    frame = counts.assign(hydrated=(counts.n_solvent > 0).astype(float)).groupby(
        class_keys + ["time_ps"], as_index=False).agg(
        target_count=("target_id", "nunique"),
        mean_n_solvent=("n_solvent", "mean"),
        total_target_counts=("n_solvent", "sum"),
        hydrated_fraction=("hydrated", "mean"))
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
        frame = counts.assign(hydrated=(counts.n_solvent > 0).astype(float)).groupby(
            class_keys + ["time_ps"], as_index=False).agg(
            target_count=("target_id", "nunique"),
            mean_n_solvent=("n_solvent", "mean"),
            total_target_counts=("n_solvent", "sum"),
            hydrated_fraction=("hydrated", "mean"))
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
    bridges_path = directory/'bridges.parquet'
    frame_parts = sorted(cache.glob('frames-*.parquet'))
    if bridges_path.exists() and frame_parts:
        frame_times = pd.concat([pd.read_parquet(p) for p in frame_parts])['time_ps'].drop_duplicates().sort_values().to_numpy(float)
        bridges = pd.read_parquet(bridges_path)
        for gap, name in ((0, 'bridge_residence'), (1, 'bridge_residence_intermittent')):
            atomic_parquet(directory/f'{name}.parquet', bridge_residence(bridges, frame_times, gap))
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



def bridge_residence(bridges, frame_times, gap_frames=0):
    """Same-solvent simultaneous target coordination, not chemical binding.

    Derive pair-specific episodes from existing sparse bridge records. Only
    observed target pairs enter this table; no observed pair is not proof of
    zero occupancy for an unconfigured or untracked target pair.
    """
    from ..mechanisms.motif_residence import summarize_lipid_residence
    columns = ['target_a','target_b']
    if bridges.empty:
        return pd.DataFrame(columns=columns + ['species','molecule_id','duration_ps',
            'observed_span_ps','left_censored','right_censored'])
    pairs = bridges[columns].drop_duplicates().sort_values(columns).reset_index(drop=True)
    pairs['motif'] = [f'bridge_{i}' for i in range(len(pairs))]
    events = bridges.merge(pairs, on=columns, validate='many_to_one').rename(
        columns={'solvent':'species', 'solvent_id':'molecule_id'})
    events = events[['time_ps','motif','species','molecule_id']].drop_duplicates()
    episodes = summarize_lipid_residence(events, frame_times, gap_frames=gap_frames)
    return episodes.merge(pairs, on='motif', validate='many_to_one')


def bridge_pair_occupancy(bridges, frame_times, eligible_pairs, solvents):
    """Any-same-ion occupancy over *all* sampled frames for explicit pairs.

    Different ions at the two sites never constitute a bridge. Zero means a
    configured, observed pair had no same-ion event; it is not missing data.
    Fractions count sampled frames, not time-weighted occupation or affinity.
    """
    times = np.asarray(frame_times, float)
    if not len(times) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError('Bridge occupancy requires a finite increasing observation axis')
    pairs = [tuple(sorted(pair)) for pair in eligible_pairs]
    if any(len(pair) != 2 or pair[0] == pair[1] for pair in pairs) or len(set(pairs)) != len(pairs):
        raise ValueError('Eligible bridge pairs must be distinct non-self target pairs')
    data = bridges.copy()
    if not data.empty:
        if not data.time_ps.isin(times).all():
            raise ValueError('Bridge event outside observation times')
        data['pair'] = [tuple(sorted(pair)) for pair in data[['target_a','target_b']].itertuples(index=False, name=None)]
        data = data.drop_duplicates(['time_ps','solvent','solvent_id','pair'])
    rows = []
    for solvent in solvents:
        for pair in pairs:
            subset = data.loc[data.solvent.eq(solvent) & data.pair.map(lambda value: value == pair)] if not data.empty else data
            n_events = len(subset)
            n_frames = int(subset.time_ps.nunique()) if n_events else 0
            rows.append((solvent, *pair, len(times), n_frames, n_frames/len(times),
                         n_events/len(times), int(subset.solvent_id.nunique()) if n_events else 0))
    return pd.DataFrame(rows, columns=['solvent','target_a','target_b','observed_frames',
        'bridge_frames','sampled_bridge_fraction','mean_simultaneous_bridging_ions', 'unique_bridging_ions'])


def contact_records(paths):
    """Read only contact identity columns in bounded Parquet batches."""
    import pyarrow.parquet as pq
    columns = ['frame', 'solvent', 'target_id', 'solvent_id']
    for path in paths:
        for batch in pq.ParquetFile(path).iter_batches(batch_size=65536, columns=columns):
            data = batch.to_pydict()
            yield from zip(*(data[column] for column in columns))


def write_contact_episodes(event_paths, observed_frames, destination):
    """Stream frame-ordered contacts into sampled episodes atomically.

    State contains at most two successive contact frames plus an output batch.
    Episodes are written in completion order; identity/duration semantics are
    identical to the in-memory reference, which also accepts unordered events.
    Survival and other downstream finalizers have separate memory contracts.
    """
    from pathlib import Path
    import tempfile
    import pyarrow as pa
    import pyarrow.parquet as pq
    from ..core.residence import sampled_episode_durations
    from .solvation import EPISODE_COLUMNS
    destination = Path(destination)
    frames = [int(frame) for frame, _ in observed_frames]
    times = np.asarray([time for _, time in observed_frames], float)
    if len(frames) != len(set(frames)):
        raise ValueError('Duplicate observation frames')
    # Shared duration kernel validates the time axis before touching outputs.
    if len(times):
        index = np.arange(len(times))
        _, following = sampled_episode_durations(times, index, index)
    order = {frame: i for i, frame in enumerate(frames)}
    active, rows = {}, []
    previous = -1
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix='.parquet', delete=False) as handle:
        temporary = Path(handle.name)
    writer = None
    def flush():
        nonlocal writer
        if not rows:
            return
        data = pd.DataFrame(rows, columns=EPISODE_COLUMNS)
        first = np.asarray([order[int(frame)] for frame in data.start_frame], int)
        last = np.asarray([order[int(frame)] for frame in data.end_frame], int)
        data['duration_ps'] = times[last]-times[first]
        data['sampled_duration_ps'] = data['duration_ps'].to_numpy(float)+following[last]
        table = pa.Table.from_pandas(data, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(temporary, table.schema)
        writer.write_table(table)
        rows.clear()
    def finish(key, bounds):
        first, last = bounds
        rows.append((*key, frames[first], frames[last], np.nan,
            first == 0, last == len(frames)-1, np.nan,
            'next_sample_interval_or_terminal_median'))
        if len(rows) >= 4096:
            flush()
    try:
        for frame, solvent, target, entity in contact_records(event_paths):
            if frame not in order:
                raise ValueError('Contact event outside observation frames')
            index = order[frame]
            if index < previous:
                raise ValueError('Persisted contact events must be ordered by observation frame')
            if index != previous:
                for key in [key for key, (_, last) in active.items() if last < index-1]:
                    finish(key, active.pop(key))
                previous = index
            key = (solvent, target, entity)
            first, _ = active.get(key, (index, index))
            active[key] = (first, index)
        for key, bounds in active.items():
            finish(key, bounds)
        flush()
        if writer is not None:
            writer.close()
            writer = None
        else:
            pd.DataFrame(columns=EPISODE_COLUMNS).to_parquet(temporary, index=False)
        temporary.replace(destination)
    finally:
        try:
            if writer is not None:
                writer.close()
        finally:
            temporary.unlink(missing_ok=True)


def bridge_pair_evidence(event_paths, observed_frames, eligible_pairs, solvents):
    """Sampled same-entity traces and a deterministic temporal shift reference.

    One common offset shifts target-B contacts for all entities, preserving
    their identities and target-specific occupancy patterns. This breaks
    contemporaneous alignment, not spatial crowding, and is not a chemical
    cooperativity test. Only explicitly configured targets are retained.
    """
    from collections import defaultdict
    frames = [int(frame) for frame, _ in observed_frames]
    times = np.asarray([time for _, time in observed_frames], float)
    if (not len(times) or len(set(frames)) != len(frames)
            or not np.isfinite(times).all() or np.any(np.diff(times) <= 0)):
        raise ValueError('Bridge evidence requires unique observations and increasing finite times')
    order = {frame: i for i, frame in enumerate(frames)}
    pairs = [tuple(sorted(pair)) for pair in eligible_pairs]
    targets = {target for pair in pairs for target in pair}
    contacts = defaultdict(lambda: defaultdict(set))
    for frame, solvent, target, entity in contact_records(event_paths):
        if frame not in order:
            raise ValueError('Contact event outside bridge observation frames')
        if target in targets and solvent in solvents:
            contacts[(solvent, target)][entity].add(order[frame])
    regular = bool(len(times) >= 4 and np.allclose(np.diff(times), times[1]-times[0], rtol=1e-5, atol=1e-8))
    offsets = np.unique(np.linspace(1, len(times)-1, min(100, len(times)-1), dtype=int)) if regular else []
    traces, references, summaries = [], [], []
    for solvent in solvents:
        for left, right in pairs:
            a, b = contacts[(solvent, left)], contacts[(solvent, right)]
            common = sorted(a.keys() & b.keys())
            counts = np.zeros(len(times), int)
            for entity in common:
                counts[list(a[entity] & b[entity])] += 1
            present = counts > 0
            traces.extend((float(time), solvent, left, right, int(count), int(hit))
                for time, count, hit in zip(times, counts, present))
            null = []
            for offset in offsets:
                shifted_counts = np.zeros(len(times), int)
                for entity in common:
                    overlap = np.intersect1d(np.fromiter(a[entity], int),
                        (np.fromiter(b[entity], int)+offset) % len(times), assume_unique=True)
                    shifted_counts[overlap] += 1
                fraction = float(np.mean(shifted_counts > 0))
                null.append(fraction)
                references.append((solvent, left, right, int(offset), fraction,
                                   float(shifted_counts.mean())))
            lo, hi = np.quantile(null, [.025, .975]) if null else (np.nan, np.nan)
            mean = float(np.mean(null)) if null else np.nan
            summaries.append((solvent, left, right, float(present.mean()), mean,
                float(lo), float(hi), float(present.mean()-mean), len(offsets),
                'descriptive_shift_reference' if regular else
                'unavailable_short_window' if len(times) < 4 else 'unavailable_irregular_sampling'))
    return (pd.DataFrame(traces, columns=['time_ps','solvent','target_a','target_b',
                'bridge_entity_count','bridge_present_fraction']),
            pd.DataFrame(references, columns=['solvent','target_a','target_b','offset_frames',
                'sampled_bridge_fraction','mean_simultaneous_entities']),
            pd.DataFrame(summaries, columns=['solvent','target_a','target_b',
                'observed_bridge_fraction','shift_reference_mean', 'shift_reference_low',
                'shift_reference_high','observed_minus_shift_mean','n_offsets','status']))
