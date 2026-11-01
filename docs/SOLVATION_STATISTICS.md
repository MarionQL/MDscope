## Solvation statistics: reuse existing moving-block and replicate hierarchy

The per-replicate solvation observer already writes:
- `extensions/solvation/counts.parquet`: per-frame, per-target unique solvent counts with solvent species, target IDs/classes and time.
- `extensions/solvation/species_timeseries.parquet`: per-frame chemical-class average coordination and count; now also `hydrated_fraction`, the fraction of configured targets with ≥1 solvent entity within its coordination shell.
- `extensions/solvation/target_summary.parquet`: per-residue or per-molecule mean, standard deviation and hydrated-frame fraction, descriptive **within a replicate**.
- `extensions/solvation/rdf.parquet`: bulk-number-density-normalized radial distribution by solvent and target/site class, with saved radial bins and raw pair counts.
- `extensions/solvation/{events,bridges,water_hbond_bridges,episodes,survival}.parquet` when the relevant opt-in modes are enabled; other outputs include `shell_candidates`, `shell_stability`, `target_residues`, optional `hydration_sasa`, and optional sparse `water_density_3d`.

### What has changed

`statistics.time_series.collect_time_series` now discovers **existing** solvation
`species_timeseries.parquet` alongside other MDscope time series. No new
resampling method or separate statistical implementation is added.
Each stable key is `solvation/{solvent}/{target_class}/{site_class}`,
with metrics `mean_n_solvent` and `hydrated_fraction`, where available.

With `statistics.enabled: true`, the **existing** framework calculates
FFT-based autocorrelation, correlation-informed temporal block lengths,
within-replicate moving-block bootstrap, between-replicate hierarchical
bootstrap, replicate means/CI, and difference-of-conditions CI.
Existing permutation/FDR settings apply; **with only 3 independent
replicates per condition, do not overinterpret p-values**.

Output tables, in the existing project-wide structure:
- `project_summary/statistics/time_series/replicate_statistics.parquet`
- `project_summary/statistics/time_series/system_statistics.parquet`
- `project_summary/statistics/time_series/system_comparisons.parquet`
- `project_summary/statistics/time_series/autocorrelation.parquet`
- `<replicate>/summary/{time_series_statistics,autocorrelation}.parquet`

All solvation rows are identifiable with `group` beginning `solvation/`.
Comparison difference convention is **system_2 minus system_1**. Use named
system columns and the actual direction; do not assume a sign convention.

RDFs are **curves**, not independent scalar samples. We therefore also output:
- `project_summary/statistics/solvation/replicate_rdf.parquet`: raw per-replicate curves.
- `project_summary/statistics/solvation/system_rdf.parquet`: pointwise
  replicate mean and SD, number of contributing replicates, and coverage flag.
  This is **descriptive** replicate variation, not per-bin confidence intervals
  or hypothesis tests. The RDF first-minimum cutoffs are exploratory.

### Important limitations

1. Each simulation replicate is the independent unit. Frames and multiple lipid
   molecules in the same trajectory are **not** independent replicates.
2. The class mean weights each configured molecular target equally at every
   frame. Compare the same chemically meaningful target definitions across
   conditions; `target_count` is available to audit accidental composition changes.
3. `hydrated_fraction` indicates ≥1 selected solvent inside the fixed shell.
   It is **not** the probability of a specific water identity or an
   uninterrupted residence event. An ion-specific fraction uses that ion
   species' configured selection and cutoff.
4. Residence/survival and 3D density remain descriptive by default.
   We do not bootstrap censored residence curves or every voxel.
5. A preexisting `species_timeseries.parquet` lacking `hydrated_fraction`
   still supports bootstrap of `mean_n_solvent`; regenerate the saved class
   time series from raw counts if occupancy statistics are also required.
6. Using many chemically fine target classes increases the number of tests
   and bootstrap work. Default to scientifically relevant classes/regions,
   not every target atom.

No additional trajectory pass or extra neighbor search is necessary.
