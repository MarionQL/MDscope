# Neuronal interface analyses

The canonical architecture and package locations are documented in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md). The rules for restructuring
are in [docs/REPOSITORY_RULES.md](docs/REPOSITORY_RULES.md).

This document focuses on the neuronal-interface configuration and scientific
interpretation of optional analyses, not an additional module map.

## Adding new analyses

1. Decide the output shape **before** writing code:
   - scalar time series: `time_ps` plus one column per quantity, with optional `group`;
   - frame-by-residue/profile: positions, labels, time-indexed matrix and reducer;
   - interaction: event key, event frame/time, identity and occupancy;
   - independent per-replicate parameter (e.g., diffusion slope): one row per *replicate*, not per frame.
2. Add a module in its appropriate `analysis/` subpackage with a top-level `run(rep,cfg,repout)` function for trajectory-specialized analyses. Use `MDAnalysis` and shared geometry tools; vectorize expensive distance calculations.
3. Register the module in `analysis/registry.py` within `specialized_tasks`. Persist scalar observables as `extensions/<analysis_name>/timeseries.parquet`, and register their file in `additional_time_series`.
4. This output feeds the *existing* `summary/timeseries.parquet`, `analysis/statistics.py` and `analysis/plotting.py`. **Do not implement custom hypothesis testing or frame-pooled error bars** in individual modules.
5. Profile analyses should expose arrays to `analysis/statistics.py:collect_temporal_profiles`; interaction events should use `analysis/interaction_tables.py`. Add a small dataset test and import test before production use.
6. Expensive specialized analyses should operate as separate restartable passes. Cheap per-frame observables should share `streaming.py` so large XTC files are not opened once per scalar observable.

## Diffusion coefficients from existing MSD

Keep the original YAML section:

```yaml
analyses:
  membrane:
    enabled: true
    axis: z
    headgroup_selection: "name P"
    msd:
      enabled: true
      selection: null
      include_3d: false
      slope_start_fraction: 0.2
      slope_end_fraction: 0.8
```

No `lipid_diffusion` analysis is needed or supported by the modular extension. For every replicate, `observables/msd.py` reads the already consolidated `membrane_msd_lateral_A2` series and fits the linear regime. The factor is **2.5** to convert a lateral slope measured in Å²/ps into nm²/ns, using `MSD = 4Dt`.

Outputs:
- `<system>/<replicate>/summary/diffusion_coefficients.parquet`: D, fit interval, slope and R².
- `project_summary/diffusion/replicate_coefficients.parquet`: all replicate-level coefficients.
- `project_summary/diffusion/system_coefficients.parquet`: descriptive summaries, *no significance claims*.
- `project_summary/diffusion/diffusion_coefficients.png`: independent replicate points and means.

This remains a **single-origin** MSD estimator. It is not species-resolved or a multi-time-origin diffusion estimator; Multi-origin MSD is implemented separately in `observables/multi_origin_msd.py` and `observables/multi_origin_pipeline.py`; choose the estimator explicitly. Fit-window sensitivity and equilibrated trajectory selection matter.

## Experimental neuronal-interface extensions

```yaml
analyses:
  local_environment:
    enabled: true
    protein_selection: "protein"
    anchor_selection: "protein and name CA and (resid 151-170 or resid 349-368)"
    lipid_selection: "resname POPC POPE POPI DOPC DPPC PSM CHL1 BGAL GM1a CER160"
    glycan_selection: "resname AGLCNA AIDOA BGLCA BGAL BXYL"
    ion_selection: "resname SOD CLA POT"
    annulus_A: 12.0
    contact_cutoff_A: 5.0
    motifs:
      nxip: "N.IP"
      tm_gxxxg: "G...G"
  coupling_screen:
    enabled: true
    block_frames: 50
    include_defaults: true
statistics:
  hypothesis_tests:
    enabled: false
```

Validate selections for each topology before production. Local membrane enrichment requires an orthorhombic simulation cell in this version. The motif lookup uses the topology sequence and is chain-local; motif contacts are *not* the same as motif SASA or binding capacity. Correlations are descriptive, not causal.

## Performance, data integrity and testing

- The common trajectory pass shares I/O for RMSD, RMSF, Rg, density, original membrane MSD and inexpensive contacts.
- Specialized readers are individually checkpointed by analysis name; result profiles preserve the existing Parquet/NPZ/npy contract.
- The code/configuration hash invalidates prior checkpoints. **A mismatch archives previous outputs under `.archives/` instead of deleting them.** Changing code can still trigger complete re-analysis, so use a distinct output directory if you wish to retain old figures without creating a new run.
- `python -m compileall -q analysis tests` checks syntax; `pytest -q tests` runs synthetic regression checks.
- CI runs these checks on pushes and pull requests. End-to-end tests against a real protein/lipid trajectory are still required before manuscript use.

For usage and prior YAML options, see `README.md`.
