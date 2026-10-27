# MDscope: modular analysis architecture

This branch reorganizes the existing software; the CLI and existing YAML keys remain valid. **RMSD and RMSF are separate calculations**, and **MSD is not recalculated by a second reader**.

## Source layout

| Module | Responsibility |
|---|---|
| `analysis/cli.py` | Stable CLI: init, validate, run, replicate, finalize |
| `analysis/config.py` | Resolved YAML, paths, selection and statistics validation |
| `analysis/runner.py` | Restart-safe orchestration, consolidation, project finalization |
| `analysis/registry.py` | Registered specialized analyses and common scalar output discovery |
| `analysis/streaming.py` | **One** streaming trajectory coordinator for inexpensive frame observables |
| `analysis/alignment.py` | Shared coordinate fitting and alignment validation |
| `analysis/rmsd.py` | RMSD only |
| `analysis/rmsf.py` | Streaming RMSF accumulator, aligned-coordinate profile reader |
| `analysis/radius_of_gyration.py` | Radius of gyration |
| `analysis/membrane.py` | APL, thickness and membrane state |
| `analysis/msd.py` | **Original MSD**, PBC incremental unwrapping and *descriptive* diffusion coefficient fits |
| `analysis/density_1d.py` | Density accumulation and bootstrappable profiles |
| `analysis/membrane_order.py` | Existing acyl-chain orientational order |
| `analysis/dssp.py` | DSSP assignments |
| `analysis/sasa.py` | FreeSASA calculations |
| `analysis/structure_profiles.py` | Aggregate DSSP/SASA into the shared profile/time-series schema |
| `analysis/residue_membrane_distance.py` | Exact blocked, PBC-aware residue distances |
| `analysis/salt_bridges.py` | Salt bridge events |
| `analysis/hydrogen_bonds.py` | Hydrogen bond events |
| `analysis/ligand_contacts.py` | Ligand interaction events |
| `analysis/interaction_tables.py` | **Shared** interaction occupancy and lifetime tables |
| `analysis/local_environment.py` | Experimental sequence motif and local lipid/HS observables |
| `analysis/biophysical_common.py` | Selection-independent neighbor, geometry and motif utilities |
| `analysis/coupling_screen.py` | Exploratory within-replicate biochemical associations |
| `analysis/statistics.py` | **One** framework for time-series, profile and interaction replicate-aware statistics |
| `analysis/plotting.py` | **One** figure framework for each canonical data structure |
| `analysis/timeseries.py`, `analysis/profiles.py`, `analysis/interactions.py` | Thin compatibility facades preserving previous imports |

## Adding new analyses

1. Decide the output shape **before** writing code:
   - scalar time series: `time_ps` plus one column per quantity, with optional `group`;
   - frame-by-residue/profile: positions, labels, time-indexed matrix and reducer;
   - interaction: event key, event frame/time, identity and occupancy;
   - independent per-replicate parameter (e.g., diffusion slope): one row per *replicate*, not per frame.
2. Add an individual `analysis/<analysis_name>.py` with a top-level `run(rep,cfg,repout)` function for trajectory-specialized analyses. Use `MDAnalysis` and shared geometry tools; vectorize expensive distance calculations.
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

No `lipid_diffusion` analysis is needed or supported by the modular extension. For every replicate, `msd.py` reads the already consolidated `membrane_msd_lateral_A2` series and fits the linear regime. The factor is **2.5** to convert a lateral slope measured in Å²/ps into nm²/ns, using `MSD = 4Dt`.

Outputs:
- `<system>/<replicate>/summary/diffusion_coefficients.parquet`: D, fit interval, slope and R².
- `project_summary/diffusion/replicate_coefficients.parquet`: all replicate-level coefficients.
- `project_summary/diffusion/system_coefficients.parquet`: descriptive summaries, *no significance claims*.
- `project_summary/diffusion/diffusion_coefficients.png`: independent replicate points and means.

This remains a **single-origin** MSD estimator. It is not species-resolved or a multi-time-origin diffusion estimator; to add those, extend the *same* MSD module rather than implementing a second overlapping trajectory analysis. Fit-window sensitivity and equilibrated trajectory selection matter.

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
- `python -m compileall -q analysis tests` checks syntax; `pytest -q tests` runs synthetic and compatibility-import checks.
- CI runs these checks on pushes and pull requests. End-to-end tests against a real protein/lipid trajectory are still required before manuscript use.

For usage and prior YAML options, see `README.md`.
