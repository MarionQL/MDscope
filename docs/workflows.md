# Workflows and saved results

## Run a full SDC4 analysis

`synd4_full_mechanisms.yaml` configures two membrane conditions and **six replicas**; its output directory is separate from `synd4_baseline.yaml`.

```bash
python -m pip install -e '.[sasa]'
mdscope validate synd4_full_mechanisms.yaml
mdscope run synd4_full_mechanisms.yaml

# On an array/HPC run, process replicas first; then:
mdscope finalize synd4_full_mechanisms.yaml

# Recalculate from saved replica outputs, no XTC reread:
mdscope statistics synd4_full_mechanisms.yaml
mdscope plotting synd4_full_mechanisms.yaml
```

- `run`: replica processing and numerical finalization.
- `replicate --index N`: one replica, no project finalization.
- Saved-data recomputation: [commands below](#recompute-statistics-and-figures-without-rereading-trajectories).

## Environment variables

Scientific settings belong in YAML; these variables control logging, resources, environment setup and test tooling.

| Variable | Reader | Default when unset | Effect |
|---|---|---|---|
| `MDSCOPE_DEBUG` | CLI | Disabled | `1` re-raises expected user-input errors with traceback |
| `MDSCOPE_PROGRESS_EVERY_PERCENT` | `core/progress.py`, Slurm replicate runner | `1` | Percentage steps for progress when total work is known |
| `MDSCOPE_PROGRESS_EVERY_FRAMES` | `core/progress.py` | `25` | Frame steps for progress with unknown total |
| `MDSCOPE_PROFILE` | `core/profiling.py`, interaction statistics; Slurm and real-data scripts | Off in Python; `1` in replicate/real-data scripts | Collect wall-time profiling data |
| `MDSCOPE_STATS_WORKERS` | `statistics/common.py`, statistics Slurm script | `SLURM_CPUS_PER_TASK` when available; otherwise module default | Worker count for supported statistics operations |
| `MDSCOPE_ORDER_FINALIZE_MEMORY_MB` | `pipeline/fused_membrane.py` | `256` | Memory budget (MiB) for lipid-order finalization |
| `MDSCOPE_ENV_SETUP` | `slurm/job_env.sh` | ARC example Miniconda `mdanalysis` environment | Path to a shell setup script; `none` uses inherited environment |
| `MDSCOPE_REAL_DATA` | `scripts/run_real_data_checks.sh` | `examples/data` | Real test-data directory |
| `MDSCOPE_REAL_OUTPUT` | `scripts/run_real_data_checks.sh` | `real_test_results` | Destination for real-data checks |
| `MDSCOPE_HBOND_BENCH_FRAMES` | Real-data H-bond benchmark via `scripts/run_real_data_checks.sh` | `8` | H-bond benchmark frame count |
| `MDSCOPE_BOX_FRAMES` | `scripts/run_real_data_checks.sh` | `40` (box sweep), `200` (clearance) | Frame cap for optional box-geometry checks |
| `PYTHONUNBUFFERED` | Slurm and real-data shell scripts / Python runtime | `1` in scripts | Flush Python log output promptly |

`SLURM_CPUS_PER_TASK` and `SLURM_ARRAY_TASK_*` are scheduler-provided
settings, not MDscope-specific environment variables.

## Slurm jobs

- Submit `slurm/{run,finalize,statistics,plotting}_mdscope.sh` from the checkout with validated YAML; each sources `slurm/job_env.sh`.
- Install once in the job environment with `python -m pip install -e '.[sasa]'`; jobs never install packages.
- Jobs verify importability and log version, source path and commit; setup failures exit with instructions.

| Setting | ARC example (script default) | Other sites |
|---|---|---|
| Account, partition | `#SBATCH -A bevanlab`, `-p normal_q` | `sbatch -A ACCOUNT -p PARTITION ...`; command-line flags override `#SBATCH` lines |
| Environment | Miniconda env `mdanalysis` | `MDSCOPE_ENV_SETUP=/path/setup.sh` (bash file sourced by the job) or `MDSCOPE_ENV_SETUP=none` (keep the environment exported by `sbatch`) |
| Replica array | `--array=0-5` (six replicas of `synd4_full_mechanisms.yaml`) | `--array=0-<N-1>`, N from the schema-only plan |

`slurm/run_mdscope.sh` runs one task per YAML replicate, indexed in plan order; the plan reads no topology or trajectory:

```bash
N=$(mdscope plan CONFIG.yaml --format json | python -c 'import json,sys; print(len(json.load(sys.stdin)["replicates"]))')
sbatch -A ACCOUNT -p PARTITION --array=0-$((N-1)) slurm/run_mdscope.sh CONFIG.yaml
```

Every task repeats this derivation and exits before analysis, printing the expected `--array` range, unless the submitted array is exactly `0..N-1` (`SLURM_ARRAY_TASK_MIN/MAX/COUNT`). Reprocess one failed replica with `mdscope replicate CONFIG.yaml --index I` in its own job, not with a partial array.

The statistics job reports its current stage and elapsed time in `slurm-out/mdscope/mdscope_statistics-<jobid>.out`. Long stages print a heartbeat every minute. Set `progress.enabled: false` to suppress progress messages.

Metric-level inference controls: [Configuration](configuration.md#selective-metric-inference).

## Data flow

| Stage | Operations | Extra XTC pass? |
|---|---|---|
| Shared source loop | RMSD/RMSF/Rg, bilayer geometry/APL/MSD, lipid order/normals/enrichment, density, contacts, H bonds, serial DSSP, SASA, SDC4, residue distances | No extra pass for enabled fused observers |
| Optional/specialized | Multi-origin MSD, potassium bridge audit, some DSSP backends | **Yes**, independently restartable |
| `summarize_replicate` | Merge framewise output; compute profiles and interaction derivatives | No |
| `statistics`/`finalize` | Replica/bootstrap comparisons, RDF/diffusion/coupling, shared states, coverage and findings | No |
| `plotting`/`report` | Render saved tables and figures | No |

**Source frames remain unrotated by protein fitting.** Protein H bonds/DSSP reuse MDAnalysis private `_single_frame` interfaces: rerun real-data parity tests after MDAnalysis upgrades. Per-leaflet lipid order fixes initial cohorts; normals/enrichment recalculate leaflets. See [Scientific contracts](developer/scientific-contracts.md).

## Restart and file retention

- Parts: `<replica>/cache/fused_{membrane,extensions,distance}/`, with numbered Parquet/NPZ files.
- Committed parts are validated by identity, metadata/inventory and required companion files; incomplete or corrupt products fail safely.
- Only **uncommitted** parts are pruned/replayed. Keep restart caches until the stream completion marker is written.
- `fused_*_timings.json` and `stage_timings.json` record where time is spent.
- Successful numerical finalization archives the previous project summary in `.archives/`; failed preparation preserves it. An unfinished attempt is not current publication: rerun finalization before downstream reuse. Archives can consume disk.
- Run one finalization/plotting writer per project. Completion markers alone do not establish scientific validity.

## SDC4 topology checkpoints

- Audit both topologies: water/ions, lipid anchors, HS attachments, configured motifs and analysis frames. Release-specific identities and measurements: [SDC4 analyses](analyses.md#syndecan-4-sdc4-mechanisms).
- The YAML owns the enabled measurement list and atom selectors. Preflight audits both GGIVG TM sequences through the same resolver as runtime; `require_gxxxg: true` fails early when either is missing. Validate rebuilt topologies.
- Geometry/SASA assumptions: [Scientific contracts](developer/scientific-contracts.md).
- Adding fused measurements to completed runs replays only missing outputs in one shared pass. Re-finalize statistics afterwards.

## Find and interpret results

| Path under `project.output` | Read it for |
|---|---|
| `<system>/<replica>/summary/` | Scalar metrics, diffusion and per-replica summaries |
| `<system>/<replica>/profiles/` | RMSF, density, SASA, membrane order, residue distances |
| `<system>/<replica>/interactions/` | Events, occupancy, residence/lifetimes |
| `<system>/<replica>/extensions/` | Solvation, lipid neighbors/footprint, SDC4/HS and K audits |
| `project_summary/statistics/findings.md` | Fast review of valid interval exclusions and coverage |
| `statistics/findings.parquet` | Exact finite contrasts excluding zero (exploratory until multiplicity checked) |
| `statistics/coverage/contrast_coverage.parquet` | All registered comparisons, including unavailable status |
| `statistics/output_inventory.parquet` | Saved products and Parquet row counts |
| `project_summary/shared_states/` | Model, frame assignments and replica populations |

```bash
mdscope report /path/to/project_summary/statistics
mdscope scan-residue-regions /path/to/project_summary/statistics
```

- Region scan: ≥3 consecutive residues on the same chain, finite in **every** expected replica, with nonoverlapping replica ranges and consistent direction; requires ≥2 replicas per condition. It is exploratory, **not** a simultaneous CI.
- RDF-shell contrasts average saved curves using shell volumes; multi-origin diffusion contrasts refit shared supported lag windows. Require complete replica coverage.
- Matched scalar windows: `statistics.matched_window_ps: auto` or `[start,stop]`. Coupling requires a common regular grid. HS episodes retain censoring.
- Relative membrane position and lipid density maps have distinct reference/area conventions; do not infer binding affinity from proximity.
- Three independent replicas per condition give limited permutation resolution; inspect temporal stability and each replicate before interpreting biological contrasts.

## Recompute statistics and figures without rereading trajectories

When replica products are complete, statistics, residue-response regions, lipid-contact matrices and figures are recomputed from saved data only:

```bash
mdscope statistics analysis.yaml            # archives the previous statistics export first
mdscope plotting analysis.yaml
mdscope scan-residue-regions RESULTS/project_summary/statistics   # residue-response tables from a complete export
mdscope statistics analysis.yaml --results-dir PATH   # finalize results stored elsewhere
```

Changing `analyses.lipid_contacts.classes` or `statistics.residue_response` settings affects only these offline products, not the per-frame observers. Rerun the replica pipeline only for genuinely missing measurement products. Definitions of these products: [Analyses](analyses.md#residue-response-heatmap) and [lipid contacts](analyses.md#lipid-contact-classes-and-statistics).

[Configuration](configuration.md) · [Analysis definitions](analyses.md) · [Scientific audit](developer/scientific-audit.md)
