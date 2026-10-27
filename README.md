# MDscope

MDscope is a compact MDAnalysis-based workflow for comparing replicated molecular-dynamics systems. It runs a shared streaming pass for inexpensive frame-wise observables, separate passes for analyses that require specialized algorithms, then applies the same replicate-aware statistical framework across time series, profiles, and interactions.

## Dependencies

Required for the main workflow:

```bash
pip install numpy pandas pyarrow pyyaml matplotlib MDAnalysis
```

SASA additionally requires FreeSASA:

```bash
pip install freesasa
```

To create a self-contained PDF evidence report and a compact review ZIP:

```bash
pip install reportlab pillow
```

`pyarrow` is recommended as the Parquet backend. `fastparquet` can also satisfy pandas Parquet output.

## Basic use

From the directory containing the `mdscope` package:

```bash
python -m mdscope.analysis.cli validate example_config.yaml
python -m mdscope.analysis.cli run example_config.yaml
```

Useful commands:

```bash
# Interactively create a minimal starter YAML.
python -m mdscope.analysis.cli init --output analysis.yaml

# Check files/selections and catch obvious configuration problems.
python -m mdscope.analysis.cli validate analysis.yaml

# Run every configured replicate, consolidate results, calculate statistics, and plot.
python -m mdscope.analysis.cli run analysis.yaml

# Run one zero-based systemxreplicate entry (for example from a job array).
python -m mdscope.analysis.cli replicate analysis.yaml --index 0

# Re-run only project-level statistics/plots from already persisted replicate outputs.
python -m mdscope.analysis.cli finalize analysis.yaml
```

The project output directory contains the resolved configuration and environment metadata so a completed analysis can be audited later.

For the full two-condition, six-replicate SDC4 mechanism analysis, see
[`synd4_full_mechanisms.yaml`](synd4_full_mechanisms.yaml) and
[`REPORT.md`](REPORT.md). The original system/trajectory paths are
unchanged; the expanded analysis writes to a separate results directory.

When `report.enabled: true`, finalization creates:

```text
project_summary/report/MDscope_evidence_report.pdf
project_summary/report/MDscope_review_bundle.zip
```

To regenerate just these files from existing finalized outputs (no XTC
re-read, no project summary archiving):

```bash
python -m mdscope.analysis.cli report synd4_full_mechanisms.yaml
```


# YAML configuration

See `example_config.yaml` for a complete example. Paths are resolved relative to the YAML file unless they are absolute.

## `project`

```yaml
project:
  name: syndecan_membrane
  output: ./mdscope_results
```

| Parameter | Meaning |
| --- | --- |
| `name` | Human-readable project name. |
| `output` | Root output directory. Relative paths are interpreted relative to the YAML file. |

## `trajectory`

```yaml
trajectory:
  start: null
  stop: null
  stride: 10
  chunk_frames: 1000
```

| Parameter | Meaning |
| --- | --- |
| `start` | First raw trajectory frame to analyze. `null` starts at the beginning. |
| `stop` | Python-slice-style exclusive stop frame. `null` runs to the end. |
| `stride` | Analyze every Nth raw frame. Must be at least 1. |
| `chunk_frames` | Number of analyzed frames accumulated before streamed tables/caches are flushed to disk. Larger values reduce write frequency; smaller values reduce RAM use and make restart checkpoints more frequent. |

The shared time-series/interaction pass opens the trajectory once and calculates several inexpensive observables together. This avoids repeatedly reading very large XTC/DCD files.

## `time`

```yaml
time:
  dt_ps: null
  origin_ps: 0.0
```

| Parameter | Meaning |
| --- | --- |
| `dt_ps` | Optional time spacing **between analyzed frames after `trajectory.stride`**. Leave `null` when the trajectory contains correct timestamps. |
| `origin_ps` | Time assigned to the first analyzed frame when `dt_ps` is supplied. |

For example, if the raw trajectory is saved every 10 ps and `stride: 10`, then `dt_ps` should be `100`, not `10`.

## `systems`

Each system has one topology and one or more trajectories:

```yaml
systems:
  plus_glycolipids:
    label: "+ Glycolipids"
    topology: ./plus/protein.pdb
    replicates:
      - trajectory: ./plus/rep1.xtc
        label: rep1
        plot_label: "Replicate 1"
      - trajectory: ./plus/rep2.xtc
        label: rep2
        plot_label: "Replicate 2"
```

| Parameter | Meaning |
| --- | --- |
| system key | Stable internal system ID used in output paths/tables. |
| `label` | Optional display label for plots. |
| `topology` | Topology used by every replicate in that system. |
| `replicates[].trajectory` | Trajectory path. |
| `replicates[].label` | Stable replicate ID. Use matching IDs across systems if using a paired statistical design. |
| `replicates[].plot_label` | Optional display label. |

## `selections`

Selections use normal MDAnalysis selection syntax.

```yaml
selections:
  protein: "protein"
  alignment: "protein and name CA"
  features: "protein and name CA"
  ligand: null
```

| Parameter | Meaning |
| --- | --- |
| `protein` | Main protein selection. Used as the default for several analyses. |
| `alignment` | Atoms used to align frames to the first-frame reference. |
| `features` | Atoms used for RMSD/RMSF and the optional aligned-coordinate cache. |
| `ligand` | Ligand selection used for ligand RMSD/contact analysis. Leave `null` when no ligand analysis is required. |

### Groups

Groups let the same structural observables be calculated for chains/domains/regions.

```yaml
selections:
  groups:
    chain_A:
      label: "Chain A"
      selection: "protein and chainID PROA"
      alignment_selection: "protein and chainID PROA and name CA"
```

| Parameter | Meaning |
| --- | --- |
| group key | Stable ID stored in tables. |
| `label` | Plot label. |
| `selection` | Atoms/residues belonging to the group. |
| `alignment_selection` | Optional group-specific alignment. If omitted, the globally aligned coordinates are used. |
| `linestyle` | Optional Matplotlib line style. A style is assigned automatically when omitted. |

---

# Analyses

Each analysis can be enabled with either `analysis_name: true` or a mapping containing `enabled: true` and its parameters.

## Alignment

```yaml
alignment:
  enabled: true
  selection: null
```

`selection: null` uses `selections.alignment`. Alignment centers and rigidly fits each frame to the first frame before RMSD/RMSF measurements. Analyses that must use the physical membrane coordinates still operate on the unaligned trajectory.

## RMSD

```yaml
rmsd:
  enabled: true
```

Calculates RMSD for `selections.features` after the configured alignment. Group RMSD is also calculated when groups are configured.

## RMSF

```yaml
rmsf:
  enabled: true
```

RMSF is accumulated online using Welford variance, so the ordinary RMSF calculation does not require every coordinate frame in RAM. When `cache.aligned_coordinates: true`, aligned feature coordinates are also saved in compressed chunks so RMSF can receive temporal block-bootstrap confidence intervals.

A PDB with equivalent RMSF B-factors is written when the topology supports B-factor output.

## Radius of gyration

```yaml
radius_of_gyration:
  enabled: true
  selection: "protein"
```

`selection` defaults to `selections.protein`. Group Rg is calculated from the intersection of the Rg selection and each group.

## 1D density

```yaml
density_1d:
  enabled: true
  axis: z
  bin_width_A: 1.0
  center_selection: "protein"
  selections:
    lipid_phosphates:
      label: "Lipid phosphates"
      selection: "name P"
      mode: number
```

| Parameter | Meaning |
| --- | --- |
| `axis` | `x`, `y`, or `z`. |
| `bin_width_A` | Target bin width in A. |
| `center_selection` | Optional selection used as the profile origin. If omitted, the box center is used. |
| `selections` | Named sets to profile. |
| `selection` | MDAnalysis atom selection for one profile. |
| `mode` | `number` for number density or `mass` for mass density. |
| `groups` | Optional per-target boolean; when true, group intersections are also profiled. |

The mean profile is accumulated online. Per-frame profiles are temporarily chunked as compressed `float32` arrays because the hierarchical profile bootstrap needs frame-resolved information without keeping it all in RAM.

### Named density runs

Density supports multiple named configurations:

```yaml
density_1d:
  enabled: true
  runs:
    membrane_z:
      axis: z
      bin_width_A: 1.0
      selections: {...}
    protein_x:
      axis: x
      bin_width_A: 2.0
      selections: {...}
```

## Membrane metrics

```yaml
membrane:
  enabled: true
  axis: z
  headgroup_selection: "name P"
  area_per_lipid:
    enabled: true
  thickness:
    enabled: true
  msd:
    enabled: true
    selection: "name P"
    include_3d: false
    slope_start_fraction: 0.2
```

| Parameter | Meaning |
| --- | --- |
| `axis` | Bilayer-normal axis. |
| `headgroup_selection` | Atoms used to classify lipid residues into upper/lower leaflets and calculate leaflet COMs. |
| `area_per_lipid.enabled` | Calculate box cross-sectional area divided by the number of lipids in each leaflet. |
| `thickness.enabled` | Calculate the periodic separation between upper/lower headgroup centers of mass. |
| `msd.enabled` | Calculate lipid/residue MSD. |
| `msd.selection` | Atoms used to define each tracked residue center. Defaults to the headgroup selection. |
| `msd.include_3d` | Also output full 3D MSD; lateral MSD is always calculated when MSD is enabled. |
| `msd.slope_start_fraction` | Start fraction of analyzed MSD interval used for the replica-level lateral diffusion slope fit; must be less than `slope_end_fraction`. |

Leaflet identity is assigned from the first analyzed configuration and kept fixed. MSD accumulates minimum-image step vectors so crossing a periodic boundary does not create a false displacement jump.

## Residue-to-membrane distance

```yaml
residue_membrane_distance:
  enabled: true
  selection: "protein"
  lipid_resnames: [POPC, POPE, POPI, CHL1, BGAL, GM1a]
  membrane_selection: null
  exclude_hydrogens: true
  expand_residue_atoms: true
  initial_search_cutoff_A: 6.0
  search_growth_factor: 2.0
  fallback_block_atoms: 4096
  contact_cutoff_A: 5.0
  heatmap_max_A: 20.0
  groups: true
```

| Parameter | Meaning |
| --- | --- |
| `selection` | Protein/query selection. |
| `lipid_resnames` | Used to build `resname ...` when `membrane_selection` is omitted. |
| `membrane_selection` | Explicit MDAnalysis membrane selection; overrides `lipid_resnames`. |
| `exclude_hydrogens` | Remove hydrogens from query and membrane selections. |
| `expand_residue_atoms` | Expand selected residues to all residue atoms before distance calculation. |
| `initial_search_cutoff_A` | Initial capped-neighbor search radius. |
| `search_growth_factor` | Multiplier used only for residues that still lack a membrane neighbor. Must exceed 1. |
| `fallback_block_atoms` | Membrane atom block size for the exact full-distance fallback. Controls memory use. |
| `contact_cutoff_A` | Optional cutoff used to report fraction of frames in membrane contact. |
| `heatmap_max_A` | Optional display ceiling for distance heatmaps. |
| `groups` | Also write group-filtered residue profiles. |

The fast path uses `capped_distance` rather than constructing the full proteinxmembrane matrix. Only unresolved residues expand the search radius. A blocked exact distance calculation is used as a fallback, so the final minimum remains exact. The framexresidue output is written as a `float32` memory-mapped array to avoid holding the full matrix in RAM.

This analysis also supports named `runs:` in the same way as `density_1d`, which is useful for comparing different membrane selections or protein regions.

## DSSP

```yaml
dssp:
  enabled: true
  selection: "protein"
  guess_hydrogens: true
  backend: serial
```

| Parameter | Meaning |
| --- | --- |
| `selection` | Structure used for DSSP. Defaults to `selections.protein`. |
| `guess_hydrogens` | Let MDAnalysis infer backbone hydrogens when needed. |
| `backend` | MDAnalysis DSSP backend. `serial` is the simplest/default mode. |
| `n_workers`, `n_parts` | Optional MDAnalysis backend parameters when a non-serial backend is explicitly selected. |
| `heavyatom_names`, `hydrogen_name` | Advanced overrides for DSSP atom naming. |

Raw per-frame residue assignments are retained so secondary-structure occupancies can be block bootstrapped.

## SASA

```yaml
sasa:
  enabled: true
  selection: "protein"
  algorithm: shrake_rupley
  probe_radius_A: null
  n_points: null
  n_slices: null
  n_threads: null
```

| Parameter | Meaning |
| --- | --- |
| `selection` | Atoms passed to FreeSASA. Defaults to `selections.protein`. |
| `algorithm` | `shrake_rupley` or `lee_richards`. |
| `probe_radius_A` | Optional FreeSASA probe-radius override. |
| `n_points` | Optional Shrake-Rupley point count. |
| `n_slices` | Optional Lee-Richards slice count. |
| `n_threads` | Optional FreeSASA thread count. |

Per-atom areas are summed to residues with `numpy.bincount`, avoiding a Python residue loop for aggregation. Stored SASA values are in nm^2.

## Salt bridges

```yaml
salt_bridges:
  enabled: true
  cutoff_A: 4.0
```

Detects acidic ASP/GLU oxygen to basic LYS/ARG nitrogen pairs within the cutoff. The detector uses `capped_distance`, so only contacts inside the cutoff are materialized rather than a full pairwise matrix.

## Hydrogen bonds

```yaml
hydrogen_bonds:
  enabled: true
  mode: protein_protein
  d_h_cutoff_A: 1.2
  d_a_cutoff_A: 3.0
  angle_cutoff_deg: 150.0
  update_selections: false
```

Uses `MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis`.

| Parameter | Meaning |
| --- | --- |
| `mode` | Interaction scope used by the configured detector (`protein_protein` by default). |
| `d_h_cutoff_A` | Donor-hydrogen cutoff. |
| `d_a_cutoff_A` | Donor-acceptor cutoff. |
| `angle_cutoff_deg` | Donor-hydrogen-acceptor angular cutoff. |
| `update_selections` | Whether MDAnalysis should update dynamic selections during the trajectory. |

Atom/residue labels are precomputed once and NumPy-indexed into detected events rather than rebuilt for every H-bond row.

## Ligand contacts

```yaml
ligand_contacts:
  enabled: false
  cutoff_A: 4.0
  save_atom_events: true
```

Requires `selections.ligand`. Outputs minimum ligand-protein distance, residue-contact count, and interaction event tables at ligand-residue / ligand-atom-residue / ligand-atom-atom resolution where available.

## Interaction lifetimes

```yaml
lifetimes:
  enabled: true
  gap_factor: 1.5
```

`gap_factor x sampled_dt` is the largest time gap considered continuous within one interaction event. Events touching the first or final analyzed frame are marked left/right censored. Mean and median lifetime statistics use completed events; occupancy still uses all observed frames.

---

# Statistics

```yaml
statistics:
  enabled: true
  bootstrap_draws: 5000
  ci: 0.95
  seed: 2026

  block:
    method: moving
    length: auto
    tau_multiplier: 2.0
    max_acf_lag: null

  hierarchy:
    block_within_replicate: true
    resample_replicates: true

  hypothesis_tests:
    enabled: true
    design: independent
    method: permutation
    permutation_draws: 10000
    exact_max_assignments: 100000
    multiple_testing: fdr_bh
```

## Explanation

MD trajectory frames are correlated. Treating every frame as an independent observation would make confidence intervals much too narrow. MDscope therefore separates **temporal uncertainty inside each replicate** from **replicate-to-replicate uncertainty**.

### 1. Estimate the correlation time

For each scalar time series (or an appropriate scalar summary signal for a profile/interaction), MDscope computes the autocorrelation function using an FFT. FFT evaluation is used because it scales approximately as `N log N` instead of evaluating every lag directly.

The integrated autocorrelation time is estimated from the initial positive part of the autocorrelation curve:

```text
tau_int ~= 0.5 + sum(ACF at positive lags)
```

The table reports:

- `tau(fr)`: integrated autocorrelation time in analyzed frames.
- `tau(ps)`: the same quantity in ps.
- `N_eff`: approximate effective number of independent observations, `N / (2*tau)`.

### 2. Choose a block length

With `length: auto`, the default is approximately:

```text
block_length = ceil(tau_multiplier * tau_int)
```

With the default `tau_multiplier: 2.0`, blocks are roughly twice the estimated correlation time.

`blocks/pseudo` in the bootstrap tables is:

```text
ceil(number_of_analyzed_frames / block_length)
```

It tells you how many contiguous blocks are sampled to build one bootstrap pseudo-study for that replicate.

### 3. Bootstrap inside each replicate

A moving-block bootstrap samples contiguous trajectory blocks rather than isolated frames. This preserves short-range serial correlation while generating a distribution for the replicate statistic.

Examples:

- RMSD/Rg/time-series metrics: recompute mean/median from the sampled blocks.
- SASA/DSSP/distance/density profiles: recompute the profile across sampled frames.
- RMSF: recompute RMSF from sampled aligned-coordinate frames.
- Interaction occupancy: bootstrap the binary present/absent trajectory.

### 4. Bootstrap across replicates

For each system-level bootstrap draw, MDscope resamples whole replicates and then selects one temporal-bootstrap realization from each selected replicate. This is the second hierarchical level:

```text
correlated frames -> temporal block bootstrap -> replicate estimates
replicates        -> replicate bootstrap      -> system estimate
```

The system confidence interval therefore reflects both within-trajectory temporal uncertainty and between-replicate variability.

### 5. System differences

When two systems are compared, matching system-bootstrap draws are subtracted:

```text
delta = system_2 - system_1
```

Outputs include the observed/bootstrapped difference, CI, and bootstrap direction probabilities such as `P(delta > 0)`. These direction probabilities are useful descriptive uncertainty measures but are **not** labeled as p-values.

### 6. P-values

P-values come from replicate-level randomization, not frame-level permutations.

- `design: independent`: replicate labels are permuted between systems.
- `design: paired`: matching replicate IDs are treated as pairs and sign-flipped. Use this only when the simulation design genuinely created matched pairs.
- Exact enumeration is used when the number of assignments is below `exact_max_assignments`; otherwise `permutation_draws` random permutations are used.

With only three independent replicates per system, exact permutation p-values are necessarily coarse because there are only `C(6,3)=20` possible group assignments. A non-significant p-value in that setting does not mean the bootstrap CI/effect size is uninformative; report the effect, CI, replicate behavior, and p/q-value together.

### 7. Multiple testing

`multiple_testing: fdr_bh` applies Benjamini-Hochberg FDR correction within logical analysis/comparison families. Both raw `p_value` and adjusted `q_value` are written where applicable.

## Statistics parameters

| Parameter | Meaning |
| --- | --- |
| `bootstrap_draws` | Number of temporal/hierarchical bootstrap realizations. Increase for final analyses if needed. |
| `ci` | Confidence level, e.g. `0.95`. |
| `seed` | Reproducible random-number seed. |
| `block.method` | Currently moving-block bootstrap. |
| `block.length` | `auto` or an explicit analyzed-frame block length. |
| `block.tau_multiplier` | Multiplier converting estimated tau to automatic block length. |
| `block.max_acf_lag` | Optional maximum autocorrelation lag in frames. `null` chooses a bounded trajectory-dependent value. |
| `hierarchy.block_within_replicate` | If false, temporal block length becomes 1. Normally leave true for MD. |
| `hierarchy.resample_replicates` | If true, system CIs include replicate-level resampling. |
| `hypothesis_tests.enabled` | Turn replicate-level tests on/off. |
| `hypothesis_tests.design` | `independent` or `paired`. |
| `hypothesis_tests.method` | Currently permutation/randomization testing. |
| `hypothesis_tests.permutation_draws` | Monte-Carlo permutation count when exact enumeration is too large. |
| `hypothesis_tests.exact_max_assignments` | Maximum assignment count for exact enumeration. |
| `hypothesis_tests.multiple_testing` | `fdr_bh` for Benjamini-Hochberg correction. |

---

# Cache

```yaml
cache:
  aligned_coordinates: true
  dtype: float32
  keep: true
```

| Parameter | Meaning |
| --- | --- |
| `aligned_coordinates` | Store aligned feature coordinates in compressed chunks. Required for temporal block-bootstrap RMSF. |
| `dtype` | Coordinate-cache dtype. `float32` is recommended; it roughly halves storage versus float64 and is more than adequate for typical RMSF uncertainty calculations. |
| `keep` | Retain caches after the run. Currently the cache is not automatically deleted when false; treat this as a reserved cleanup preference. |

---

# Plots

```yaml
plots:
  enabled: true
  time_series: true
  distributions: true
  autocorrelation: true
  bootstrap_tables: true
  profiles: true
  interactions: true
  min_occupancy: 0.05
  replicate_lightness_range: [0.28, 0.72]
  font_sizes:
    base: 13
    title: 16
    suptitle: 18
    legend: 10
    ticks: 11
  dpi: 250
```

| Parameter | Meaning |
| --- | --- |
| `enabled` | Master plotting switch. |
| `time_series` | Time-series figures. |
| `distributions` | Distribution/histogram/violin figures where implemented. |
| `autocorrelation` | ACF figures. |
| `bootstrap_tables` | Render statistics as Matplotlib `plt.table` figures. |
| `profiles` | Profile figures and profile statistics tables. |
| `interactions` | Interaction occupancy/lifetime figures and tables. |
| `min_occupancy` | Hide extremely rare interactions from occupancy heatmaps below this fraction. Raw tables are still retained. |
| `replicate_lightness_range` | Lightness range used to distinguish replicate lines within a system color. |
| `font_sizes` | Plot typography. |
| `dpi` | Raster figure resolution. |

Bootstrap table figures include replicate-level `tau(fr)`, `tau(ps)`, effective N, block length, and `blocks/pseudo`. System rows summarize replicate-specific values rather than pretending one temporal correlation scale applies to every trajectory.

---

# Main outputs

A typical output tree is:

```text
mdscope_results/
+-- resolved_config.yaml
+-- run_metadata.json
+-- .state/
+-- SYSTEM/
|   \-- REP/
|       +-- raw/
|       +-- cache/
|       +-- summary/
|       +-- profiles/
|       \-- interactions/
\-- project_summary/
    +-- statistics/
    |   +-- time_series/
    |   +-- profiles/
    |   \-- interactions/
    \-- figures/
        +-- time_series/
        +-- profiles/
        \-- interactions/
```

Important numeric outputs are stored primarily as Parquet tables, with large framexresidue matrices stored as NumPy arrays/memmaps when that is substantially more efficient.

---

# Current limitations

1. **No built-in scheduler submission.** `run` is serial; external arrays can call `replicate --index N`, then `finalize` must be run once after all replicate jobs succeed. MDscope does not submit or monitor scheduler jobs.
2. **Small replicate counts limit p-value resolution.** With 3 vs 3 independent replicates, exact permutation tests have only 20 assignments; p-values are therefore coarse. Hierarchical CIs and effect sizes should be interpreted alongside p/q-values.
3. **Block bootstrap assumes approximate stationarity.** Strong uncorrected drift or state transitions can make one global autocorrelation time/block length a poor summary. Inspect time series and ACF plots before relying on the CI.
4. **Automatic tau uses an initial-positive-sequence estimate.** This is simple and robust for routine use, but it is not a full convergence diagnostic and may be unstable for short/noisy trajectories.
5. **Interaction lifetime inference is not survival analysis.** Left/right-censored events are identified and excluded from completed-event mean/median lifetimes, but Kaplan-Meier or other censored survival models are not implemented.
6. **Leaflet assignment is fixed from the initial analyzed structure.** Lipid flip-flop or major membrane reorganization is not dynamically reclassified.
7. **Membrane thickness is global.** It is based on upper/lower headgroup COM separation, not a local thickness surface/map.
8. **1D density bins are established from the initial box length.** Per-frame density normalization uses the current box volume, but very large changes in the profiled box dimension can make fixed bin edges undesirable.
9. **Residue-to-membrane distances can be storage-heavy.** The `float32` memmap reduces RAM pressure, but the full framexresidue matrix is intentionally retained for hierarchical profile statistics and heatmaps.
10. **RMSF hierarchical CIs require the aligned-coordinate cache.** Turning that cache off saves disk space but removes frame-level RMSF resampling information.
11. **SASA requires FreeSASA and is an additional trajectory pass.** It is not part of the shared streaming pass because the external surface calculation dominates its cost.
12. **Hydrogen-bond quality depends on topology/atom naming.** Missing/incorrect hydrogens or unusual donor/acceptor naming can affect MDAnalysis H-bond detection.
13. **Salt-bridge chemistry is intentionally simple.** The built-in detector covers common ASP/GLU acidic atoms and LYS/ARG basic atoms at a distance cutoff; protonation-state-specific or alternative definitions require code/config extension.
14. **No automatic equilibration detection.** `trajectory.start` must be chosen by the user if equilibration should be excluded.
15. **No automatic cache cleanup yet.** `cache.keep` is retained as a preference, but automatic deletion is not currently implemented.
16. **System topology is shared within a system.** All replicates under one system use the same configured topology path.
17. **The pipeline does not replace convergence assessment.** Replicate agreement, trajectory length, state sampling, and physical interpretation still need to be checked explicitly.

---

# Practical recommendations

For production MD comparisons:

- Use at least three independent replicates when possible, and more when feasible.
- Exclude equilibration explicitly with `trajectory.start` rather than relying on the statistics to absorb drift.
- Keep `block_within_replicate: true` and inspect the ACF/tau outputs.
- Keep `resample_replicates: true` for system-level confidence intervals.
- Use `design: paired` only when replicate IDs are genuinely paired by the simulation design.
- Keep `aligned_coordinates: true` when RMSF confidence intervals matter.
- Treat effect size + hierarchical CI + replicate behavior as primary scientific evidence; use p/q-values as complementary inferential summaries.


## Modular architecture and diffusion fits (experimental branch)

RMSD and RMSF now live in separate modules (`analysis/rmsd.py` and `analysis/rmsf.py`). Other structural and interaction analyses are similarly split into single-purpose modules while `analysis/streaming.py` reads the XTC once for inexpensive framewise measurements. `analysis/registry.py` is the extension point for restartable specialized analyses. The existing statistics and plotting backends are **shared by output data structure**, not duplicated in individual analysis scripts. See [NEURONAL_EXTENSION.md](NEURONAL_EXTENSION.md) for a module map and contribution contract.

The existing `membrane.msd` code is retained in `analysis/msd.py` and now optionally produces **replicate-level diffusion coefficients** as well. Configure `analyses.membrane.msd.slope_start_fraction` (default `0.2`) and `slope_end_fraction` (default `0.8`). This is a descriptive **single-origin MSD** slope fit; it does not claim species-resolved or multi-origin diffusion, and it does not treat time-correlated frames as independent replicates.

Prior results with a mismatched code/configuration hash are archived in `.archives/` instead of silently deleted. A significant refactor may still trigger expensive reanalysis: choose a new project.output path if you do not want the old results invalidated.


## Lipid analyses from issues #1–#5

The optional membrane extensions include leaflet-resolved C–H order profiles, terminal/Slurm progress indicators, local bilayer-normal deviation, protein-centered lipid enrichment/depletion, and lipid–lipid hydrogen-bond events grouped by lipid type. See [ISSUE_ANALYSES.md](ISSUE_ANALYSES.md) for YAML examples, precise definitions, saved tables and limitations. They use the existing MDscope replicate-aware statistical engine and plot conventions.


### Monitoring progress in Slurm

In a noninteractive Slurm job, MDscope writes newline-delimited progress to
**stderr** (the `#SBATCH --error` file), not a terminal `tqdm` bar.
It reports every **1% or 25 completed frames**, whichever comes first,
including elapsed time, average frames/second and estimated time remaining
for the current trajectory stage. The estimate may fluctuate while the
analysis warms up and does not include final plotting/statistics.

```bash
tail -f slurm-out/mdscope/mdscope_synd4-<jobid>_<taskid>.err
```

Optional environment overrides in the Slurm script (positive integers):

```bash
export MDSCOPE_PROGRESS_EVERY_FRAMES=25
export MDSCOPE_PROGRESS_EVERY_PERCENT=1
```

These settings do not alter which frames are analyzed or how the scientific
observables are computed. Already-running jobs continue using the Python
code loaded at job startup; submit a new job to use an updated checkout.


### SLURM progress and source checkout

`run_mdscope.sh` runs the checked-out `MDscope.analysis.cli`, not the separately
installed lowercase `mdscope` package. Submit from this repository after
`git pull --ff-only origin main`. The default configuration is
`synd4_full_mechanisms.yaml`; provide a YAML path as the first `sbatch`
script argument to override it. Before submitting, create the SLURM log
directory with `mkdir -p slurm-out/mdscope` (SLURM opens log files before
executing the script). The script combines stderr into the `.out` file so
batch progress and exceptions are visible together. It prints imported source
paths and the checkout commit for reproducibility. Batch progress reports
every 1% by default; `MDSCOPE_PROGRESS_EVERY_PERCENT` overrides the interval.

### Membrane MSD coordinates and drift

Membrane MSD is computed directly from **raw coordinates** in the one-pass
trajectory reader, before fitting/alignment for protein RMSD/RMSF. The aligned
feature-coordinate cache is **not** used for MSD. Supply the original
unrotated, non-protein-centered MD trajectory and its compatible topology.
Protein centering performed **before writing the input trajectory** can change
lab-frame displacement and cannot be reliably reconstructed by MDscope.
Minimum-image frame-to-frame steps account for periodic wrapping, provided
per-sampled-frame displacement stays below half the box length. The membrane
MSD is a single-origin descriptive estimator, not a multi-origin estimator.

Optional `analyses.membrane.msd.remove_drift: true` subtracts the instantaneous
mean displacement of selected lipid residue centers in each step (the
lipid-frame reference). This reduces common translation artifacts but changes
the physical reference frame and introduces a small finite-population bias;
leave it disabled when laboratory-frame MSD is desired. The option is included
in the configuration hash and is restart-safe. For scientific interpretation,
check membrane drift, trajectory preprocessing, and the fit regime.

### Performance profiling

Set `MDSCOPE_PROFILE=1` to write `streaming_phase_profile.json` in each
replicate output directory. The profile separates trajectory-frame decoding,
fused membrane observers, extensions, distance observers, membrane metrics,
alignment, density, and contact updates. Each family also writes its own
more detailed timing JSON where supported. The shipped SLURM script enables
this low-overhead profiling by default; set `MDSCOPE_PROFILE=0` to disable.
Compare profiles across replicates by **seconds per analyzed frame**, while
also recording CPU allocation, system atom count, storage and sampled frames.
A GPU allocation does not accelerate ordinary NumPy/MDAnalysis analyses.
Synthetic microbenchmarks can be run with
`pytest -q -s tests/test_performance_benchmarks.py`; these check numerical
equivalence and print measured timings without brittle fixed timing thresholds.
