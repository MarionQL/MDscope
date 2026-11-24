# Configuration

Use [`synd4_full_mechanisms.yaml`](../synd4_full_mechanisms.yaml) for the complete, maintained SDC4 configuration. The examples below show **patterns**, not validated selections for another topology.

## Config paths, output location and relocation

- Relative `project.output`, topology, trajectory and `msd_trajectory` paths resolve from the **YAML directory**, not the shell directory. Absolute paths stay unchanged.
- `mdscope validate` and `mdscope plan` show the resolved output; JSON plans expose `resolved_output_directory` without reading trajectories.
- **Existing runs:** keep the original resolved inputs and raw `project.output` string. Path-string edits can change the saved configuration hash; incompatible resume fails.
- Invoke the original YAML by absolute path when working elsewhere. A symlink retains the original resolved parent where supported.
- An absolute output path survives YAML relocation only if relative inputs still resolve correctly; edits to those strings may still invalidate the hash.
- Move a YAML with its relative directory tree only for a **new project in fresh output**. Changing `project.output` is not a saved-product migration; no automated migration is provided.

## Frames, time and coordinate cache

```yaml
trajectory: {start: null, stop: null, stride: 1, chunk_frames: 1000}
time: {dt_ps: null, origin_ps: 0.0}
cache: {aligned_coordinates: true, dtype: float32}
```

| Key | Meaning |
|---|---|
| `trajectory.start` / `stop` / `stride` | Frame-index window, Python-slice semantics (`stop` exclusive, `null` = trajectory end); `stop` must exceed `start`, `stride` ≥ 1. Choose the production window explicitly: equilibration is never trimmed automatically |
| `trajectory.chunk_frames` | Frames per committed chunk (≥ 1). Sets the restart granularity and peak memory, not the science; solvation `shell_stability` is the one diagnostic that treats chunks as time windows |
| `time.dt_ps`, `time.origin_ps` | When `dt_ps` is set, time = `origin_ps + i · dt_ps` for the *i*-th **analyzed** frame (after `stride`); otherwise the trajectory's own times are used and a trajectory without usable times fails |
| `cache.aligned_coordinates`, `cache.dtype` | Save fitted feature coordinates (`float32` by default) for RMSF and reuse. Enabling it requires the primary streaming pass |

## Analysis options at a glance

| Feature | YAML key | Output / rule |
|---|---|---|
| Protein groups | `selections.groups` | Shared frame pass; `summary/group_timeseries.parquet` |
| Multiple density/distance/SASA instances | `analyses.<name>.runs` | Each named run has distinct selection, ID and files |
| Optional figures | `plots.enabled`, `plots.families` | Offline `mdscope plotting`; no trajectory reread |
| Statistics families | `statistics.families` | Override specific estimators and CI settings |
| Interaction screening | `statistics.interaction_min_occupancy` | Filter *statistical identities*, not saved events |
| Matched analysis window | `statistics.matched_window_ps` | Additional matched-time statistics |
| Automatic HS identification | `analyses.sdc4_mechanisms.hs_auto` | Requires covalent bond topology |
| Glycan chains for hydration | `selections.glycans` + `analyses.hydration.glycans` | Requires covalent bond topology; [below](#shared-glycan-definitions) |

## YAML contract for generalization

- Expose new user-selectable science in the shared validated YAML and CLI; Studio uses the same configuration. Internal extraction alone needs no new option.
- Declare defaults, schema, units, ranges and semantic validation. Preserve omitted versus explicit-empty selections; reject unsupported options.
- Configure selections, motif/chain scope, lipid naming, glycan attachments and supported physical references through implemented fields. Named runs work only for analyses that support them.
- Preserve SDC4 CLI/YAML/output contracts, configured motifs, domains, HS identities, lipid anchors/classes and reference models.
- Resolve presets/overrides once. Expose effective settings, prerequisites and provenance; preserve Studio import/export. Groups/presets must not enable unrelated work.
- Verify non-SDC4 behavior and applicable SDC4 numerical references. Keep all three `synd4_*.yaml` configurations loading and round-tripping; document necessary migrations and equivalent requested work.
- Examples must use implemented keys. Schema validity and naming presets do not validate chemistry, geometry or numerical parity.

Existing controls include groups, analysis-specific selections, membrane-order tail mappings, SDC4 domain/HS/motif/lipid definitions, metric selection and exact figure IDs. These are not universal keys. Shared motifs/states and residue-response controls are documented below.

### Groups vs named runs

```yaml
selections:
  groups:
    chain_a: {selection: "segid PROA", alignment_selection: "segid PROA and name CA"}
analyses:
  rmsd: {enabled: true, groups: [chain_a]}
  rmsf: {enabled: true}
  density_1d:
    enabled: true
    axis: z
    bin_width_A: 1.0
    selections: {water: "resname TIP3 and name OH2"}
    runs:
      membrane_relative:
        membrane_center_selection: "resname POPC DPPC POPE POPI POPS PSM DOPC and name P"
        half_width_A: 100.0
```

- RMSD/Rg `groups` chooses a subset of shared groups; omitted means all. Global metrics stay in `summary/timeseries.parquet`.
- Density, residue–membrane distance and SASA support named `runs`. Do **not** add `runs` to RMSD/Rg.
- Membrane-relative density uses the midpoint of circular leaflet centers. Use **one comparable anchor per lipid**, a half-width divisible by the bin width and no `center_selection`. The half-width must fit within half the current normal box length. A different reference needs a new run ID.

## Lipid-tail naming for -S_CH

`analyses.membrane_order` owns the acyl-tail naming; SDC4 `local_order` uses the
same `tails`/`lipid_overrides` with its own lipid selection.

| Assumption | Kind |
|---|---|
| Tail atom names, residue names, `exclude_resnames` (`CHOL`, `CHL1`), headgroup anchor names (`P`, `O3`, `NF`) | Naming only: configure per topology |
| Explicit, bonded C–H hydrogens; tail carbons form one bonded chain numbered from the carbonyl | Scientific constraint: checked, failures are errors |
| One topology residue per lipid; fixed initial leaflets; global axis as bilayer normal | Scientific model: reported, not converted |

- **Corrected default naming (issue #272):** `sn1: ^C3(\d+)$` maps the
  CHARMM36 glycerophospholipid **sn-1** chain (C31…), and `sn2: ^C2(\d+)$`
  maps **sn-2** (C21…). The matching `tail_preset: charmm36` is optional,
  never inferred from residue names; the two maintained full SDC4 YAMLs select
  it explicitly. These mappings were checked against bonded CHARMM36 POPC
  fixtures and the SDC4 release TPRs.
- **Historical results generated before #272:** saved `lipid_order__sn1*`
  and SDC4 local-order `sn1` products measured the **C2x/sn-2** chain;
  saved `sn2` products measured **C3x/sn-1**. Their numeric order parameters
  are not corrected by renaming files. Older statistics, figures and findings
  inherit those mislabeled identities. Archive them as legacy and **rerun
  replicate analysis and statistical finalization/plotting in a new output
  directory**. Do not resume or mix files from the prior output directory:
  changing the resolved default/preset changes the configuration hash.
- The opt-in `charmm36` preset also maps CHARMM36 sphingolipids (PSM, CER*)
  using `sphingosine: ^C(\d+)S$` and `n_acyl: ^C(\d+)F$` as
  separate tail identities. Sphingoid and N-acyl carbons are never pooled
  into glycerophospholipid `sn1`/`sn2` profiles. Pooled scalar observables
  may change because they include those additional sphingolipid tails.
- Precedence: `lipid_overrides.<RESNAME>.<tail>` naming > `tails.<tail>` >
  preset tail > default. With a preset, an explicit `tails` entry replaces that
  tail; without one it merges into the default tail (existing behavior).
  `carbon_pattern` needs an integer capture group (the carbon position);
  `carbon_names` (optional matching `positions`) wins over a merged pattern.
- `mdscope validate` prints, per species and tail, mapped lipids, positions,
  C–H bonds and exclusions (no matching carbons, fewer than
  `min_carbons_per_tail`, no C–H). It fails when bonds contradict the mapping
  (bonded carbons of different tails, non-consecutive positions, split chains,
  non-carbon elements, overlapping tails), when no explicit hydrogens exist
  (united-atom), or when `leaflet_anchor_selection` misses an order lipid.
  Without topology bonds, mappings are reported as unverified.
- Multi-residue lipids (for example AMBER Lipid21 head/tail residues) are
  analyzed per residue; validate warns when one bonded molecule holds several
  order residues.

User-defined example (not a preset): Maciejewski–Rog POPC with sequential atom
numbers, names from the NMRlipids `mappingPOPCMaciejewskiRog1.yaml`; verify
them against your topology.

```yaml
analyses:
  membrane_order:
    enabled: true
    lipid_selection: "resname POPC"
    leaflet_anchor_selection: "resname POPC and name P8"
    tails:
      sn1:
        carbon_names: [C36, C38, C39, C40, C41, C42, C43, C44, C45, C46, C47, C48, C49, C50, C51, C52]
        positions: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
      sn2:
        carbon_names: [C15, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26, C27, C28, C29, C30, C31, C32, C33]
        positions: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
```

## Shared motifs and saved-feature states

A shared definition discovers topology membership; it never enables an analysis
or assigns biological function. Solvation is the first generic measurement
consumer. Existing inline solvation sequences and SDC4 motif settings keep their
identities and numerical semantics. Other consumers must explicitly adopt the
shared resolver before accepting motif references.

```yaml
selections:
  motifs:
    acidic_patch:
      sequence: DE                       # or pattern: "D.E" (Python regex)
      selection: protein                 # residue domain; CA atoms identify sequence
      chains: null                       # all chains; [] selects none
      overlap: false
      missing: error                     # error or allow
analyses:
  solvation:
    enabled: true
    water_selection: "resname SOL and name OW"
    targets:
      motif:
        patch: {motif: acidic_patch, atom_selection: "not name H*"}
  shared_states:
    enabled: true
    preset: null                         # sdc4 is the backward-compatible default
    features: [patch_water]
    feature_sources:
      patch_water:
        source: solvation
        group: solvation/water/motif/patch
        metric: mean_n_solvent
        units: molecules
    alignment: exact
    scaling: standardize                  # median/std, the supported scaling rule
    n_states: 2
```

- Merge the fragment with project/systems settings; disable unrelated analyses. Match water names and sequences to the topology.
- Preview with `mdscope validate` or Studio preflight. `mdscope plan` exposes motif definitions, state inputs and source prerequisites.
- Hydration/bridge/lipid identities remain templates until saved readers verify them. Preflight reports ordered residue indices/numbers/names, segments and chain labels.

### Motif identity and consumers

- Membership uses topology residue indices within contiguous selected fragments/segments. Unknown/modified residues cannot match wildcards.
- `sequence`: uppercase standard amino acids. `pattern`: consuming Python regex. Chain labels use chainID, then segid.
- Omitted chains means all; `chains: []` with `missing: allow` previews unavailable. `selections.motifs` must be a mapping; `{}` defines none.
- Solvation requires every requested motif target to match. Absence is not zero hydration. `overlap: true` permits discovery previews, but solvation rejects shared target atoms.
- Solvation counts each instance/frame, then averages instances per solvent/target-class/site-label in `species_timeseries.parquet`. Use separate chain-scoped labels for chain-specific state features.
- `local_environment.motifs` accepts regexes within `protein_selection` (nonoverlapping, zero matches allowed) or `{motif: name}` references preserving the shared definition's scope/policies.
- Local-environment defaults: `nxip: N.IP`, `gxxxg: G...G`, `hs_candidate: SG`; `{}` measures none. Without chainID/segid labels, set `{}`. Group names remain `motif_<label>_<segid>_<first resid>`.
- References preserve aggregation, RDF normalization, water identity and heavy-atom selection. Solvation saves membership in `metadata.json` and target tables.

| State source | Saved scalar identity (`group`, `metric`) |
|---|---|
| `summary` | `null` or named group; numerical column such as `rmsd_A` |
| `solvation` | `solvation/<solvent>/<target_class>/<site_class>`; `mean_n_solvent` or `hydrated_fraction` |
| `hydration` | `hydration/<compartment>/<level>/<label>`; `unique_water_count` or `unique_water_per_site` |
| `solvation_bridges` | Existing bridge-pair group identity from the saved scalar catalog; supported numerical bridge column |
| `lipid_contacts` | `lipid_contacts/<leaflet>/<species_a>/<species_b>`; `neighbors_per_molecule` or `enrichment` |

### State fitting and outputs

- Generic mode (`preset: null`) requires nonempty features exactly matching `feature_sources` aliases; it does not enable SDC4. Each source declares product, exact group/metric and units. Units are declarations, not conversions or chemical validation.
- Readers match saved time-series statistics. Features within a replica require the **same finite, unique time grid**; no interpolation, nearest-time or row-position joins. Replica lengths/windows may differ.
- Missing values retain observation positions; fit complete vectors only. `min_complete_fraction: 0.90` rejects excessive loss. Irregular/gapped data disable temporal block intervals; transitions never bridge dropped rows.
- One median/std scaler and seeded k-means model train on equal frame counts per replica (at most `training_frames_per_replicate: 2000`). Correlated features can overweight one property despite standardization.
- Constant/all-missing features, missing products, mismatched times and insufficient distinct vectors fail. Disabled analysis fits nothing.
- SDC4 preset: omitted/null `features` uses six legacy columns; `features: []` fails.
- `project_summary/shared_states/model.json` records sources, units, preset, scaling, complete windows and training counts. Successful reruns archive old results under `project_summary/archive/shared_states*`; preparation failures preserve them.
- SDC4 conditioning products require the SDC4 preset. Clusters do not establish metastability, causation or convergence.

## Shared glycan definitions

A definition discovers covalently attached glycan chains from **topology bonds**;
like a motif, it enables nothing. `analyses.hydration` is the generic consumer.

```yaml
selections:
  glycans:
    o_glycans:
      glycan_selection: "resname BGALNA BGAL ANE5AC"   # every residue of each chain
      attachment_selection: "protein and resname THR SER and name OG1 OG"
      min_glycan_residues: 1      # default 1
      expected_chains: null       # default: no count check
analyses:
  hydration:
    enabled: true
    water_oxygen_selection: "resname TIP3 and name OH2"
    cutoff_A: 3.5
    glycans: [o_glycans]          # omitted or [] selects none
```

- Each bonded component of `glycan_selection` must bond to exactly **one**
  `attachment_selection` atom (the allowed linkage chemistry); a bond to any other
  atom outside the glycan selection, an unattached or undersized fragment, two
  attachment atoms/residues and an `expected_chains` mismatch all fail.
- Chains are labeled `<chain>_<residue letter><resid>` (`A_T12`, `PROA_S39`), using
  chainID then segid. Sugars are numbered from the attachment sugar outward along
  bonds; **branched chains fail** because they have no unique position.
- Hydration sites are each sugar's oxygen atoms (`O*` names). Each definition is a
  hydration compartment with `chain`, `sugar_type` (`<chain>:<resname>`) and
  `total` unions in `aggregate_counts.parquet`; per-sugar identities are in
  `target_metadata.parquet`, linkages in `glycan_attachments.parquet`.
- Names are letters, digits, `_`, `.` or `-`, excluding `protein`, `membrane` and `hs`.
  Definitions sharing atoms, or overlapping protein/membrane targets, fail.
- Residue names and bonds come from the force field (TPR/PSF); bare PDBs cannot
  define chains. Check `mdscope validate`, which lists each chain's linkage and
  sugars, before running. SDC4 keeps `sdc4_mechanisms.hs_auto` (compartment `hs`,
  `min_glycan_residues: 3`, `expected_chains: 6`).

## SASA environments

```yaml
analyses:
  sasa:
    enabled: true
    runs:
      isolated: {target_selection: "protein", occluder_selection: "protein", n_points: 200}
      membrane:
        target_selection: "protein"
        occluder_selection: "protein or resname POPC POPE DOPC DPPC PSM CHL1"
        reference_run: isolated
        periodic: true
        n_points: 200
```

| Contract | Rule |
|---|---|
| Target SASA | FreeSASA target atom areas with selected occluders; target must be included in occluders |
| Units and files | nm²; `profiles/sasa/<run_id>/residue_sasa.parquet`; burial in `residue_burial*.parquet` |
| Burial | Isolated minus environment SASA; fraction unavailable when isolated area is zero |
| Comparable runs | Same target identities, radii, probe, algorithm and sampling |
| Radii | Geometric element radii, **not** CHARMM Lennard-Jones radii; override `element_radii_A` when needed |
| Periodic mode | Compact **whole** target within half-box; minimum-image local occluders, **not** molecule reconstruction |
| Validation | Confirm element inference, target selection, radii, box and a native-library reference on new topologies |

Default geometric radii (Å): H 1.20; C 1.70; N 1.55; O 1.52; P/S 1.80; F 1.47; Cl 1.75; Br 1.85; I 1.98. Ambiguous atom-name inference must be resolved in the topology.

## Select individual figures

When `plots.include` is absent, the existing family-level plot switches are unchanged. To produce **only** specific supported figures, provide exact plot IDs:

```yaml
plots:
  enabled: true
  include:
    - time_series:rmsd_A:timeseries
    - time_series:rmsd_A__group_chain_A:distribution
    - solvation:shell_rdf
  options:
    time_series:rmsd_A:timeseries:
      dpi: 300
```

- **Dynamic IDs:** `time_series:<saved_metric_or_group_stem>:<timeseries|distribution|autocorrelation|bootstrap_table|comparison_table>`; `profiles:<saved_profile_stem>:<profile|distribution|system_difference|bootstrap_table|comparison_table>`; `interactions:<kind>:<occupancy_heatmap|lifetime_distribution>`, with an additional `:<metric>:<bootstrap_table|comparison_table>` for interaction statistics tables.
- **Distance heatmaps:** `profiles:residue_membrane_distance__<run_id>:<system_id>:distance_heatmap`.
- **Named static figures:** `profiles:residue_response:heatmap`, `solvation:shell_rdf`, `solvation:water_and_ions`, `solvation:contact_survival`, `diffusion:diffusion_coefficients`, `lipid_enrichment:lipid_enrichment_heatmap`, `lipid_enrichment:lipid_area_density_heatmap`.
- **Shared-state figures:** `shared_states:state_populations`, `shared_states:state_sampling_diagnostics`.
- **SDC4 condition figures:** `sdc4:<system_id>:<figure_stem>`; supported stems include `nxip_height_occlusion_joint`, `vvhplvpl_proline_ramachandran`, `lipid_density_maps`, `lipid_enrichment_maps`, and the five `lipid_region_*` charts. A single PNG is saved under `project_summary/figures/sdc4/<system_id>/`; each independent replica appears as an explicitly labelled panel or trace. The *numerical* Parquet/NPZ measurements remain per replica in `<system_id>/<replicate>/extensions/sdc4_mechanisms/`. Previous per-rep diagnostic PNGs from older runs are not automatically deleted; archive them when regenerating plots. No trajectories or statistical CIs are pooled to draw these condition diagnostics.
- **Lipid-contact figures:** `lipid_contacts:<figure_stem>` for each saved `lipid_contact_matrix_*`, `lipid_contact_difference_*`, `lipid_contact_class_matrix_*` or `lipid_contact_class_difference_*` figure; valid endings are `_pair_count_undirected`, `_neighbors_per_molecule` or `_enrichment`.

Per-figure `plots.options.<plot ID>.dpi` is supported for resolution overrides (integer 50–600); it does not change scientific results, axis definitions or statistical estimands. Unsupported per-plot keys fail configuration validation.

- IDs match exactly; no wildcards/substrings. Invalid formats and unknown static IDs fail validation.
- Saved inputs and family prerequisites are still required; enabling a figure does not compute missing statistics.
- Data-dependent metrics/runs/species/systems are verified against existing products, not syntax alone.
- Explicit solvation IDs need no extra family opt-in. `plots.include: []` disables all figures.

### Compact statistics-table PDFs

- `plots.table_reports_pdf: true` publishes one paginated PDF and full-precision CSV per analysis under `project_summary/figures/statistical_reports/<analysis>/statistics_tables.*`. Reports read saved Parquet statistics only; plot-only runs need no trajectories and perform no inference.
- Historical `*:bootstrap_table` / `*:comparison_table` IDs select exact source rows within the consolidated report. Selecting multiple IDs produces one report with all selected sources; it does not expand to sibling metrics. `bootstrap_tables: false` suppresses descriptive/bootstrap and diagnostic sections while keeping comparisons. `table_reports_pdf: false` (the unchanged default) retains legacy per-metric PNGs. Enable PDF mode for consolidated reports.
- Comparison tables join condition means/CIs and available replica support on exact saved identities. Separate sections retain condition/replica summaries and saved diagnostics. Integrated autocorrelation, effective N and block summaries are included; full lag curves remain in their original Parquet and dedicated autocorrelation plots rather than being copied into PDF tables. Unregistered metrics remain visible under `unmapped`; no top-N filter is applied. Lipid order has pointwise six-column tail/leaflet contrast-CI matrices, with N-acyl/sphingosine kept separate from sn1/sn2. Absent/unselected entries are N/A, not zero.
- Red contrast intervals strictly exclude zero (finite, ordered bounds); bold q values fall below `statistics.hypothesis_tests.alpha` (default `0.05`). These criteria are distinct. PDF values are rounded for display; classification and CSV use full precision. The interval's saved `ci_level` is included; no fixed confidence level is assumed.
- PDF, CSV and manifest publish together with rollback; prior reports and selected old table PNGs are archived after rendering succeeds. The manifest records every source ID and input digest. Inventory rejects changed/missing sources or damaged companions. Incomplete/invalid statistics cannot replace existing reports. Legacy statistics without a publication manifest retain the existing unverified policy.
- `plots.bootstrap_tables` controls bootstrap pages; comparison pages remain eligible independently. Exact includes and family switches still apply; saved statistics are required.
- After successful PDF publication, superseded table PNGs move to `.archives/figure_tables/`. Other figures and full-precision Parquet remain unchanged.
- Each PDF's `statistics_tables.manifest.json` records ordered page identities, count, completion and digest. Inventory includes PDFs and PNGs; `formats=('png',)` restricts it.
- Rerenders archive each superseded PDF/manifest together under `.archives/figure_tables/pdf_reports/<figure directory>/<UTC timestamp>/`.
- Empty/unselected/disabled sessions archive stale PDFs and leave a manifest with status, old identities and archive path. Repeated runs preserve that state without archiving twice; `include_rejected=True` exposes it. Inventory paths are project-relative POSIX paths.
- Publication failure preserves the previous PDF/metadata. Use one table-report writer per project; overlapping plotting sessions are unsupported.
- `plots.enabled: false` renders nothing; inventory marks its products disabled. PDFs are presentation, not numerical exports.

### Suppressing uninformative lipid-footprint traces

A lipid-footprint metric whose available sampled values are **all exactly
zero** across every configured condition and replicate produces no descriptive
time-series, distribution or autocorrelation PNG. This avoids empty-looking
flat-line figures (for example, an unobserved BGAL leaflet/region), including
after applying the trailing rolling window. A metric with any nonzero finite
observation is still plotted, and non-footprint metrics retain their zero-valued
traces. Missing values are not interpreted as zeros.

Previously generated PNGs for the suppressed metric are moved into
`<project_output>/.archives/zero_only_footprint_figures/`, not deleted.
Full-precision numerical statistics, Parquet outputs and bootstrap/comparison
table reports remain available regardless of figure suppression.

### Time-series rolling averages

Use a trailing window measured in saved `time_ps`, not a fixed frame count:

```yaml
plots:
  enabled: true
  time_series: true
  time_series_rolling_window_ns: 100
```

- `null` (default): raw traces. Positive finite ns: mean over `(t - window, t]`, using available initial samples.
- Smooth each replica independently; exclude nonfinite values. Entirely missing windows stay missing.
- Display only: distributions, ACFs, bootstrap/CIs and saved measurements are unchanged.

## Plot and statistics controls

```yaml
plots:
  enabled: true
  families: {time_series: false, interactions: false, solvation: true}
  interaction_exclude_kinds: [lipid_hydrogen_bonds]
statistics:
  families:
    time_series: {enabled: true, bootstrap_draws: 1000}
    profiles: false
    interactions: {enabled: true, bootstrap: false, hypothesis_tests: false}
    solvation_rdfs: true
  interaction_min_occupancy: 0.20
  hypothesis_families:
    primary: [sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction]
```

| Setting | Meaning |
|---|---|
| `families` | Per-family `enabled`, `bootstrap`, `bootstrap_draws` (≥100), `ci`, `hierarchy`, `hypothesis_tests` |
| `bootstrap: false` | Keep point estimates; intervals and uncomputed directional probabilities become NaN, draws = 0 |
| `hierarchy.enabled: false` | Retain replica intervals/description; omit system and contrast bootstrap intervals |
| `block_within_replicate: false` | Block length 1: **assumes** independent frames within each replica |
| `interaction_min_occupancy: 0.20` | Retain an interaction if occupancy ≥20% in **any** replica; keep zeros from other replicas |
| `interaction_include_kinds` / `exclude_kinds` | Filter before Parquet loading; exclusions win |
| `hypothesis_families` | Predefine BH correction families; unlisted metrics are exploratory |
| `matched_window_ps` | `auto` or `[start, stop]`; preserve a separate full-window series |

Project-wide statistics settings (each family can override `bootstrap_draws`, `ci` and `hierarchy`):

| Key | Default | Meaning |
|---|---|---|
| `statistics.enabled` | `true` | Run offline statistics at finalization |
| `bootstrap_draws` | 2000 | Bootstrap draws (≥ 100) |
| `ci` | 0.95 | Interval level, strictly between 0 and 1 |
| `seed` | 2026 | Seed for bootstrap and permutation draws; reproducible only with the same NumPy version and `MDSCOPE_STATS_WORKERS` |
| `block.method` / `length` | `moving` / `auto` | Moving-block bootstrap within each replica; `auto` derives the length from the autocorrelation time, or set an integer ≥ 1 |
| `block.tau_multiplier`, `block.max_acf_lag` | 2.0, `null` | Block length = multiplier × integrated autocorrelation time; optional ACF lag cap (≥ 1) |
| `hierarchy.block_within_replicate`, `hierarchy.resample_replicates` | `true`, `true` | Two-stage resampling: frames within replicas, then replicas |
| `hypothesis_tests.design` | `independent` | `paired` requires identical replicate labels across systems |
| `hypothesis_tests.method`, `permutation_draws`, `exact_max_assignments` | `permutation`, 10000, 100000 | Exact enumeration when the assignments fit under the limit, Monte Carlo permutation (≥ 100 draws) otherwise |
| `hypothesis_tests.multiple_testing` | `fdr_bh` | Benjamini–Hochberg within each `hypothesis_families` family, or `none` |

Optional `statistics.families.time_series.adaptive_hierarchy`:

```yaml
min_draws: 1000
max_draws: 8000
batch_size: 500
stable_checks: 2
relative_tolerance: 0.05
```

- It monitors **both** mean/median CI endpoints, with simultaneous rank bands. Require positive absolute or relative tolerance.
- Status and actual draws are saved; `not_converged` is not silently accepted.
- This measures Monte Carlo endpoint precision **conditional on the fixed replica bootstrap pools**; it does not create additional independent replicas.

**Statistical safeguards**

- Frames are correlated. Replicas—not frames, bins, ions or lipid molecules—are the condition-level units.
- Irregular/gapped/nonfinite traces keep finite descriptive means but lose unsupported ACF/temporal CIs; invalid replicas cannot be silently removed from comparisons.
- Profile positions align by labels; unsupported/missing positions remain unavailable, not zero.
- RDF pointwise means/SD remain descriptive, regardless of `bootstrap` settings.
- A 3+3 exact two-sided unpaired permutation test has minimum p = 0.10. CI exclusion is **not** FDR-adjusted significance.
- Rerunning statistics archives old family tables/manifests under `.archives/statistics`. The active manifest is `running` until complete.

## Selective metric inference

Use exact **saved** metric identities. Omitting `statistics.metric_selection`
preserves legacy inference; a family set to `[]` skips its statistics without
disabling raw analyses.

```yaml
statistics:
  metric_selection:
    time_series:
      - {group: PROA, metric: rmsd_A, ci: true, hypothesis_tests: false}
      - {metric: membrane_thickness_A, inference: false}
    interactions:
      - {kind: salt_bridges, metric: occupancy}
    profiles:
      - {profile: rmsf}
    solvation_rdfs:
      - {solvent: water, target_class: protein, site_class: acidic}
    diffusion:
      - {metric: DOPC}
    coupling:
      - {mechanism: membrane_coupling}
plots:
  include: ["time_series:rmsd_A:bootstrap_table"]
```

| Family | Identity | Selectable inference |
|---|---|---|
| `time_series` | `group`, `metric` | CI, hypothesis tests |
| `interactions` | `kind`, optional `interaction`, optional `metric` | CI, hypothesis tests |
| `profiles` | `profile` | CI, hypothesis tests |
| `solvation_rdfs` | `solvent`, `target_class`, `site_class` | Shell contrast |
| `diffusion` | `metric` | Replicate contrast |
| `coupling` | `mechanism` | Replicate contrast |
| `lipid_contacts` | `level`, `leaflet`, `species_a`, `species_b`, `metric` | CI |

- Default `inference: true`; `false` keeps descriptive statistics with null CIs/p/q values. The first three families allow separate `ci`/`hypothesis_tests` switches; unsupported combinations fail.
- Requested identities must exist in saved data; example IDs are illustrative. Omitted interaction identity selects every identity of that kind; omitted metric selects all five interaction statistics. Occupancy screening still applies.
- `mdscope plan` previews candidates, identity fields, estimator eligibility, conditional files and figure patterns without XTCs. `mdscope.planning.inspect_saved_statistics(root)` reads actual saved schemas/IDs with a bounded identity limit.
- Analyses, metric inference and figures are independent controls. Replica-level uncertainty and pairing remain unchanged.
- Selective Monte Carlo streams are identity-stable; changing the hypothesis family can change BH q-values. Reruns archive previous statistics.

## Preview the work without analyzing a trajectory

```bash
mdscope plan analysis.yaml                    # schema-only; no TPR/XTC required
mdscope plan analysis.yaml --format json      # machine-readable preview
mdscope plan analysis.yaml --check-topology   # real topology selections/counts
mdscope plan analysis.yaml --check-inputs     # require TPR and trajectory paths
```

- Plans list registered artifacts, named runs, extra passes, statistics and figure families. Paths are replica-relative; nothing runs or is created.
- Schema-only planning validates options/constraints without checking files, selections, chemistry, PBC or optional backends.
- Topology mode inspects topology and required backends, never XTC frames.
- The catalog is partial: some raw/cache outputs, exact schemas and plot IDs are unregistered. Listed families remain conditional on saved inputs; they do not guarantee results or CIs.
- Run `mdscope validate` before production.

## Validate before long runs

Unknown YAML option names, including nested analysis and statistics options, now fail with a dotted path and a suggested correction when available. This check also applies to disabled analyses. System, group, species, motif, solvation-site and named-run identifiers remain user-defined; their option fields are still validated.

```bash
mdscope validate analysis.yaml
mdscope run analysis.yaml
# Optional profiling:
pytest -q -s -m benchmark tests/test_performance_benchmarks.py
```

Preflight first confirms optional backends required by enabled analyses (FreeSASA for `sasa`: `python -m pip install 'mdscope[sasa]'`) and exits nonzero with that install hint before reading a topology. It then checks topology selections, group/alignment sizes, lipid-anchor ownership, SASA references, paths, run IDs, valid frame ranges and paired replica labels **before** replacing outputs. Overlapping protein groups are allowed; competing lipid species anchors are not.

- Expected YAML/options/path/backend/preflight/index errors and interrupted `init` prompts print `mdscope: error: ...` to stderr and exit 2. Unexpected faults retain tracebacks.
- `mdscope --debug <command>` or `MDSCOPE_DEBUG=1` also traces expected errors.
- `report STATISTICS_DIR` requires a complete manifest and contrast coverage.
- Use `statistics.block.*`; retired aliases, `analyses.autocorrelation` and `analyses.distribution_statistics` are unsupported.

[Scientific estimators](analyses.md) · [Workflow and restart](workflows.md) · [Scientific contracts](developer/scientific-contracts.md)

## Residue-response discovery

Cross-metric residue-response discovery is a separate exploratory saved-data
stage. Configure it independently of profile inference and figure selection:

```yaml
statistics:
  residue_response:
    enabled: true
    profiles: [sasa_mean__isolated, residue_occlusion__any]
    chains: [PROA, PROB]
    min_length: 4
    minimum_cosine: 0.8
```

- Omitted settings preserve the earlier behavior: enabled, all supported saved
  residue profiles/chains, minimum length 3 and cosine similarity 0.7.
- `profiles` and `chains` are exact identity lists. `null` includes all supported
  available identities; `[]` includes none. Unsupported profile families fail
  rather than acquiring a fabricated residue interpretation.
- `min_length` is an integer at least 2; `minimum_cosine` is finite in [0, 1].
  Changing these settings does not alter the underlying profile contrasts, CIs
  or physical units. Two complete metrics are still required to form a region.
- Disabling discovery or selecting no identities archives previous response
  matrices/regions and records why no current product was produced. The separate
  older single-profile `residue_regions.parquet` screen retains its own behavior.
- `residue_response_manifest.json` records effective settings, replica identities,
  source table names/row counts and completion status. Running/stale/disabled
  response products are not advertised as current heatmaps. Legacy outputs
  without this new manifest remain explicitly unverified.
- Planning exposes effective settings and prerequisites; Studio's Advanced
  statistics YAML edits the same fields and preserves import/export semantics.
  Your configured SDC4 motifs and predefined regions remain independent.

Regenerate from saved statistics using the saved manifest's settings, or an
explicit validated YAML override (no topology/trajectory is opened):

```bash
mdscope scan-residue-regions RESULTS/project_summary/statistics
mdscope scan-residue-regions RESULTS/project_summary/statistics --config analysis.yaml
```

These regions remain exploratory groupings of the same estimated effects.
Pointwise CI exclusions do not validate discovered region boundaries or provide
selection-adjusted region intervals.

## Raw within-replica residue correlations

```yaml
statistics:
  residue_correlations:
    enabled: true
    min_frames: 30
    pairs: null  # curated hydration/packing/helix/SASA association pairs
```

This optional saved-data stage reads exact matching `time_ps` and
`topology_residue_index` records from SDC4 residue occlusion, whole-protomer
contact, DSSP and configured isolated/environment SASA runs. No trajectory pass
or bootstrap is performed. Output:
`project_summary/statistics/residue_correlations/replicate_correlations.parquet`.

- Both raw Pearson `r` and rank-based Spearman `rho` are saved per residue,
  metric pair and **replica**, with paired-frame counts and time bounds.
- Missing series, insufficient paired frames and constant measurements are
  reported explicitly; neither frame interpolation nor zero imputation occurs.
- `pairs` can be a list of two-element supported metric names, or `[]` for
  an explicitly empty selection. `null` selects curated defaults.
- Correlations are descriptive and potentially autocorrelation-driven; they
  are **not** independent-replica tests, condition effects or causal evidence.
  Interpret individual replicate time traces before drawing conclusions.

## Existing results and readiness

- Studio and offline finalization share registered-product, time-grid, checkpoint identity/frame-count and prerequisite checks.
- Reuse requires compatible scientific settings and completed stages. Reporting labels/options or source revision alone do not invalidate data; changed sampling, source identity, selections or run settings do.
- Saved-only sources enter downstream catalogs only after **every configured replica** is verified. Reuse neither schedules trajectory work nor edits YAML; statistics records effective inputs.
- Existence cannot verify unregistered/unavailable products. Inspection can be expensive and does not certify chemistry, convergence or source content hashes.

The CLI rechecks inputs at execution, before replacing provenance or archived
results. Keep all intended source paths, condition/replica/run identities,
selections and downstream statistics/figures in the YAML:

```bash
mdscope statistics analysis.yaml --results-dir PATH
mdscope plotting analysis.yaml --results-dir PATH
```
