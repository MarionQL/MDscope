# Configuration

Use [`synd4_full_mechanisms.yaml`](../synd4_full_mechanisms.yaml) for the complete, maintained SDC4 configuration. The examples below show **patterns**, not validated selections for another topology.

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
- **Named static figures:** `solvation:shell_rdf`, `solvation:water_and_ions`, `solvation:contact_survival`, `diffusion:diffusion_coefficients`, `lipid_enrichment:lipid_enrichment_heatmap`, `lipid_enrichment:lipid_area_density_heatmap`.
- **Shared-state figures:** `shared_states:state_populations`, `shared_states:state_sampling_diagnostics`.
- **SDC4 condition figures:** `sdc4:<system_id>:<figure_stem>`; supported stems include `nxip_height_occlusion_joint`, `vvhplvpl_proline_ramachandran`, `lipid_density_maps`, `lipid_enrichment_maps`, and the five `lipid_region_*` charts. A single PNG is saved under `project_summary/figures/sdc4/<system_id>/`; each independent replica appears as an explicitly labelled panel or trace. The *numerical* Parquet/NPZ measurements remain per replica in `systems/<replicate>/extensions/sdc4_mechanisms/`. Previous per-rep diagnostic PNGs from older runs are not automatically deleted; archive them when regenerating plots. No trajectories or statistical CIs are pooled to draw these condition diagnostics.
- **Lipid-contact figures:** `lipid_contacts:<figure_stem>` for each saved `lipid_contact_matrix_*` or `lipid_contact_difference_*` figure; valid endings are `_neighbors_per_molecule` or `_enrichment`.

Per-figure `plots.options.<plot ID>.dpi` is supported for resolution overrides (integer 50–600); it does not change scientific results, axis definitions or statistical estimands. Unsupported per-plot keys fail configuration validation.

Plot IDs are matched exactly, not by wildcard or substring. Unrecognized figure IDs fail YAML validation. An enabled ID does **not** generate missing statistical inputs; source data and the corresponding family prerequisites remain necessary. The plot ID format is checked during YAML validation; data-dependent metrics, run IDs, species and systems are verified only when persisted inputs exist. An unknown data-derived identifier cannot be rejected merely from its syntax. When selecting solvation by an explicit ID, no separate `plots.families.solvation` opt-in is needed. `plots.include: []` disables all figures.

### Compact statistics-table PDFs

Enable `plots.table_reports_pdf: true` to collect existing bootstrap and
between-condition comparison tables into one multipage
`statistics_tables.pdf` **per figure directory** (time_series, profiles,
and each interaction type). When `plots.bootstrap_tables: true`, bootstrap
pages are included; comparison pages still depend on available statistics.

The report replaces individual `*__bootstrap_table.png` and
`*__comparison_table.png` outputs. Previously generated PNG tables in
the same directories are moved to
`<project_output>/.archives/figure_tables/` only after the new PDF
has been written successfully. Other PNG figures and the full-precision
Parquet statistics are retained. `plots.include` and disabled families
continue to control which tables appear. The PDF is for reading and
sharing; it is not a substitute for numerical data exports.
### Time-series rolling averages

The standard per-replica time-series traces can show a trailing, time-based
rolling mean. The window is measured using the saved `time_ps` timestamps,
not a fixed frame count, so irregular sampling does not change its duration.

```yaml
plots:
  enabled: true
  bootstrap_tables: true
  table_reports_pdf: true
```

`table_reports_pdf: false` (default) preserves individual table PNGs.
  time_series: true
  time_series_rolling_window_ns: 100
```

`null` (default) plots the unsmoothed values. A positive finite number
enables the trailing average over `(t - window, t]`; the initial part of
each trace uses the available samples. Missing/nonfinite observations do
not contribute to the mean; an entirely missing window is not plotted.
Each replicate is smoothed separately. This is **display-only**:
distributions, autocorrelation, bootstrap statistics, CIs and saved
trajectory-derived quantities remain unchanged.

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

## Preview the work without analyzing a trajectory

```bash
mdscope plan analysis.yaml                    # schema-only; no TPR/XTC required
mdscope plan analysis.yaml --format json      # machine-readable preview
mdscope plan analysis.yaml --check-topology   # real topology selections/counts
mdscope plan analysis.yaml --check-inputs     # require TPR and trajectory paths
```

The plan reads the existing observer output registry and identifies required per-replica artifacts, named analysis runs, optional extra trajectory passes, enabled statistics families and figure families. Paths are **relative to each replica**; nothing is created or executed. Schema-only mode still validates YAML options and analysis constraints but does not prove that selections, chemistry, files, or PBC assumptions are valid. Topology mode inspects the topology only, not XTC frames.

This is a **partial** catalog: many raw/cache outputs, exact Parquet/NPZ schemas, and individual plot IDs are not yet registered. A listed statistics or plot family is conditional on available saved inputs, not a promise of a nonempty result or a confidence interval. Use `mdscope validate` before production runs.

## Validate before long runs

Unknown YAML option names, including nested analysis and statistics options, now fail with a dotted path and a suggested correction when available. This check also applies to disabled analyses. System, group, species, motif, solvation-site and named-run identifiers remain user-defined; their option fields are still validated.

```bash
mdscope validate analysis.yaml
mdscope run analysis.yaml
# Optional profiling:
pytest -q -s -m benchmark tests/test_performance_benchmarks.py
```

Preflight checks topology selections, group/alignment sizes, lipid-anchor ownership, SASA references, paths, run IDs, valid frame ranges and paired replica labels **before** replacing outputs. Overlapping protein groups are allowed; competing lipid species anchors are not.

`mdscope report STATISTICS_DIR` requires a **complete** manifest and contrast coverage table. Block settings live under `statistics.block.*`; retired aliases and `analyses.autocorrelation` / `analyses.distribution_statistics` must not be used.

[Scientific estimators](analyses.md) · [Workflow and restart](workflows.md) · [Scientific contracts](developer/scientific-contracts.md)
