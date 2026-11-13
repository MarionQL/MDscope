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

## Validate before long runs

```bash
mdscope validate analysis.yaml
mdscope run analysis.yaml
# Optional profiling:
pytest -q -s -m benchmark tests/test_performance_benchmarks.py
```

Preflight checks topology selections, group/alignment sizes, lipid-anchor ownership, SASA references, paths, run IDs, valid frame ranges and paired replica labels **before** replacing outputs. Overlapping protein groups are allowed; competing lipid species anchors are not.

`mdscope report STATISTICS_DIR` requires a **complete** manifest and contrast coverage table. Block settings live under `statistics.block.*`; retired aliases and `analyses.autocorrelation` / `analyses.distribution_statistics` must not be used.

[Scientific estimators](analyses.md) · [Workflow and restart](workflows.md) · [Scientific contracts](developer/scientific-contracts.md)
