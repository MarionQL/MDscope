# Configuration and analysis registration

This guide consolidates the authoritative technical details formerly scattered across topic and implementation-history files. Sections retain their scientific caveats and validation limitations; consult the source code for the executable contract.


---

### Analysis instances and grouping

MDscope uses **one shared group definition**: `selections.groups`. A group specifies
its atom selection and optional independent alignment selection. The streaming
reader computes group RMSD, RMSF and radius of gyration on the same trajectory
frames and writes rows distinguished by `group`; it does not create separate readers.

```yaml
selections:
  groups:
    chain_a:
      selection: "segid PROA"
      alignment_selection: "segid PROA and name CA"
    chain_b:
      selection: "segid PROB"
      alignment_selection: "segid PROB and name CA"
analyses:
  rmsd:
    enabled: true
    groups: [chain_a, chain_b]
  rmsf:
    enabled: true
  radius_of_gyration:
    enabled: true
    groups: [chain_a]
```

For RMSD and radius of gyration, optional `groups` filters which configured
groups receive that metric; omitted means all groups. The global default
metric remains in `summary/timeseries.parquet`; group measurements are in
`summary/group_timeseries.parquet`. RMSF accumulators and group outputs
continue to use the shared group definitions.

Other multi-instance analyses (currently density and residue–membrane distance)
use `analysis_runs(cfg, name)` with named runs, because their instances have
distinct configuration and output structures. Do not add a second RMSD/Rg
`runs` configuration on top of `selections.groups`.

For fused observers, `registry.py` describes expected persisted artifacts.
The shared trajectory pass remains the execution engine. Every independently
configured instance needs its own output identity and verifiable completion
state. Existing algorithms, PBC rules and coordinate spaces must be preserved.

### SASA: named targets and environments

SASA supports named runs that all share the existing trajectory pass.
Each run computes **per-frame, per-residue target SASA**, with all selected
occluders included in the FreeSASA geometric surface calculation.
The target atoms **must** also be included in the occluder selection.

```yaml
analyses:
  sasa:
    enabled: true
    runs:
      isolated:
        target_selection: "protein"
        occluder_selection: "protein"
        n_points: 200
      membrane:
        target_selection: "protein"
        occluder_selection: "protein or resname POPC POPE POPI DOPC DPPC PSM CHL1 BGAL GM1a CER160"
        periodic: true
        reference_run: isolated
        n_points: 200
```

Output for each named run: `profiles/sasa/<run_id>/residue_sasa.parquet`,
with columns `time_ps`, `topology_residue_index`, `residue`, `sasa_nm2`.
A legacy single unnamed run retains `profiles/sasa/residue_sasa.parquet`.
With `reference_run: isolated`, MDscope also writes `residue_burial.parquet`
(per-frame isolated/environment SASA, buried nm² and buried fraction) plus
`residue_burial_summary.parquet` (means and standard deviations by residue).
The target atom identities, radii, probe and algorithm settings must match.
Buried fractions are NaN when isolated SASA is zero.
FreeSASA parameters and probe radius may be configured on each run.

**Scientific conventions:** Nonstandard atoms are included using a documented
element-wise geometric van der Waals radius map (Angstroms): H 1.20,
C 1.70, N 1.55, O 1.52, P/S 1.80, F 1.47, Cl 1.75, Br 1.85, I 1.98.
Override with `element_radii_A` for your force-field/radius conventions.
When topology elements are absent (common in PSF/GROMACS files), basic
H/C/N/O/P/S/F/I and Cl/Br identities are inferred from atom names. Ambiguous
names and unrecognized elements must be resolved explicitly in the topology;
check atom typing before trusting heterogeneous-membrane SASA. These are geometric SASA radii,
not CHARMM force-field Lennard-Jones radii.

For `periodic: true`, occluder coordinates are minimum-imaged around the
target centroid using the trajectory cell. The target must be a compact
whole molecule and fit within half the box along each Cartesian dimension;
otherwise the calculation raises an error. This does not reconstruct a
protein split across PBC. Prepare a whole-molecule trajectory first.
Environment-occluded SASA is **not equivalent** to direct water hydration.

For large occluder selections, each frame uses a radius/probe-bounded neighbor
search to omit atoms too far from any target SASA sphere to contribute to the
result. The single trajectory pass is preserved; this is not a second frame
reader. Always benchmark and compare with a full-occluder reference frame
before using a new topology/radius convention in production.


### Independent statistics families

`statistics.families` accepts `time_series`, `profiles`, `interactions`, and
`solvation_rdfs`. Values may be booleans or mappings containing `enabled`,
`bootstrap_draws` (>=100), `ci`, and `hypothesis_tests` (boolean). Omitted
families retain the historical behavior. Overrides apply to one family and
never mutate the trajectory configuration. Example:

```yaml
statistics:
  families:
    time_series: {enabled: true, bootstrap_draws: 1000}
    profiles: false
    interactions: {enabled: true, hypothesis_tests: false}
    solvation_rdfs: true
  interaction_include_kinds: [salt_bridges, hydrogen_bonds]
  interaction_exclude_kinds: [lipid_hydrogen_bonds]
  interaction_progress_every: 1000
  hypothesis_families:
    primary: [sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction]
```

Inclusion/exclusion occurs before interaction Parquet reads; exclusions take
precedence. Scalar hypothesis-family membership partitions BH adjustment by
condition contrast and family. All unlisted endpoints remain exploratory.
Define families before follow-up inference, never after inspecting p-values.
Without this configuration, the original global scalar BH family is retained.
A two-sided independent exact test with 3+3 replicas has minimum p=2/20=0.1;
a narrow bootstrap interval is not FDR-adjusted significance.

Gapped/nonfinite or irregular scalar traces retain finite-observation point
estimates but have unavailable ACF/effective-sample estimates and temporal CIs.
System temporal intervals are also unavailable when any constituent replica
has an invalid time axis; an invalid replica is not silently dropped.
`statistics/manifest.json` records methods and enabled/disabled families.
This does not repair missing-not-at-random sampling or establish convergence.

Run optional timing benchmarks with
`pytest -q -s -m benchmark tests/test_performance_benchmarks.py`.
Normal pytest excludes timings; direct COM/MSD and planar-normal references
remain in the correctness suite. Packaging smoke checks are consolidated in
`test_package_contract.py`; wheel/sdist and outside-checkout checks remain CI
installation gates. Existing config persistence/provenance exports remain
supported and their atomic-write/rollback tests are retained.


Invalid forward frame windows are rejected before opening trajectories;
negative Python frame indexing remains supported. With invalid/gapped scalar
series, omitted temporal draws intentionally change subsequent scalar RNG
progression compared with the old compressed-time behavior. Complete regular
data retains its draw ordering; the interaction grouping optimization and its
RNG sequence are unaffected by this scientific correction.
