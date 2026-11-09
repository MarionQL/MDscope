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


### Plot families, PDF sections and figure freshness

Plot generation and PDF inclusion are independent from raw collection and
statistics. Omitted controls retain enabled behavior:

```yaml
plots:
  enabled: true
  families:
    time_series: true
    profiles: false
    interactions: true
    sdc4: true
    shared_states: true
    lipid_enrichment: true
    diffusion: true
  interaction_exclude_kinds: [lipid_hydrogen_bonds]
report:
  enabled: true
  sections:
    overview: true
    time_series: true
    shared_states: true
    solvation: true
    sdc4: true
    lipid_enrichment: true
    methods: true
    atlas: true
  interaction_exclude_kinds: [lipid_hydrogen_bonds]
```

A false plot family prevents its plot entry point from reading inputs;
interaction kinds are excluded before raw identity tables are read. Existing
broad `plots.time_series`, `profiles`, `interactions`, and `enabled` switches
remain honored. Diffusion summaries still compute persisted coefficients when
only their *plot* is disabled. Shared-state numerical modeling remains governed
by its analysis setting, independently of SDC4/shared-state plot switches.

Plotting records source paths, sampling, units and generation time for newly
written PNGs in `project_summary/report/figure_provenance.json`. Offline PDF
and ZIP export can use readable older PNGs and marks unknown provenance as
`legacy_unverified`. It excludes disabled families, corrupt PNGs and figures
whose recorded sources are missing. Changing code or render settings does
not force regeneration. A readable figure and available sources do not prove
that the figure was made from the current versions of those sources.

PDF tables use readable SDC4 labels and units while the bundle preserves raw
identifiers. `Unavailable` denotes missing/nonfinite data; physical zeros stay
zero. Scalar comparison labels use adjusted q-values for support, never the
sign of a bootstrap interval. Disabled PDF sections do not leave separator-only
pages. `report_provenance.json` and `figure_manifest.csv` support offline review;
no trajectory is reopened to regenerate the report.
### Independent statistics families

`statistics.families` accepts `time_series`, `profiles`, `interactions`, and
`solvation_rdfs`. Values may be booleans or mappings containing `enabled`,
`bootstrap` (boolean), `bootstrap_draws` (>=100), `ci`, `hierarchy`, and
`hypothesis_tests` (boolean). Omitted
families retain the historical behavior. Overrides apply to one family and
never mutate the trajectory configuration. Example:

```yaml
statistics:
  families:
    time_series: {enabled: true, bootstrap_draws: 1000}
    profiles: false
    interactions: {enabled: true, bootstrap: false, hypothesis_tests: false}
    solvation_rdfs: true
  interaction_include_kinds: [salt_bridges, hydrogen_bonds]
  interaction_exclude_kinds: [lipid_hydrogen_bonds]
  interaction_progress_every: 1000
  hypothesis_families:
    primary: [sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction]
```

Optional adaptive Monte Carlo precision is available for scalar
`time_series` system and contrast intervals. Set
`statistics.families.time_series.adaptive_hierarchy` to a mapping such as:

```yaml
min_draws: 1000
max_draws: 8000
batch_size: 500
stable_checks: 2
relative_tolerance: 0.05
```

Each row checks both endpoints of its mean and median intervals. The tolerance
is the optional `absolute_tolerance` in the metric's units plus
`relative_tolerance` times the corresponding interval width at the first look.
At least one tolerance must be positive. Simultaneous binomial rank bands
account for the scheduled looks and endpoints within each row. Rows report
actual draws, achieved endpoint bounds, settings, seed and status. A capped
row reports `not_converged` and retains its estimated interval.

The bands measure second-stage Monte Carlo precision conditional on the
fixed-size replicate bootstrap pools (`bootstrap_draws`). Replicate intervals,
other families, and the scientific interpretation of uncertainty retain their
existing contracts. Adaptive mode uses a separate seeded RNG stream and is
reproducible for a fixed input order and configuration.

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

`statistics.bootstrap: false` retains point estimates and disables temporal
bootstrap draws. A family may override it independently. `hypothesis_tests`
remains independent: point estimates and permutation p/q values can be
requested without computing bootstrap intervals. Uncomputed intervals and
directional bootstrap probabilities are NaN, with `bootstrap_draws: 0` in
their output rows; no point estimate is relabeled as a confidence interval.

Each family's `hierarchy` mapping may override `enabled`,
`block_within_replicate`, and `resample_replicates`. Setting `enabled: false`
disables system/contrast bootstrap intervals while retaining replica-level
intervals and descriptive equal-replica system estimates. Setting
`block_within_replicate: false` selects block length 1, an explicit temporal
independence assumption; replicas remain the independent units for tests.
The RDF family is descriptive mean/SD only, so it does not acquire bootstrap
intervals or per-bin hypothesis tests from these settings.

Profiles use the same regular-time requirement as scalar series. Missing or
nonfinite columns keep finite-observation descriptive points but cannot produce
temporal intervals; a system/contrast interval is unavailable at a position if
any requested replica is missing or invalid there. A paired contrast with no
complete pairs is unavailable, rather than an empty-array failure. RMSF coordinates with missing
samples likewise cannot provide a temporal interval. Labels align columns
across replicas before aggregation; they are not silently truncated by width.
Complete regular-data defaults retain their numerical results and RNG order.
Fixed configurations are reproducible. Opting out of a calculation omits its
random draws, so later Monte Carlo work within that family may use a different
subsequence; other families retain their separately seeded streams.

Statistical reruns archive older family tables and method manifests under
`.archives/statistics`, including when a family is disabled or now has no data.
The active manifest is marked `running` until every requested family finishes;
offline reports reject incomplete statistical products. Raw replica outputs
are retained. `report_provenance.json` embeds the active statistical methods.

### Failure checks before trajectory analysis

`validate`, `run`, and `replicate` share topology-only selection validation.
It checks protein/features, alignment size, groups, missing/overlapping lipid
anchors, density and distance selections, and SASA target/occluder/reference
identity using the same SASA setup code without requiring FreeSASA for
validation. `replicate` preflights only its selected system. Failed checks
precede project initialization, output replacement and trajectory processing.
Offline statistics/report regeneration does not invoke this topology preflight.

Frame indices, stride and chunk size must be integers. YAML sections must be
mappings; input paths must be files, named analysis IDs must be filesystem-safe,
and paired condition designs require matching replica labels. Overlapping
analysis groups remain valid distinct outputs (e.g. a protein and a domain);
overlapping lipid species anchors are invalid molecule ownership.

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

### Offline evidence export integrity and layout

`mdscope report` retains the existing PDF/ZIP filenames and does not rerun
trajectory analyses or statistics. Report orchestration, persisted-table
inventory, PDF presentation and ZIP export now have separate modules. ReportLab
is imported only when rendering a PDF; label formatting and bundle code do not
import plotting or trajectory readers merely to load those modules.

Statistics manifests record completion, enabled families and methods. Offline
reports use readable persisted tables and the methods recorded with them.
Missing method metadata is labeled unverified, and an unrecorded confidence
level is left unspecified. Incomplete or malformed manifests are excluded.
To apply new statistical settings, explicitly rerun statistics. A software
upgrade or changed report configuration does not require that rerun.

`table_manifest.csv` records numerical table eligibility;
`report_content.json` records displayed file-backed figures and project tables.
The figure manifest also identifies byte-identical duplicate PNGs. The PDF
and ZIP use the same eligibility policy, and the ZIP contains this export's
manifests, original eligible project Parquet tables and CSV exports (CSV
materialization is skipped above 50 MB; the Parquet is still preserved).
A PDF or bundle construction failure preserves the previous successful
PDF/ZIP and report sidecars. Publication is staged locally; interruption during
the final individual file replacements is not a multi-file transaction.

Scalar tables separate scientific labels from units, keep raw identifiers in
companion data, show the configured CI level, and place estimate, interval,
p/q and inference status on separate lines. The font is 8 pt with controlled
columns and no midword splitting. The contrast names the actual ordered
systems; support still depends on adjusted q, never interval sign. Unequal
replica windows/frame counts are flagged in the sampling section. Original
six-replica production pages and all plot families still require human review;
a synthetic layout check does not close that scientific QA requirement.
