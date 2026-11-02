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

The plotting finalization phase stamps **only newly written PNGs** in
`project_summary/report/figure_provenance.json`. Offline PDF and ZIP generation
reject unverified legacy PNGs, disabled families/kinds, mismatched render
configurations, changed source tables, missing sources and modified images.
The complete figure manifest retains rejection reasons, source metric/units,
generation configuration hash, timestamp, family-level source paths and actual
replica frame/time coverage. PNGs have SHA256 content hashes; source freshness
uses path/size/mtime fingerprints, explicitly not content hashes. This keeps
validation cheap without scanning multi-million-row interaction tables.

Migration: regenerate the plotting phase to obtain verified figures before
exporting an old result directory. Files are not deleted. Changing PDF section
selection alone does not invalidate figures; changing source code or plot/analysis/statistical
settings requires rerendering. Do not stamp an old directory by hand to bypass
source verification. Original trajectory-parameter provenance remains in the
original checkpoints/run records; the current rendering configuration does not
prove the original analysis configuration. This limitation is explicit in
`report_provenance.json`, rather than inventing provenance for legacy tables.

PDF tables use readable SDC4 labels and units while the bundle preserves raw
identifiers. `Unavailable` denotes missing/nonfinite data; physical zeros stay
zero. Scalar comparison labels use adjusted q-values for support, never the
sign of a bootstrap interval. Disabled PDF sections do not leave separator-only
pages. `report_provenance.json` and `figure_manifest.csv` support offline review;
no trajectory is reopened to regenerate the report.
