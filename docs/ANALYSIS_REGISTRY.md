# Analysis instances and grouping

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

## SASA: named targets and environments

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
