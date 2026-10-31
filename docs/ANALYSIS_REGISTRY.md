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
