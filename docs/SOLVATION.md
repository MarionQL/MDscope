# Fused solvation analysis (optional)

The solvation observer shares MDscope's raw-frame fused processing. It runs one
periodic capped-distance query per solvent class and frame, writes atomic chunks,
and reuses them for offline counts, bulk-normalized radial distributions and
sampled-frame contact episodes.

```yaml
analyses:
  solvation:
    enabled: true
    preset: standard  # standard | detailed | kinetics
    water_selection: "resname TIP3 and name OH2"
    ion_selection: "resname SOD and name SOD"
    targets:
      protein:
        acidic_oxygen:
          selection: "protein and (resname ASP GLU) and name OD1 OD2 OE1 OE2"
          group_by: residue
      hs:
        sulfate:
          selection: "resname AGLCNA AIDOA BGLCA and name S1O S2O S3O"
          group_by: residue
      lipid:
        headgroup:
          selection: "resname POPC POPE and name O11 O12 O13 O14"
          group_by: residue
    search_radius_A: 8.0
    rdf_bin_width_A: 0.2
    shell_cutoffs_A: {water: 3.5, ions: 4.0}
```

**Selections are examples, not validated for your topology.** Check CHARMM atom
names for water, ions, sulfate groups, lipids and glycolipids. For lipids
spanning several topology residues choose `group_by: fragment` when bonds
identify intact molecules; multiple molecules in a disconnected fragment or
missing bonds must be corrected in the topology, not guessed spatially.
Target selections must be disjoint.

- Standard: compact per-target hydration/ion counts and site-class RDFs.
- Detailed: also saves selected solvent identity events and candidate bridges.
- Kinetics: same identity tracking for sampled-frame episode analysis.
- `options: {events: false, bridges: false}` can override preset behavior.

Outputs under `extensions/solvation/`: `counts.parquet`,
`rdf.parquet`, optional `events.parquet`, `bridges.parquet`,
`episodes.parquet`, and `metadata.json`. A bridge is the same solvent
entity near two distinct targets in one sampled frame; it is **not** a
hydrogen-bond-qualified bridge. Episodes store sampled-frame contact runs;
duration zero for a single observed frame does not mean zero physical lifetime.
Censor flags show whether contacts touched analyzed boundaries. Genuine
high-frequency exchange cannot be inferred from coarsely saved trajectories.

The RDF is normalized to the **bulk** number density computed from total
selected waters/ions and periodic cell volume. This is **not corrected**
for the solvent-excluded volume of protein, glycans and bilayer. Buried sites
can exhibit lower radial density solely because water is inaccessible.
Interpret RDF alongside SASA and occlusion, and calibrate first-shell
boundaries using stable first minima and cutoff sensitivity analysis.
Shell cutoffs are explicit and are *not automatically fitted* in this version.
The pooled site-class RDF must not be mistaken for per-residue density.

This version does not yet calculate 3D local water density, hydrogen-bond
angles, solvent-excluded-volume-normalized RDF, intermittent contact
correlations, or Kaplan–Meier survival estimates. Those require additional
validated definitions; episodes are the sufficient sampled-contact input
for later continuous survival estimation.
