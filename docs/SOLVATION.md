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
    ion_selections:
      Na: "resname SOD"
      K: "resname POT"
      Cl: "resname CLA"
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
- Detailed: also saves solvent identity events and candidate bridges.
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

Current limitations include no solvent-excluded-volume normalization, no full
intermittent contact correlations, and no general water-acceptor bridge detector.

## Postprocessing and optional 3D water density

After the fused pass, `shell_candidates.parquet` records a candidate first
minimum of each bulk-normalized RDF with an explicit quality status. It does
not automatically change already computed coordination cutoffs: calibrate on
representative frames and fix the same chemical-class cutoff across compared
conditions before conducting kinetics. Ambiguous RDFs yield missing candidates.

`target_residues.parquet` maps target IDs to topology residue indices. When
a matching SASA run exists, `hydration_sasa.parquet` combines the **saved**
per-residue SASA with per-target solvent counts using exact frame time; missing
constituent residues produce NaN rather than partial SASA. Set `sasa_run` to
the selected run ID if needed. SASA exposure is a contextual covariate, not an
automatic correction to bulk-normalized RDFs.

When identity tracking is enabled, `survival.parquet` provides Kaplan–Meier
sampled-contact survival estimates per solvent/target. Left-censored episodes
are excluded; right-censored episodes remain at risk. Assess uncertainty across
simulation replicates, not across individual correlated contact episodes.

For optional locally aligned 3D hydration maps, specify exact target IDs and
**three ordered, noncollinear atoms** defining a moving coordinate frame:

```yaml
    grid_width_A: 12.0
    grid_spacing_A: 0.75
    grid_frames:
      protein/acidic_oxygen/42: "index 100 101 102"
```

The numeric ID and atom indices above are placeholders. Maps use the same
periodic candidate-water pairs as the radial calculations and accumulate
co-moving local density in `water_density_3d.parquet`. Only the inscribed
sphere (radius half the grid width) is sampled, and the frame origin must be
part of the target site selection. Missing voxels *inside* this sphere have
zero observed density; voxels outside are undefined, not zero. By design, this
is a sparse observed-density map: absent voxels correspond to zero observed
counts. Its `g_bulk` uses bulk density and does not adjust for inaccessible
local volume. A local map cannot be interpreted as an equilibrium fixed-frame
binding site if its neighboring molecular environment varies substantially.

## Reliability and performance validation

`tests/test_solvation.py` validates PBC, unique solvent counts and saved products.
`tests/test_solvation_postprocess.py` checks RDF minima, censored survival,
SASA joins, and chunk-level shell stability. `shell_stability.parquet` is a
convergence diagnostic computed entirely from already persisted histogram
chunks; the result is not an uncertainty interval from independent replicates.
A stable minimum needs four resolvable chunks with ≤0.5 Å cutoff spread.
Use a common validated cutoff for each chemical site class across systems.

For production workloads, first run standard mode with representative frames,
then measure throughput and peak memory using normal system profiling tools
on a single replicate. The optional event/bridge presets can produce large
outputs for all-water tracking; enable them only for selected chemical sites.
No production-scale performance benchmark has been measured in this PR.

## Optional geometry-qualified water bridges

Set `options.water_hbond_bridges: true` in the detailed or kinetics preset
(or any preset). The selected water oxygen must have exactly two covalently
bonded hydrogen atoms in the topology; the solvation target atoms must be
**hydrogen-bond acceptor atoms**. The implementation screens the existing
site–water neighbor pairs and checks each water O–H···acceptor angle (default
150°) and oxygen–acceptor separation (default 3.5 Å). It reports a
`water_hbond_bridges.parquet` event when **two different water hydrogens**
simultaneously donate to distinct targets. This intentionally does not cover
target-donor→water-acceptor hydrogen bonds or mixed-direction bridges.
Consequently, it is a *specific, conservative subset* of water hydrogen
bonding, not a comprehensive solvent-bridge detector. Incorrect acceptor
selections can produce false assignments; validate against topology chemistry.

```yaml
    options:
      water_hbond_bridges: true
      hbond_da_cutoff_A: 3.5
      hbond_angle_deg: 150
```


## Molecular aggregation and ion species

Set `ion_selections` to an explicit disjoint mapping such as
`{K: "resname POT", Cl: "resname CLA"}` to obtain separate `ions/K`
and `ions/Cl` outputs from separate small ion searches. The legacy
`ion_selection` aggregates all selected ions; use one or the other.
Explicit species selections permit distinct per-species RDF interpretation,
with the same search machinery and checkpointing.

`species_timeseries.parquet` reports per-frame mean coordination per
molecule/target in each configured class and the total of target counts;
`target_summary.parquet` reports each molecule/residue's mean, SD and
hydrated fraction. These are **not counts of distinct waters across an
entire chemical class**, because one water may coordinate two molecules.
Per-HS-chain grouping uses `group_by: hs_chain` and the already available
covalent HS-discovery implementation, not chainID or fragment grouping.
It requires full glycan bond topology and the `hs_auto` configuration;
MDscope reuses `sdc4_mechanisms.hs_auto` when provided.

See `docs/SOLVATION_SDC4.md` for an old-build SDC4 selection
template and `scripts/benchmark_solvation_real.py` for a **real observer**
preflight and benchmark against existing HPC XTC+TPR inputs. The GitHub
repository does not contain the full HPC trajectories, so publication-scale
timing and memory cannot be claimed from CI or synthetic benchmarks.
