# Open-issue analyses (#1–#5)

All additions are opt-in and share MDscope's existing replicate-aware statistics, profile bootstrap and visualization. You can run all selected passes with `python -m mdscope.analysis.cli run config.yaml`.

## Issue #1: all-atom C–H order

Existing `membrane_order.py` computes conventional **-S_CH**, sn-1/sn-2 carbon-position profiles, averages and per-species metrics. This PR extends it to upper/lower leaflets (assigned once by the first analyzed configuration), with position-3 series and true bonded C–H connectivity where bonds are present. The original `membrane_order` configuration remains valid, plus:

```yaml
analyses:
  membrane_order:
    enabled: true
    axis: z
    per_leaflet: true
    per_lipid_type: true
    leaflet_anchor_selection: "(resname POPC POPE POPI DOPC DPPC PSM and name P)"
```

When a lipid residue is absent from the optional `leaflet_anchor_selection`, a proximal tail carbon is used. In heterogeneous glycolipid systems, specify a meaningful anchor for every analyzed lipid where possible. Profile plots and block-bootstrap statistics use the same existing pipeline.

## Issue #2: progress

All shared streaming, lipid-order, leaflet-normal, and enrichment loops report progress. An interactive terminal uses `tqdm` when installed; noninteractive Slurm output prints periodic percentage milestones without requiring an extra dependency. Disable using:

```yaml
progress:
  enabled: false
```

## Issues #3 and #4: normals and local enrichment

These analyses share one optimized atom-to-lipid anchor reduction. Define **one anchor per physical lipid**; for ceramide-based glycosphingolipids use an identifiable ceramide anchor rather than selecting every glycan residue. The example below is illustrative: adapt residue and chain names to the actual CHARMM-GUI topology and validate with a short trajectory.

```yaml
analyses:
  bilayer_normal:
    enabled: true
    axis: z
    neighbor_radius_A: 15.0
    min_neighbors: 5
    species_selections:
      POPC: "resname POPC and name P"
      POPE: "resname POPE and name P"
      CHL1: "resname CHL1 and name O3"
  lipid_enrichment:
    enabled: true
    axis: z
    protein_selection: "protein and name CA and resid 146-170 344-368"
    inner_radius_A: 0.0
    outer_radius_A: 15.0
    species_selections:
      POPC: "resname POPC and name P"
      POPE: "resname POPE and name P"
      CHL1: "resname CHL1 and name O3"
```

By default every selected species must be present. For paired system configurations where GM1a/BGAL may be intentionally absent in the control, use `allow_absent_species: true` to skip missing species while still requiring at least one lipid species; missing categories are not silently interpreted as zero enrichment. The lipid selection must contain a representative anchor for each physical molecule. If lipids in the topology consist of multiple bonded residues, designate only the anchoring residue/atom. Use your actual membrane normal axis.

**Normals**: fit a local plane to neighboring anchor atoms within the same leaflet using periodic coordinates and report its angular deviation from the configured global axis in degrees. This measures **local interfacial tilt/undulation**, not a lipid acyl-chain order parameter or a mechanical modulus. Lipids with insufficient neighbors are recorded as missing.

**Enrichment**: use a cylindrical (projected lateral) annulus about the protein selection. In each leaflet and frame, compute

```
fold enrichment = (N_species,local / N_all,local) /
                  (N_species,bulk  / N_all,bulk)
```

where bulk means all selected lipids in the same leaflet, not an external experimental bulk sample. Values above 1 indicate enrichment; below 1 indicate depletion. A zero observed count yields fold 0, while an empty local annulus yields NaN. Outputs preserve numerator and denominator for audit. Local and bulk populations overlap by design, so this is a **compositional** enrichment, not a binding free energy.

- `extensions/bilayer_normal/timeseries.parquet` and optional chunked `per_lipid/angles_part_*.parquet` (enable `save_lipid_angles`)
- `extensions/lipid_enrichment/timeseries.parquet`, `counts.parquet`, `selection_audit.parquet`

Scalar metrics automatically enter the existing `summary/timeseries.parquet`, project statistics and standard time-series figures.

## Issue #5: lipid–lipid H bonds

```yaml
analyses:
  lipid_hydrogen_bonds:
    enabled: true
    lipid_selection: "resname POPC POPE POPI DOPC DPPC PSM CHL1"
    exclude_same_residue: true
    exclude_same_molecule: true
    d_a_cutoff_A: 3.0
    angle_cutoff_deg: 150.0
    # hydrogens_sel: "..." # use explicit selection if auto-detection fails
    # acceptors_sel: "..."
```

This performs **one** MDAnalysis HydrogenBondAnalysis over all selected lipids, then distinguishes unordered lipid-type pairs (e.g. POPC–POPE and POPC–POPC). It now excludes **all bonds within the same covalently connected lipid molecule** by default, even when that molecule spans multiple glycan residues. The molecule IDs are derived from topology bonds; chain IDs identify glycolipid *species*, not individual molecules. Select complete lipid molecules, not only their carbohydrate headgroups, so the connected-component graph is complete. The output includes `molecule_audit.parquet` with counts of excluded same-residue and cross-residue intramolecular events.

Events live in `interactions/lipid_hydrogen_bonds/events.parquet`; the common occupancy/lifetime statistics and interaction plot pipeline processes them identically to other H bonds.

## Testing and limitations

Tests include circular leaflet assignment across periodic z, an ideal tilted membrane plane, enrichment denominator correctness, carbon-profile selection, pair canonicalization, and progression fallback. Automated CI checks imports and the synthetic streaming trajectory too.

These tests **do not** establish physical validity on your large, heterogeneous neuronal bilayers. In particular: verify anchors, leaflet assignment, coordinate wrapping, headgroup chemistry and acceptor/donor guesses using TPR+XTC before treating output as manuscript-ready. With three independent replicates, interpret confidence intervals and between-replicate variability, not frame-pooled statistical significance.

## Validated Syndecan-4 TPR selections (test-data-v1, October 2026)

Real GROMACS TPR/XTC smoke tests on both published reference systems confirmed that **the glycolipid molecule names are chain IDs, not residue names**. Do not use `resname GM1a` to select ganglioside molecules: it matches zero atoms in the glycolipid-containing topology. Likewise, `resname BGAL` includes carbohydrate residues outside the BGAL glycolipid chains (e.g., HS glycans), and is not a molecule-specific lipid selection.

For the published test systems, the following ceramide anchors each identify one physical glycosphingolipid/ceramide molecule:

```yaml
analyses:
  bilayer_normal:
    enabled: true
    allow_absent_species: true
    species_selections:
      POPC: "resname POPC and name P"
      POPE: "resname POPE and name P"
      POPS: "resname POPS and name P"
      POPI: "resname POPI and name P"
      DOPC: "resname DOPC and name P"
      DPPC: "resname DPPC and name P"
      PSM: "resname PSM and name P"
      CHL1: "resname CHL1 and name O3"
      GM1a: "chainID GM1a and resname CER160 and name NF"
      BGAL: "chainID BGAL and resname CER160 and name NF"
      CER160: "chainID CER160 and resname CER160 and name NF"
```

Reuse the same `species_selections` for `lipid_enrichment`. In the **+headgroup** test topology, there are 24 GM1a and 58 BGAL ceramide anchors and no free CER160 chain anchors. In the **no-headgroup** control, there are 82 CER160 chain anchors and no GM1a/BGAL chain anchors. `allow_absent_species: true` is required when using one shared selection map across both systems.

This mapping is specific to the uploaded TPRs and should be revalidated for other CHARMM-GUI builds. The lipid–lipid H-bond routine now uses bonded connected components to exclude intramolecular events. For full glycosphingolipids, use a complete selection such as `(resname POPC POPE POPI POPS DOPC DPPC PSM CHL1) or (chainID GM1a BGAL CER160)`; **do not** use `resname BGAL` alone, because it also selects HS-associated carbohydrate residues. A topology with explicit covalent bonds (e.g. TPR) is required.

GitHub Actions runs the optional `tests/real_trajectory_smoke.py` against the two reference systems. The initial eight-frame smoke tests confirmed runtime compatibility and finite numeric outputs; they do **not** validate physical convergence, long-lag kinetics, or statistical significance.

## Full-length physical audit and performance profiling

`tests/real_trajectory_physics.py` reads **all** frames of both published TPR/XTC test systems for leaflet continuity, local normal fitting, PBC invariance, and lipid enrichment. It explicitly shifts subsets of lipid and protein atoms by entire box vectors and checks that the physical observables do not change. It compares the naive protein center to the circular, periodic-aware center and records how often their enrichment counts differ. A separate 20-frame hydrogen-bond test uses complete GM1a/BGAL/CER160 chains and checks bonded molecule counts and removal of within-molecule H-bonds.

Per-analysis timing is recorded in `<system>/<replicate>/stage_timings.json`. Full-frame audit timings and peak resident memory are in `physical_audit.json` and `<system>/physics/geometry_per_frame.parquet`. These are single-run wall times on GitHub-hosted CPU hardware, not portable throughput benchmarks. They are useful for prioritizing improvements but should be remeasured on the intended HPC system.

**Interpretation:** a 200-frame audit tests input geometry, selection stability and PBC invariance, not independent-replicate statistical significance or convergence. For protein-centered enrichment, prefer a transmembrane-domain selection rather than the whole extracellular protein when the biological question is membrane recruitment. The lipid-order module still assigns leaflet profiles from the first analyzed frame; check the leaflet-change audit before interpreting those profiles on trajectories with actual flip-flop.
