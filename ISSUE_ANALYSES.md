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

Each selected species must be present in the corresponding system; omit absent species or use system-specific YAML files. The lipid selection must contain a representative anchor for each physical molecule. If lipids in the topology consist of multiple bonded residues, designate only the anchoring residue/atom. Use your actual membrane normal axis.

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
    d_a_cutoff_A: 3.0
    angle_cutoff_deg: 150.0
    # hydrogens_sel: "..." # use explicit selection if auto-detection fails
    # acceptors_sel: "..."
```

This performs **one** MDAnalysis HydrogenBondAnalysis over all selected lipids, then distinguishes unordered lipid-type pairs (e.g. POPC–POPE and POPC–POPC), excluding within-the-same-residue bonds by default. This includes interspecies pairs as well as within-species bonds. If a glycolipid is represented by several covalently linked residues, exclusion by residue index is not equivalent to exclusion by whole lipid molecule. Such systems require molecule-specific selections/IDs for a truly intermolecular interpretation.

Events live in `interactions/lipid_hydrogen_bonds/events.parquet`; the common occupancy/lifetime statistics and interaction plot pipeline processes them identically to other H bonds.

## Testing and limitations

Tests include circular leaflet assignment across periodic z, an ideal tilted membrane plane, enrichment denominator correctness, carbon-profile selection, pair canonicalization, and progression fallback. Automated CI checks imports and the synthetic streaming trajectory too.

These tests **do not** establish physical validity on your large, heterogeneous neuronal bilayers. In particular: verify anchors, leaflet assignment, coordinate wrapping, headgroup chemistry and acceptor/donor guesses using TPR+XTC before treating output as manuscript-ready. With three independent replicates, interpret confidence intervals and between-replicate variability, not frame-pooled statistical significance.
