# Scientific analyses and methods

This guide consolidates the authoritative technical details formerly scattered across topic and implementation-history files. Sections retain their scientific caveats and validation limitations; consult the source code for the executable contract.


---

### Neuronal interface analyses

The canonical architecture and package locations are documented in
[docs/ARCHITECTURE.md](developer/architecture.md). The rules for restructuring
are in [docs/REPOSITORY_RULES.md](developer/contributing.md).

This document focuses on the neuronal-interface configuration and scientific
interpretation of optional analyses, not an additional module map.

### Adding new analyses

1. Decide the output shape **before** writing code:
   - scalar time series: `time_ps` plus one column per quantity, with optional `group`;
   - frame-by-residue/profile: positions, labels, time-indexed matrix and reducer;
   - interaction: event key, event frame/time, identity and occupancy;
   - independent per-replicate parameter (e.g., diffusion slope): one row per *replicate*, not per frame.
2. Add a module in its appropriate `analysis/` subpackage with a top-level `run(rep,cfg,repout)` function for trajectory-specialized analyses. Use `MDAnalysis` and shared geometry tools; vectorize expensive distance calculations.
3. Register the module in `analysis/registry.py` within `specialized_tasks`. Persist scalar observables as `extensions/<analysis_name>/timeseries.parquet`, and register their file in `additional_time_series`.
4. This output feeds the *existing* `summary/timeseries.parquet`, `analysis/statistics.py` and `analysis/plotting.py`. **Do not implement custom hypothesis testing or frame-pooled error bars** in individual modules.
5. Profile analyses should expose arrays to `analysis/statistics.py:collect_temporal_profiles`; interaction events should use `analysis/interaction_tables.py`. Add a small dataset test and import test before production use.
6. Expensive specialized analyses should operate as separate restartable passes. Cheap per-frame observables should share `analysis/pipeline/streaming.py` so large XTC files are not opened once per scalar observable.

### Diffusion coefficients from existing MSD

Keep the original YAML section:

```yaml
analyses:
  membrane:
    enabled: true
    axis: z
    headgroup_selection: "name P"
    msd:
      enabled: true
      selection: null
      include_3d: false
      slope_start_fraction: 0.2
      slope_end_fraction: 0.8
```

No `lipid_diffusion` analysis is needed or supported by the modular extension. For every replicate, `observables/msd.py` reads the already consolidated `membrane_msd_lateral_A2` series and fits the linear regime. The factor is **2.5** to convert a lateral slope measured in Å²/ps into nm²/ns, using `MSD = 4Dt`.

Outputs:
- `<system>/<replicate>/summary/diffusion_coefficients.parquet`: D, fit interval, slope and R².
- `project_summary/diffusion/replicate_coefficients.parquet`: all replicate-level coefficients.
- `project_summary/diffusion/system_coefficients.parquet`: descriptive summaries, *no significance claims*.
- `project_summary/diffusion/diffusion_coefficients.png`: independent replicate points and means.

This remains a **single-origin** MSD estimator. It is not species-resolved or a multi-time-origin diffusion estimator; Multi-origin MSD is implemented separately in `observables/multi_origin_msd.py` and `observables/multi_origin_pipeline.py`; choose the estimator explicitly. Fit-window sensitivity and equilibrated trajectory selection matter.

### Experimental neuronal-interface extensions

```yaml
analyses:
  local_environment:
    enabled: true
    protein_selection: "protein"
    anchor_selection: "protein and name CA and (resid 151-170 or resid 349-368)"
    lipid_selection: "resname POPC POPE POPI DOPC DPPC PSM CHL1 BGAL GM1a CER160"
    glycan_selection: "resname AGLCNA AIDOA BGLCA BGAL BXYL"
    ion_selection: "resname SOD CLA POT"
    annulus_A: 12.0
    contact_cutoff_A: 5.0
    motifs:
      nxip: "N.IP"
      tm_gxxxg: "G...G"
  coupling_screen:
    enabled: true
    block_frames: 50
    include_defaults: true
statistics:
  hypothesis_tests:
    enabled: false
```

Validate selections for each topology before production. Local membrane enrichment requires an orthorhombic simulation cell in this version. The motif lookup uses the topology sequence and is chain-local; motif contacts are *not* the same as motif SASA or binding capacity. Correlations are descriptive, not causal.

### Performance, data integrity and testing

- The common trajectory pass shares I/O for RMSD, RMSF, Rg, density, original membrane MSD and inexpensive contacts.
- Specialized readers are individually checkpointed by analysis name; result profiles preserve the existing Parquet/NPZ/npy contract.
- Older outputs can be reused when the required files and checkpoint structure are compatible. Use a distinct output directory for an intentional recomputation with changed scientific settings.
- `python -m compileall -q analysis tests` checks syntax; `pytest -q tests` runs synthetic regression checks.
- CI runs these checks on pushes and pull requests. End-to-end tests against a real protein/lipid trajectory are still required before manuscript use.

For usage and prior YAML options, see `README.md`.

---

### Lipid and interface analyses

All additions are opt-in and share MDscope's existing replicate-aware statistics, profile bootstrap and visualization. You can run all selected passes with `python -m mdscope.cli run config.yaml`.

### Issue #1: all-atom C–H order

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

### Issue #2: progress

All shared streaming, lipid-order, leaflet-normal, and enrichment loops report progress. An interactive terminal uses `tqdm` when installed; noninteractive Slurm output prints periodic percentage milestones without requiring an extra dependency. Disable using:

```yaml
progress:
  enabled: false
```

### Issues #3 and #4: normals and local enrichment

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
    footprint_selection: "protein and (resid 146:174 or resid 416:444)"
    shell_edges_A: [0, 5, 10, 20]
    exclusion_radius_A: 3.0
    slab_half_width_A: 8.0
    area_grid_size: 64
    contact_cutoff_A: 5.0
    species_selections:
      POPC: "resname POPC and name P"
      POPE: "resname POPE and name P"
      CHL1: "resname CHL1 and name O3"
```

By default every selected species must be present. For paired system configurations where GM1a/BGAL may be intentionally absent in the control, use `allow_absent_species: true` to skip missing species while still requiring at least one lipid species; missing categories are not silently interpreted as zero enrichment. The lipid selection must contain a representative anchor for each physical molecule. If lipids in the topology consist of multiple bonded residues, designate only the anchoring residue/atom. Use your actual membrane normal axis.

The circular leaflet observers (normals, enrichment and molecular lipid
contacts) require the configured Cartesian membrane normal to be orthogonal
to the lateral lattice vectors. In-plane triclinic shear is supported (for
example, `alpha=beta=90` with arbitrary valid `gamma` for `axis: z`). Shear
through that normal fails explicitly, including when it appears in a later
frame. Fixed initial per-leaflet order profiles use the same geometry check
when their cohorts are assigned. SDC4 maps and
footprint enrichment retain their stricter orthorhombic requirements.

**Normals**: fit a local plane to neighboring anchor atoms within the same leaflet using periodic coordinates and report its angular deviation from the configured global axis in degrees. This measures **local interfacial tilt/undulation**, not a lipid acyl-chain order parameter or a mechanical modulus. Lipids with insufficient neighbors are recorded as missing.

**Footprint enrichment** uses the minimum lateral distance to the edge of
the projected TM and juxtamembrane heavy-atom disks. Each leaflet has its own
disk union, built from protein atoms within the configured slab of its circular
mean anchor height. A single representative lipid anchor assigns its molecule to
a band. The same distance function assigns periodic grid samples to bands and
gives the accessible area of each band and leaflet. Anchors inside the disk union
are excluded from the local bands. The bulk population includes all anchors in the
same leaflet.

The composition fold is the band mole fraction divided by the same-leaflet mole
fraction. The area-density ratio is the species count per accessible band area
divided by its count per accessible leaflet area. Empty or undefined denominators
produce missing ratios. Direct proximity is reported separately as the fraction
of each leaflet/species population within the 3D heavy-atom contact cutoff and
the mean and closest of the per-molecule minimum 3D distances to the selected
protein atoms.
Bonded lipid fragments include linked sugars for the 3D measurement. A bondless
topology is suitable only for single-residue lipids.

- `extensions/bilayer_normal/timeseries.parquet` and optional chunked `per_lipid/angles_part_*.parquet`
- `extensions/lipid_footprint/` contains framewise metrics, counts, selection audit and reference settings

The offline `composition_summary.parquet` inventories the sampled species
population and availability by leaflet and band. Scalar metrics enter the
ordinary replicate and bootstrap statistics. Grid
resolution, disk radius, slab width and lipid anchor chemistry define the model.
Validate their sensitivity before interpreting a small compositional difference.

### Issue #5: lipid–lipid H bonds

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

### Testing and limitations

Tests include circular leaflet assignment across periodic z, an ideal tilted membrane plane, enrichment denominator correctness, carbon-profile selection, pair canonicalization, and progression fallback. Automated CI checks imports and the synthetic streaming trajectory too.

These tests **do not** establish physical validity on your large, heterogeneous neuronal bilayers. In particular: verify anchors, leaflet assignment, coordinate wrapping, headgroup chemistry and acceptor/donor guesses using TPR+XTC before treating output as manuscript-ready. With three independent replicates, interpret confidence intervals and between-replicate variability, not frame-pooled statistical significance.

### Validated Syndecan-4 TPR selections (test-data-v1, October 2026)

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

### Full-length physical audit and performance profiling

`tests/real_trajectory_physics.py` reads **all** frames of both published TPR/XTC test systems for leaflet continuity, local normal fitting, PBC invariance, and lipid enrichment. It explicitly shifts subsets of lipid and protein atoms by entire box vectors and checks that the physical observables do not change. It compares the naive protein center to the circular, periodic-aware center and records how often their enrichment counts differ. A separate 20-frame hydrogen-bond test uses complete GM1a/BGAL/CER160 chains and checks bonded molecule counts and removal of within-molecule H-bonds.

Per-analysis timing is recorded in `<system>/<replicate>/stage_timings.json`. Full-frame audit timings and peak resident memory are in `physical_audit.json` and `<system>/physics/geometry_per_frame.parquet`. These are single-run wall times on GitHub-hosted CPU hardware, not portable throughput benchmarks. They are useful for prioritizing improvements but should be remeasured on the intended HPC system.

**Interpretation:** a 200-frame audit tests input geometry, selection stability and PBC invariance, not independent-replicate statistical significance or convergence. For protein-centered enrichment, prefer a transmembrane-domain selection rather than the whole extracellular protein when the biological question is membrane recruitment. The lipid-order module still assigns leaflet profiles from the first analyzed frame; check the leaflet-change audit before interpreting those profiles on trajectories with actual flip-flop.

---

### SDC4 mechanism analyses: eight opt-in endpoints

The `analyses.sdc4_mechanisms` observer uses the **existing single raw-frame
trajectory pass**. It never aligns or reopens the XTC. The
`analyses.shared_states` project finalizer reads saved per-replicate tables
and fits **one** state model across all conditions and replicates.

**Do not enable these analyses before validating topology selections.**
`chainID`, `resid`, glycolipid anchor names and HS attachment mappings are
specific to the actual TPR and simulation build. The example below uses
the older SDC4 PROA/PROB numbering as an *illustration*, not a guarantee for
rebuilt mature-protein systems.

### Config example

Add this under `analyses:` in `synd4_baseline.yaml`, adjusting all selectors:

```yaml
  sdc4_mechanisms:
    enabled: true
    axis: z
    protein_selection: "protein"
    tm_chains:
      A: "protein and chainID PROA and resid 146:170 and name CA"
      B: "protein and chainID PROB and resid 344:368 and name CA"
    ectodomains:
      A: "protein and chainID PROA and resid 1:145 and name CA"
      B: "protein and chainID PROB and resid 270:343 and name CA"

    # One lipid-anchor *residue* per molecule. For multi-residue GM1a
    # select the ceramide residue, not the separate sugar residues.
    lipid_anchors:
      POPC: "resname POPC and name P"
      POPE: "resname POPE and name P"
      POPI: "resname POPI and name P"
      DOPC: "resname DOPC and name P"
      DPPC: "resname DPPC and name P"
      PSM: "resname PSM and name P"
      CHL1: "resname CHL1 and name O3"
      GM1a: "chainID GM1a and resname CER160 and name NF"
      BGAL: "chainID BGAL and resname CER160 and name NF"
      CER160: "chainID CER160 and resname CER160 and name NF"
    allow_absent_species: true
    # Thickness/curvature use comparable P-group surface anchors only,
    # avoiding artifactual height shifts from cholesterol O3 or ceramide NF.
    surface_anchor_species: [POPC, POPE, POPI, DOPC, DPPC, PSM]
    lipid_contact_selection: >-
      resname POPC POPE POPI DOPC DPPC PSM CHL1
      or chainID GM1a BGAL CER160
    water_oxygen_selection: "resname TIP3 and name OH2"
    ion_selections:
      Na: "resname SOD"
      K: "resname POT"
      Ca: "resname CAL"

    # Explicitly supply a one-to-one HS chain -> protein attachment mapping.
    # Or enable hs_auto above for bond-verified glycan component discovery.
    # Empty maps without hs_auto permit geometry-only testing.
    hs_auto:
      enabled: true
      glycan_selection: "resname AGLCNA AIDOA BGLCA BGAL BXYL and chainID PROA PROB"
      attachment_selection: "protein and resname SER and name OG"
      min_glycan_residues: 3
      expected_chains: 6
    hs_chains: {}
    hs_attachment_sites: {}
    # Example only after checking the TPR topology:
    # hs_chains:
    #   HS_A1: "chainID abA"
    # hs_attachment_sites:
    #   HS_A1: "protein and chainID PROA and resid 39 and name CA"

    motif_selection: "protein and name CA"
    motifs:
      NxIP: "N.IP"
    require_motifs: true
    require_gxxxg: false

    contact_cutoff_A: 5.0
    hydration_cutoff_A: 3.5
    ion_cutoff_A: 4.0
    map_extent_A: 40.0
    map_bin_A: 2.0
    local_radius_A: 20.0
    local_min_lipids: 8

    # Requires explicit lipid C–H hydrogens and bonded all-atom topology.
    # The same -S_CH convention as the existing membrane_order analysis.
    local_order:
      enabled: true
      lipid_selection: "resname POPC POPE POPI DOPC DPPC PSM"

  shared_states:
    enabled: true
    n_states: 3
    seed: 2026
    training_frames_per_replicate: 2000
    bootstrap_draws: 300
    min_complete_fraction: 0.90
    # Default features require ectodomain A/B and a populated local
    # bilayer. You can instead supply the same explicit features in
    # every system; do NOT include species-absent metrics.
    features:
      - sdc4_ectodomain_A_height_A
      - sdc4_ectodomain_B_height_A
      - sdc4_ectodomain_A_tilt_deg
      - sdc4_ectodomain_B_tilt_deg
      - sdc4_tm_crossing_angle_deg
      - sdc4_local_thickness_A
```

**Important:** `allow_absent_species` skips missing lipid species (e.g. GM1a
in a no-headgroup condition). To compare density maps across conditions,
use the same species list and report absent species as biologically absent,
not as zero-density enrichment. Density map coordinates are rotated so the
x-axis points from TM chain A to chain B. Leaflet-resolved density is in
lipids/nm² and fold enrichment uses each leaflet's mean species density.
Maps are **not** corrected for the protein-excluded accessible area.

### Observable definitions and scientific limitations

| Endpoint | Output/interpretation |
|---|---|
| Lipid enrichment and density maps | Per-leaflet GM1a/BGAL/cholesterol/phospholipid maps centered on TM midpoint, plus the existing annular `lipid_enrichment` fold values when enabled |
| Ectodomain height/tilt/contact | CA-center height relative to extracellular leaflet, terminal-to-terminal tilt to the outward membrane normal, binary membrane-contact state and fraction of ectodomain atoms contacting selected lipid atoms |
| NxIP accessibility | Regex matches actual CA-derived topology sequence by chain; reports fraction of motif heavy atoms with protein/HS/lipid steric neighbors within cutoff. This is an **occlusion proxy**, not a ligand-binding free energy. When `sasa.enabled`, the observer also joins the existing FreeSASA motif residue SASA, calculated with that analysis's atom selection |
| Individual HS orientation/contacts | Attachment-to-distal-terminus orientation cosine, extension, direct HS–lipid contacts and separate attachment-site contacts; distal end chosen by attachment proximity on the first analyzed frame |
| TM crossing and GXXXG | PCA helix-axis crossing angle, handedness and cross-chain GXXXG heavy-atom/residue contacts; no GXXXG motif produces zero contacts unless `require_gxxxg` is true |
| HS hydration/ions | Unique water oxygens within hydration cutoff per HS heavy atom and unique ions within coordination cutoff; contact-conditioned within-replicate means. These are coordination proxies, **not** electrostatic free energies |
| Local membrane structure | Local leaflet separation, per-leaflet quadratic height-field curvature approximation, and optional local minus-S_CH averaged equally per lipid |
| Shared conformational states | Replicate-balanced standardized k-means across all conditions, per-replicate state populations, replicate-bootstrap system means, within-replicate moving-block intervals, transitions, and first/second-half JS divergence |

When interpreting motif SASA, note that a protein-only SASA calculation does not
include HS/lipid shielding; the steric occlusion metric reports these
separately. Adjacent sequence neighbors (±2 residues) are excluded from
protein occlusion to avoid trivial backbone contacts.

The ectodomain terminal direction assumes the supplied CA selection is in
**N-to-C order**, with the N terminus distal to the membrane. For branched
HS chains, the first/last-residue terminal approximation may be unsuitable:
inspect the `selection_audit.json` and supply a more appropriate selection
or a dedicated branched-glycan end-to-end definition before interpreting
orientation.

The local membrane model requires an **orthorhombic, z-normal** box. It
uses only `surface_anchor_species` for thickness/curvature to avoid mixing
chemically different headgroup heights. It fits a quadratic height field and reports a small-slope curvature proxy in
Å⁻¹; it is not a full triangulated-surface mean-curvature calculation.
For local order, tail sites use explicit bonded C–H pairs and the same
minus-S_CH convention as the original lipid order analysis. Curvature and
thickness require at least eight lipid anchors per local leaflet.

**State modeling:** training uses equal numbers of frames from each
replicate, one scaler, and one model across all conditions. It does not
claim that different states imply different thermodynamic free energies.
Frame bootstrap intervals are **not** independent-replicate confidence
intervals. For biological condition comparisons, independent simulations
are the unit of replication; first/second-half divergence is a sampling
diagnostic, not proof of convergence.

### Outputs

Per replicate:

- `extensions/sdc4_mechanisms/timeseries.parquet`: all frame-level scalar
  geometry, contact, hydration and membrane endpoints
- `extensions/sdc4_mechanisms/lipid_density_maps.npz`: per-leaflet 2D
  counts, density and fold maps, grid edges and frame count
- `extensions/sdc4_mechanisms/selection_audit.json`: exact matched
  motifs, GXXXG residue indices, HS termini and lipid counts
- `extensions/sdc4_mechanisms/hs_conditioned_solvation.parquet`:
  HS hydration and ion coordination conditioned on direct lipid contact
- `extensions/sdc4_mechanisms/lipid_{density,enrichment}_maps.png`
- `extensions/sdc4_mechanisms/timings.json`

Project-level (when `shared_states.enabled`):

- `project_summary/shared_states/model.json`
- `state_centroids.parquet`, `frame_assignments.parquet`
- `replicate_populations.parquet`, `system_populations.parquet`
- `sampling_diagnostics.parquet`, `state_conditioned_hs_solvation.parquet`
- `state_populations.png`, `state_sampling_diagnostics.png`

The scalar `timeseries.parquet` columns are also merged into MDscope's
existing `summary/timeseries.parquet` and analyzed by the existing
replicate-aware statistics framework. No new trajectory pass is added.


### Automatic covalent HS chain discovery

In the original carb/no-carb TPRs, each SDC4 protomer has three
covalently attached HS chains but the sugar residues share the same
`chainID` as the core protein. MDscope's `hs_auto` uses the TPR bond
graph to separate glycan-only connected components and requires exactly
one covalent protein attachment per glycan. The two audited topologies
have attachments PROA S39/S61/S63 and PROB S309/S331/S333. The
`selection_audit.json` includes per-chain protein attachment details.
The auto-discovery refuses missing or ambiguous attachments; it does
not infer attachments from spatial proximity.

See `synd4_full_mechanisms.yaml` for the complete configuration and
the [execution guide](workflows.md) for findings review.

---

### Opt-in multi-time-origin lipid diffusion

Enable in `analyses.membrane.msd`:

```yaml
msd:
  enabled: true
  multi_origin: true
  multi_origin_max_lag_frames: 1000  # null = all sampled lags
  multi_origin_particle_batch: 64
  slope_start_fraction: 0.2
  slope_end_fraction: 0.8
```

Runs a **separate trajectory pass** per replicate after the standard streaming
stage. This is deliberately opt-in and preserves the existing single-origin
timeseries and diffusion estimates. The optional `msd_trajectory` is used if
specified; otherwise the main trajectory is used. Neither can be assumed
unmodified merely from filenames.

For each selected lipid, the reader accumulates minimum-image frame-to-frame
center displacements (and optionally subtracts common drift), then writes
these into a temporary on-disk float64 NumPy memmap. The FFT MSD is computed
over small batches of lipids and averaged across lipids and all available time
origins. Disk usage is approximately
`24 * sampled_frames * selected_lipids` bytes, plus temporary FFT memory
bounded by `multi_origin_particle_batch` (the final tables are much smaller).
Temporary coordinates are deleted when processing finishes or errors.

The maximum lag is in **sampled frames** after trajectory stride. The fitting
window is expressed as fractions of maximum retained lag time, excluding lag
zero, and must contain at least four lag points. Sampling must be regular.
Output units are Å² for MSD and nm²/ns for D, using the 2D Einstein relation
`D = slope/4`. Reported D values are descriptive slopes, not uncertainty
estimates; independent replicates remain the inference units.

Replicate outputs:
- `summary/multi_origin_msd.parquet`: lag time, MSD, number of time origins
- `summary/multi_origin_diffusion.parquet`: per-replicate fit by all/upper/lower lipids

Project outputs:
- `project_summary/diffusion/multi_origin_replicates.parquet`
- `project_summary/diffusion/multi_origin_systems.parquet`

Set `analyses.membrane.msd.multi_origin_fit_diagnostics: true` to add
`project_summary/diffusion/multi_origin_fit_diagnostics.parquet` during
finalization (all or plotting phase). It reads saved curves and works with
existing replicas offline. It compares the configured fit interval with
its earlier and later overlapping two-thirds. Each fit needs four positive
lag points and at least 10 origins or 10% of the sampled frames, whichever is
larger, at every included lag. These are overlapping origins, so the count is
only a support guard. The table reports observed lag bounds and origin counts,
D, linear fit R², and the log–log MSD slope by replica and lipid cohort.

The descriptive screen flags a spread in D greater than 25% of the configured
window's D or a log–log slope more than 0.20 from unity in any window.
These are pragmatic sensitivity thresholds, with no inferential calibration.
Missing curves and insufficient or nonpositive fits are marked unavailable.
A result within the screen cannot establish asymptotic diffusion, rule out
drift/finite-size effects, or supply a confidence interval. Where support is
limited or windows disagree, inspect longer trajectories before interpreting D.

**Limits:** Minimum-image reconstruction assumes sampled displacements are less
than half the relevant box dimension in an unambiguous image sense. The FFT
kernel assumes stable lipid identity, a fixed analysis frame, and that the
appropriate raw trajectory has been supplied. Long-lag MSD estimates have
fewer time origins and correlated lag values; high R² is not evidence that the
diffusive regime was reached. The second pass increases trajectory I/O but
does not alter standard streaming checkpoints.

---

### Fused solvation analysis (optional)

The solvation observer shares MDscope's raw-frame fused processing. It runs one
periodic capped-distance query per solvent class and frame, writes atomic chunks,
and reuses them for offline counts, full-cell-density-normalized radial distributions and
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

The RDF uses a **full periodic cell average** number density: total
selected waters/ions divided by periodic cell volume. It does not measure
a separate aqueous bulk reference; `g_bulk` is retained as a legacy column name. This is **not corrected**
for the solvent-excluded volume of protein, glycans and bilayer. Buried sites
can exhibit lower radial density solely because water is inaccessible.
Interpret RDF alongside SASA and occlusion, and calibrate first-shell
boundaries using stable first minima and cutoff sensitivity analysis.
Shell cutoffs are explicit and are *not automatically fitted* in this version.
The pooled site-class RDF must not be mistaken for per-residue density.

Current limitations include no solvent-excluded-volume normalization, no full
intermittent contact correlations, and no general water-acceptor bridge detector.

### Postprocessing and optional 3D water density

After the fused pass, `shell_candidates.parquet` records a candidate first
minimum of each full-cell-density-normalized RDF with an explicit quality status. It does
not automatically change already computed coordination cutoffs: calibrate on
representative frames and fix the same chemical-class cutoff across compared
conditions before conducting kinetics. Ambiguous RDFs yield missing candidates.

`target_residues.parquet` maps target IDs to topology residue indices. When
a matching SASA run exists, `hydration_sasa.parquet` combines the **saved**
per-residue SASA with per-target solvent counts using exact frame time; missing
constituent residues produce NaN rather than partial SASA. Set `sasa_run` to
the selected run ID if needed. SASA exposure is a contextual covariate, not an
automatic correction to full-cell-density-normalized RDFs.

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

The numeric ID and atom indices above are placeholders. To specify atom
order explicitly, use a list of three selectors that each match one atom,
in origin, axis, and plane order. Maps use the same
periodic candidate-water pairs as the radial calculations and accumulate
co-moving local density in `water_density_3d.parquet`. Only the inscribed
sphere (radius half the grid width) is sampled, and the frame origin must be
part of the target site selection. Missing voxels *inside* this sphere have
zero observed density; voxels outside are undefined, not zero. By design, this
is a sparse observed-density map: absent voxels correspond to zero observed
counts. Its legacy-named `g_bulk` uses full periodic cell average density and does not adjust for inaccessible
local volume. A local map cannot be interpreted as an equilibrium fixed-frame
binding site if its neighboring molecular environment varies substantially.

### Reliability and performance validation

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

### Optional geometry-qualified water bridges

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


### Molecular aggregation and ion species

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

### Bounded identity tracking

The standard preset computes RDFs and per-target counts for **all** configured
targets without writing solvent identities. Detailed/kinetics mode can limit
expensive solvent identities and bridge enumeration to selected classes:

```yaml
    preset: kinetics
    options:
      track_target_classes: [protein, hs]
```

This retains water/ion coordination counts and RDFs for lipids, but
persists interaction identities and bridges only when their constituent
targets are among the specified classes. An omitted filter tracks all
classes (potentially very expensive); avoid full-system water tracking
without first benchmarking a representative frame window.

---

### Solvation analysis for the original syndecan-4 CHARMM-GUI models

The repository's `synd4_full_mechanisms.yaml` uses old-build topology paths
on Tinkercliffs. **Those files are not present in this GitHub repository.**
Consequently, atom selections below are candidates inferred from the existing
configuration, not independently validated on the TPRs. Do not enable a full
multi-microsecond run until the preflight succeeds for **both conditions**.

### Minimal common-target configuration

Add this under `analyses:` in a copy of the existing YAML. It is deliberately
disabled, and no costly kinetic preset is enabled by default.

```yaml
  solvation:
    enabled: false
    preset: standard
    water_selection: "resname TIP3 and name OH2"
    ion_selections:
      K: "resname POT"
      Cl: "resname CLA"
    targets:
      protein:
        acidic_sidechain:
          selection: "protein and ((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))"
          group_by: residue
        basic_sidechain:
          selection: "protein and ((resname LYS and name NZ) or (resname ARG and name NE NH1 NH2))"
          group_by: residue
      hs:
        glycan_oxygen:
          selection: "(resname AGLCNA AIDOA BGLCA BGAL BXYL) and (chainID PROA PROB) and name O*"
          group_by: hs_chain
      lipid:
        phospholipid_phosphate:
          selection: "resname POPC POPE POPI POPS DPPC DOPC PSM and name O11 O12 O13 O14"
          group_by: residue
    search_radius_A: 8.0
    rdf_bin_width_A: 0.2
    shell_cutoffs_A: {water: 3.5, ions: 4.0}
```

The `hs_chain` grouping reuses MDscope's bond-based
`discover_covalent_hs` and, when present, the existing
`sdc4_mechanisms.hs_auto` settings. It identifies six attached HS chains
as separate glycan-only bonded components (rather than merging them through
the shared PROA/PROB protein backbone). This works only if glycan atom
selections cover the **complete glycan bond graph** and bond connectivity
is present in the TPR/PSF. The example `name O*` defines measured hydration
sites, not the glycan component discovery selection.

The LYS/ARG selections above are appropriate for *water proximity*, but
ion coordination requires **chemical interpretation by ion sign**.
Cl⁻ near the protein basic groups and K⁺ near acidic groups are not
interchangeable. Sulfate-specific versus generic glycan oxygen groups should
be defined using atom names audited from the actual HS topology.
Likewise, this example's phosphate oxygen names must be checked against
the actual CHARMM lipid topology.

GM1a and BGAL headgroups are missing in the nonglycolipid condition.
Define separate **condition-specific** solvation target selections for those
species or analyze them individually. Do not silently treat a missing
species as having zero hydration. GM1a's glycan residues can span multiple
topology residues; do not assume one residue represents an entire lipid.

### Preflight and real observer benchmark on Tinkercliffs

Run the tool from an editable MDscope installation or with its package on
`PYTHONPATH`. Use a new YAML file whose `solvation.enabled` is true.
You may temporarily switch off other analyses in that YAML to avoid
unrelated validation requirements.

```bash
python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --output solvation_benchmark_standard.json

python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset detailed --output solvation_benchmark_detailed.json
```

Repeat for `synd4_no_glycolipids`. Compare walltime, peak RSS,
per-frame processing rate, and output sizes. The benchmark invokes the
actual solvation observer (including chunk flush and finalization), but
not the rest of MDscope's complete analysis pipeline. Memory numbers
include the loaded Universe and dependencies.

**Scientific limitations:** RDF normalization uses selected-water count / full periodic cell volume, not
the solvent-accessible volume around a site. Analyze saved SASA with water
counts and RDFs, and evaluate cutoff stability across trajectory chunks
and independent replicates. Unresolved RDF minima should not be accepted
as automatic shell cutoffs. With an XTC saved every 100 ps and a trajectory
stride of 10, adjacent analyzed frames would be 1 ns apart; water-contact
episodes then quantify sampled persistence, not picosecond water residence.
Check actual `ts.time` rather than assume the interval.

### Validation criteria

1. The preflight outputs the expected topology atom names, nonzero target
   counts and per-solvent counts for each condition; no unexpected overlaps.
2. Run with standard mode and inspect RDFs and available SASA for burial.
3. Check candidate minimum stability for each chemical class, with matched
   cutoff criteria between conditions.
4. Only enable detailed or kinetics modes on selected target sites once
   output-size and physical frame-interval limitations are understood.
5. Compare per-replicate results; do not use individual frames or hydration
   episodes as independent biological replicates.

---

### Solvation statistics: reuse existing moving-block and replicate hierarchy

The per-replicate solvation observer already writes:
- `extensions/solvation/counts.parquet`: per-frame, per-target unique solvent counts with solvent species, target IDs/classes and time.
- `extensions/solvation/species_timeseries.parquet`: per-frame chemical-class average coordination and count; now also `hydrated_fraction`, the fraction of configured targets with ≥1 solvent entity within its coordination shell.
- `extensions/solvation/target_summary.parquet`: per-residue or per-molecule mean, standard deviation and hydrated-frame fraction, descriptive **within a replicate**.
- `extensions/solvation/rdf.parquet`: full-cell-number-density-normalized radial distribution by solvent and target/site class, with saved radial bins and raw pair counts.
- `extensions/solvation/{events,bridges,water_hbond_bridges,episodes,survival}.parquet` when the relevant opt-in modes are enabled; other outputs include `shell_candidates`, `shell_stability`, `target_residues`, optional `hydration_sasa`, and optional sparse `water_density_3d`.

### What has changed

`statistics.time_series.collect_time_series` now discovers **existing** solvation
`species_timeseries.parquet` alongside other MDscope time series. No new
resampling method or separate statistical implementation is added.
Each stable key is `solvation/{solvent}/{target_class}/{site_class}`,
with metrics `mean_n_solvent` and `hydrated_fraction`, where available.

With `statistics.enabled: true`, the **existing** framework calculates
FFT-based autocorrelation, correlation-informed temporal block lengths,
within-replicate moving-block bootstrap, between-replicate hierarchical
bootstrap, replicate means/CI, and difference-of-conditions CI.
Existing permutation/FDR settings apply; **with only 3 independent
replicates per condition, do not overinterpret p-values**.

Output tables, in the existing project-wide structure:
- `project_summary/statistics/time_series/replicate_statistics.parquet`
- `project_summary/statistics/time_series/system_statistics.parquet`
- `project_summary/statistics/time_series/system_comparisons.parquet`
- `project_summary/statistics/time_series/autocorrelation.parquet`
- `<replicate>/summary/{time_series_statistics,autocorrelation}.parquet`

All solvation rows are identifiable with `group` beginning `solvation/`.
Comparison difference convention is **system_2 minus system_1**. Use named
system columns and the actual direction; do not assume a sign convention.

RDFs are **curves**, not independent scalar samples. We therefore also output:
- `project_summary/statistics/solvation/replicate_rdf.parquet`: raw per-replicate curves.
- `project_summary/statistics/solvation/system_rdf.parquet`: pointwise
  replicate mean and SD, number of contributing replicates, and coverage flag.
  This is **descriptive** replicate variation, not per-bin confidence intervals
  or hypothesis tests. The RDF first-minimum cutoffs are exploratory.

### Important limitations

1. Each simulation replicate is the independent unit. Frames and multiple lipid
   molecules in the same trajectory are **not** independent replicates.
2. The class mean weights each configured molecular target equally at every
   frame. Compare the same chemically meaningful target definitions across
   conditions; `target_count` is available to audit accidental composition changes.
3. `hydrated_fraction` indicates ≥1 selected solvent inside the fixed shell.
   It is **not** the probability of a specific water identity or an
   uninterrupted residence event. An ion-specific fraction uses that ion
   species' configured selection and cutoff.
4. Residence/survival and 3D density remain descriptive by default.
   We do not bootstrap censored residence curves or every voxel.
5. A preexisting `species_timeseries.parquet` lacking `hydrated_fraction`
   still supports bootstrap of `mean_n_solvent`; regenerate the saved class
   time series from raw counts if occupancy statistics are also required.
6. Using many chemically fine target classes increases the number of tests
   and bootstrap work. Default to scientifically relevant classes/regions,
   not every target atom.

No additional trajectory pass or extra neighbor search is necessary.

---

### Solvation hot-path audit and benchmark (2026-11-01)

This reviews the existing merged PR #78 implementation, rather than
introducing an additional analysis framework.

### Code footprint

PR #78: **1,608 added lines in 11 files** (1 deletion), including 298 lines
of documentation, 196 lines of benchmark scripts, and 355 lines of tests.
The numerical/postprocessing implementation was 505 + 203 = 708 lines;
the remaining lines were configuration and fused-pipeline wiring. Removing
documentation or independent tests would reduce line count without eliminating
runtime overhead. No new abstraction layer is justified.

### Locally measured synthetic benchmark

Measured with MDAnalysis 2.10.0 and NumPy, 150 Å cubic periodic box,
2,400 selected atoms grouped into 800 targets, 80,000 water oxygens,
8 Å RDF search and 3.5 Å coordination shell. Uniform random coordinates;
**not** a biomolecular trajectory or a validation of physical hydration.

| Operation | Baseline | Optimized | Notes |
| --- | ---: | ---: | --- |
| Periodic search, 8 Å | 85.97 ms | unchanged | Median over 5 invocations; 122,166 pairs |
| Deduplicate + count, 3.5 Å | 1.61 ms | 1.08 ms | Median over 150 iterations; 10,276 shell pairs |
| Neighbor search backend auto vs nsgrid | 123.1 ms | 98.9 ms | Different short benchmark, 4 repetitions; **not a verified production speedup** |
| Neighbor search backend pkdtree | 277.2 ms | — | Same short benchmark, 4 repetitions |

The optimized coordination counts match the original Python set algorithm
exactly, including repeated site–water pairs and unoccupied target IDs.
**The full workload gain will be small** because the periodic neighbor
search dominates this synthetic scenario. The code change is valuable because
it avoids per-pair Python objects entirely in the standard preset. It should
not be sold as a large end-to-end speedup.

### Code changes

- `unique_target_solvent_pairs` performs lexicographic sort, adjacent-pair
  deduplication and per-target bincount, retaining exact target/water identity
  when detailed analysis needs it.
- Standard mode **does not create** the large Python `dict[int, set[int]]`
  that detailed/kinetics mode requires for solvent identity output.
- Optional `neighbor_method` defaults to `auto` to preserve existing
  behavior. Alternatives are `nsgrid`, `pkdtree`, `bruteforce`.
- The existing *real* trajectory benchmark takes `--neighbor-method`;
  compare options before setting any HPC-wide default.

### HPC benchmark to perform before tuning further

On Tinkercliffs, with validated solvation target selections, run the same
representative trajectory window with each neighbor backend:

```bash
python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --neighbor-method auto \
  --output solvation-auto.json

python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --neighbor-method nsgrid \
  --output solvation-nsgrid.json
```

Repeat with `--preset kinetics` and an explicit
`options.track_target_classes` for interesting protein/HS sites, then
repeat in the nonglycolipid condition. Compare **total observer throughput**,
maximum RSS and bytes of output, not just the native PBC search speed.

### Remaining performance risks

1. Full-system water identity tracking and all-target bridge enumeration
   can still produce enormous tables; use the built-in class filter.
2. `finalize()` concatenates all identity-event chunks in memory before
   writing outputs, which may become limiting on very long trajectories.
   A subsequent bounded-memory streaming finalizer is worthwhile **only
   after measuring the affected realistic kinetics preset**.
3. The radial search radius determines a large fraction of candidate-pair
   work. Preserve a radius sufficient for RDF structure/3D targets;
   do not shorten it purely to improve benchmark speed.
4. The current benchmark does not profile other fused MDscope observers
   (SASA, DSSP, protein/lipid mechanisms) running concurrently.
5. Profile a representative multi-thousand-frame run for queue-time
   extrapolation and space usage before comparing complete 2 µs replicates.

**No production TPR/XTC was accessible in this environment; performance
numbers are synthetic and do not establish real MDscope throughput.**

### Shared PBC implementation

The review also found duplicate calls to MDAnalysis
`minimize_vectors` inside the solvation H-bond kernel and local 3D
water-density alignment. Both now use MDscope's existing
`core.geometry.minimum_image_displacements`, which preserves a single
triclinic minimum-image validation contract. A dedicated triclinic wrapped
H-bond regression was added. The sparse neighbor search still calls
MDAnalysis `capped_distance` directly because its atom-pair + distance
output is not equivalent to the existing `biophysical_common.atom_contacts`
(pair-only API); forcing one API onto the other would add complexity and
unnecessary distance computation for existing callers.

RDF histogram values at *exact* bin boundaries remain sensitive to
floating-point rounding between native distance backends. This is a
histogram numerical-convention issue, not a distinct PBC convention;
ordinary nonboundary backend parity is tested separately.


## Molecular lipid neighbor matrices and bridge episodes

Optional `analyses.lipid_contacts` shares the raw fused membrane pass. Configure
`enabled: true`, `selection` (all heavy atoms of the desired complete molecules),
`species_selections` (one anchor residue per molecule, as for lipid enrichment),
`cutoff_A` (default 5), `axis` (default z), and `batch_atoms` (default 4096).
Bond connectivity is required. Every selected molecule must have exactly one
anchor; the mapping is validated rather than guessing glycolipid residue IDs.
This uses the circular planar leaflet policy described above, including
supported in-plane lattice shear.

Each distinct bonded molecule pair is counted once irrespective of atom count,
with self pairs removed and same-species pairs retained. AB and BA rows repeat
the same undirected pair count. Neighbors per A molecule are C_AB/N_A for A!=B
and 2C_AA/N_A for A==B. Possible pairs are N_A*N_B or N_A*(N_A-1)/2. Enrichment
is observed contact share divided by possible-pair share among all distinct
molecules in that instantaneous leaflet. This finite-size random-mixing null
has enrichment 1 for a complete contact graph; it does not correct for lipid
size, steric packing or accessible membrane area and is not an affinity.
Absent populations and undefined reference populations produce NaN, while
an available pair class with zero contacts produces zero enrichment when the
leaflet has other contacts. Cross-leaflet pairs are counted separately in the
scalar output and excluded from within-leaflet normalization.

Outputs are `extensions/lipid_contacts/matrix.parquet`, `timeseries.parquet`
and `selection_audit.json`. Searches batch reference atoms and retain only
unique encoded molecule pairs between batches, not atom event tables. Chunk
and checkpoint numbers follow the existing fused observer. Matrix group time
series use the existing temporal/hierarchical statistics and standard plots;
condition/leaflet heatmaps average replicas equally and are descriptive.
Difference heatmaps subtract condition means in configuration order and show
a cell only when both species and a finite estimate are present in every
configured replica on each side. The companion replicate CSV includes mean
undirected pair counts and mean leaflet species populations. Raw matrices
and audits remain available as saved products.
The optional module is off by default; enabling it adds a necessary neighbor
query, not a second trajectory pass.

Solvation sequence-defined targets now accept `atom_selection` within the
mapped motif. For the acidic sequence, `atom_selection: "name OD1 OD2 OE1 OE2"`
restricts sites to explicitly named carboxylate atoms after topology mapping.
HS carboxylate/sulfate selections must likewise be explicit and validated for
the force field. Do not substitute whole-motif proximity for ion coordination.

When sparse `bridges.parquet` records exist, offline derivatives add
`bridge_residence.parquet` and `bridge_residence_intermittent.parquet` without
rereading trajectories. A bridge requires the **same** solvent entity to
coordinate both named targets simultaneously. Duplicate records are removed;
continuous/intermittent episodes share the validated sampled-duration kernel,
retain ion identities and both censor flags, and report observed spans as
well as interval-imputed durations. The intermittent gap is one sampled frame,
not a fixed physical duration. Only observed target pairs enter these outputs;
an absent pair is not automatically a measured zero. These geometric episodes
do not establish cooperative binding, contact-ion versus solvent-separated
mechanisms, or free energies. Null controls and replica-level production
comparisons remain necessary to validate these outputs on the full production runs.

HS outputs additionally report equal-heavy-atom Rg, shape anisotropy (0 for
isotropic, 1 for linear), vertical span/extension and attachment position
relative to the TM center. Coordinates are reconstructed by minimum image
around the attachment, assuming the chain fits in that periodic image. These
reuse current coordinates and do not fit additional conformational states.
State robustness, matched windows and survival comparisons remain research
validation requirements for #104.

Batched-search benchmark (`python scripts/benchmark_lipid_contacts.py` in fresh
processes): 1,024 molecules / 12,288 atoms / three frames, batch 256 took
0.750 s at 170.5 MiB peak RSS; batch 4096 took 0.619 s at 174.7 MiB. Matrix
digests were identical. This demonstrates a memory/runtime tradeoff and
batch invariance, not a production trajectory speedup prediction.

### Optional aqueous-region solvation RDF reference

`analyses.solvation.aqueous_reference` adds `g_aqueous_reference` to the
existing `rdf.parquet`; it never changes the historical full-cell `g_bulk`.
For an axis-aligned orthorhombic bilayer, choose a nonempty membrane anchor
selection and an explicit excluded slab half-width that includes the membrane
and a justified interfacial buffer. For example (the width is illustrative,
not a recommended physical boundary):

```yaml
analyses:
  solvation:
    aqueous_reference:
      membrane_selection: "resname POPC DOPC and name P"
      axis: z
      half_width_A: 30.0
```

The circular membrane center is recalculated per frame. For normal length L,
cell volume V and half-width h, the reference geometric volume is
`V * (1 - 2*h/L)`. Each solvent species uses **the same selected entities**
as its RDF, counted outside that slab. The per-frame count, volume, density
(Å⁻³), center (Å), and availability are in `aqueous_reference.parquet`.
The RDF denominator sums site-atom count times frame-specific density,
rather than using a ratio of average volumes/counts. An empty reference in
any frame makes that species' entire aqueous-reference curve unavailable;
its full-cell curve remains available. Invalid boxes/regions and RDF radii
larger than half the shortest box length fail explicitly. General triclinic,
rotated and strongly curved membranes are unsupported for this option.

This is a **geometric aqueous-region reference**, not proof of bulk density.
Protein/HS volume within that region and the inaccessible part of each RDF
shell remain uncorrected. A hemisphere can still give g(r)=0.5 at the correct
aqueous density. Inspect spatial density/plateau and slab-width sensitivity
on real trajectories before comparing systems. The optional reference is
carried to replica/system descriptive RDF tables with its own replica count
and coverage flag, and shown as dashed curves. Legacy curves never acquire
an inferred reference. No production benefit is claimed from the toy tests.

### Bridge evidence and streamed sampled-contact episodes

Residence finalization reads frame-ordered contact identities in bounded Parquet
batches. It retains active contacts across two successive observation frames
and a bounded output batch; it does not concatenate the entire event history.
Observed-span `duration_ps`, interval-imputed `sampled_duration_ps`, and boundary
censoring retain their definitions. Persisted episode row order follows episode
completion, rather than first-seen entity grouping. Consumers should use identity
and frame columns, not row position. Invalid observation/event order fails before
replacing the previous episode file. Kaplan–Meier output uses aggregated tied
durations and reverse cumulative risk counts, preserving its numerical contract.
Survival, bridge residence and SASA joins have separate table-loading limits;
this does **not** make the entire finalizer bounded-memory.

For explicitly configured `bridge_target_pairs`, finalization also writes:

- `bridge_pair_timeseries.parquet`: one row per observed frame, selected solvent
  species and configured pair, including zeros. Entity counts and presence use
  the **same entity ID** at both targets, not two different nearby ions.
- `bridge_pair_shift_reference.parquet`: at most 100 deterministic, nonzero
  circular offsets of target-B contact histories. All entities share the same
  offset, preserving their identities and each target's occupancy pattern.
- `bridge_pair_control_summary.parquet`: observed occupancy, shift-reference
  mean, 2.5/97.5 percentiles, observed-minus-reference mean and availability.

The two percentiles describe the selected shifts; they are not confidence
intervals or adjusted significance tests. No reference is calculated for fewer
than four observations or irregular times. Circular wrapping and nonstationary
trajectories can invalidate this temporal reference. It breaks simultaneous
alignment, **not spatial crowding**, and cannot establish chemical bridging or
cooperativity. Chemical site identity, shell thresholds and spatial/null geometry
still need validation on the actual system.

The frame traces enter the existing scalar block/replicate/hierarchical
statistics; independent trajectories remain the condition-level units. Group
IDs encode both targets and solvent explicitly. Existing legacy pair summaries
and ion-count field names remain readable, while new traces use neutral entity
labels that also apply to water. Bridge-record counts are checked against the
independent reconstruction from persisted target-contact events.

`target_atoms.parquet` records selected atom names, topology residue membership,
elements and topology charges (missing charges remain NaN). The metadata's
`solvent_selection_audit` records the actual selector, entity count, atom names
and residue names. These records expose the chemical selections; they do not
guess protonation or certify chemistry. Compatible older outputs without these
new optional products remain readable and acquire no inferred null evidence.

### Explicit eligible bridge-pair occupancy

For `solvation.options.bridge_target_pairs`, supply a list of pairs of exact
`target_id` strings from `target_residues.parquet` (or the selection audit).
Both targets must be distinct, tracked, and have bridges enabled. For example:

```yaml
options:
  bridges: true
  bridge_target_pairs:
    - ["protein/carboxylate/7", "hs/sulfate/200"]
```

These are illustrative topology IDs; use the actual saved IDs. The optional
`bridge_pair_summary.parquet` includes every configured pair for each selected
solvent species, even if it never bridged. `sampled_bridge_fraction` is the
fraction of **all sampled observation frames** with at least one *same ion*
simultaneously coordinating both targets. `mean_simultaneous_bridging_ions`
is the unique-ion count summed over frames divided by that same denominator.
Duplicate event records do not increase either metric. An unknown/untracked
pair fails rather than becoming zero. Occupancy is sampled-frame based,
not time-weighted; it does not establish a chemical bridge or cooperativity.
Existing per-ion episode tables and censoring conventions remain distinct.

Set `solvation.options.bridge_spatial_ion: ions/K` to additionally evaluate an
in-plane spatial reference for that explicitly selected species. For each
sampled frame and eligible target pair, eight fixed nonzero fractional cell
translations shift all selected ions together. Periodic wrapping preserves ion
height and the ion–ion arrangement. The reference uses the same contact cutoff
and same-ion requirement as the observations. It uses orthorhombic boxes and
the aqueous-reference normal if configured, with `z` as the default.

`bridge_spatial_reference.parquet` contains the per-offset unique-ion counts.
`bridge_pair_timeseries.parquet` includes the fraction of offsets with a bridge
and the observed-minus-reference contrast. The scalar statistics family can
summarize that contrast across frames and replicas. The translations may put
ions inside excluded volume, so treat this as a descriptive geometry and
crowding reference. Chemical assignments, RDF peak boundaries and sensitivity
to steric exclusion still require separate evidence.

### Per-residue water coordination

Set `residue_occlusion: true` and `residue_water_cutoff_A: 3.5` under
`sdc4_mechanisms`. The observer counts distinct water oxygens from
`water_oxygen_selection` within the cutoff of any protein N or O in each
residue. It stores the water count and the number of eligible polar atoms
alongside framewise protein, HS and lipid contact fractions in
`extensions/sdc4_mechanisms/residue_occlusion.parquet`.

Profile statistics report overall coordination, the fraction of sampled frames
in each partner state, and mean coordination during protein, HS, lipid and
no-partner contact states. The first three states can overlap. Conditional
confidence intervals resample the original time grid and require observations
in multiple blocks. Residues without a selected polar atom have no hydration
estimate. Protein contact states exclude same-chain neighbors within two
sequence positions. The count describes a geometric hydration shell and depends on the
number and arrangement of polar atoms in a residue.


### Potassium bridge audit for the Syndecan-4 release

The production YAML enables `analyses.potassium_bridge_audit`. MDscope runs
it once per replicate in a separate pass on the configured frame window and
time axis. It publishes the completed output directory at
`extensions/potassium_bridge_audit/` and includes the tables in the review
review. Finalization reads the saved products.

The standalone script also accepts a single topology and trajectory. The
example uses the glycolipid release pair. A signed-in GitHub CLI
can download the private release with `gh release download test-data-v1
--repo MarionQL/MDscope --dir real_test_data`.

```bash
python scripts/audit_potassium_bridges.py \
  --tpr real_test_data/carb-sdc4-step7_production.tpr \
  --xtc real_test_data/carb-sdc4-rep1_test.xtc \
  --stride 10 --output potassium-audit-carb-rep1
```

Omit `--xtc` for the topology audit alone. Use `--initial-frame` with
MDAnalysis 2.11 or newer to check the starting structure stored in the TPR,
including its box. A starting structure supports a geometry check only.

- `chemical_sites.parquet` records the terminal oxygens of Asp/Glu
  carboxylates and covalently attached HS carboxylates, sulfates and
  sulfamates, as present. The graph and charges come from the TPR. HS chains
  are identified by attachment bonds. The two exact acidic-patch matches and
  the unbonded +1 CHARMM POT selection are checked before analysis.
- `rdf.parquet` stores 0.1 Å radial bins, raw oxygen–K pair counts and
  full-cell-density-normalized RDFs by oxygen chemistry. The normalization
  includes membrane and solute volume. Compare peak locations across sampling
  windows and the 0.1/0.2 Å smoothing widths in
  `rdf_peak_candidates.parquet`. The first peak, intervening minimum and
  second peak are distance candidates. To evaluate solvent sharing, compare
  the shared-water pair counts against the total counts at each peak. The same
  water must be within 3.5 Å of K and donate an H bond to the target oxygen
  (donor–acceptor distance ≤3.5 Å and angle ≥150°). A pair can satisfy this
  geometric water criterion while also in direct contact.
- `bridge_frames.parquet` and `bridge_summary.parquet` keep each patch,
  HS chain and HS chemical class separate. Compare the observed fraction
  of frames with a same-K bridge with the fraction of frame–offset samples
  in the shifted reference, at each contact cutoff (3, 3.5, 4 and 5 Å).
  The same geometric exclusion rule applies to both observed and shifted
  ions (none, or minimum solute heavy-atom separation of 1.8, 2.1 or 2.4 Å).
  Check eligible ion counts alongside occupancy because steric filtering can
  change both. A minimum periodic target–target separation greater than
  twice the contact cutoff rules out a bridge by geometry. The exclusion sweep is a geometric sensitivity analysis,
  with water and mobile monatomic ions excluded from the blocker set. It
  provides neither force-field insertion energies nor chemical cooperativity.

The release provides one trajectory per condition. Frame–offset samples are
correlated references, and condition-level inference requires separate
production replicas.


### Water coordination along attached HS

Set `analyses.sdc4_mechanisms.hs_sugar_water_cutoff_A: 3.5` to count distinct water
oxygens within 3.5 Å of any oxygen on each sugar. A water can contribute once
to each of several neighboring sugars. MDscope starts at the sugar bonded to
the protein and numbers along the covalent chain toward the distal end.
Linker sugars are included. Sugar types are the exact topology residue names.

The replica stage writes `extensions/sdc4_mechanisms/hs_sugar_hydration.parquet`
with frame, chain, position, sugar type and water count. It adds this output in
a dedicated trajectory pass and reuses it on subsequent runs with the same
cutoff and frame grid. To change the cutoff, move that output aside and rerun
the replica stage. Statistics use the existing profile bootstrap. Profiles
`hs_water_by_position__*` compare like positions and sugar types within the
same attachment chain. Profiles `hs_water_by_sugar_type__*` average positions
of a given type per frame within each chain, then compare the replicas. Type
averages retain the positions at which that type occurs in the topology.
