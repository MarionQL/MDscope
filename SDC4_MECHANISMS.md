# SDC4 mechanism analyses: eight opt-in endpoints

The \`analyses.sdc4_mechanisms\` observer uses the **existing single raw-frame
trajectory pass**. It never aligns or reopens the XTC. The
\`analyses.shared_states\` project finalizer reads saved per-replicate tables
and fits **one** state model across all conditions and replicates.

**Do not enable these analyses before validating topology selections.**
\`chainID\`, \`resid\`, glycolipid anchor names and HS attachment mappings are
specific to the actual TPR and simulation build. The example below uses
the older SDC4 PROA/PROB numbering as an *illustration*, not a guarantee for
rebuilt mature-protein systems.

## Config example

Add this under \`analyses:\` in \`synd4_baseline.yaml\`, adjusting all selectors:

\`\`\`yaml
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
\`\`\`

**Important:** \`allow_absent_species\` skips missing lipid species (e.g. GM1a
in a no-headgroup condition). To compare density maps across conditions,
use the same species list and report absent species as biologically absent,
not as zero-density enrichment. Density map coordinates are rotated so the
x-axis points from TM chain A to chain B. Leaflet-resolved density is in
lipids/nm² and fold enrichment uses each leaflet's mean species density.
Maps are **not** corrected for the protein-excluded accessible area.

## Observable definitions and scientific limitations

| Endpoint | Output/interpretation |
|---|---|
| Lipid enrichment and density maps | Per-leaflet GM1a/BGAL/cholesterol/phospholipid maps centered on TM midpoint, plus the existing annular \`lipid_enrichment\` fold values when enabled |
| Ectodomain height/tilt/contact | CA-center height relative to extracellular leaflet, terminal-to-terminal tilt to the outward membrane normal, binary membrane-contact state and fraction of ectodomain atoms contacting selected lipid atoms |
| NxIP accessibility | Regex matches actual CA-derived topology sequence by chain; reports fraction of motif heavy atoms with protein/HS/lipid steric neighbors within cutoff. This is an **occlusion proxy**, not a ligand-binding free energy. When `sasa.enabled`, the observer also joins the existing FreeSASA motif residue SASA, calculated with that analysis's atom selection |
| Individual HS orientation/contacts | Attachment-to-distal-terminus orientation cosine, extension, direct HS–lipid contacts and separate attachment-site contacts; distal end chosen by attachment proximity on the first analyzed frame |
| TM crossing and GXXXG | PCA helix-axis crossing angle, handedness and cross-chain GXXXG heavy-atom/residue contacts; no GXXXG motif produces zero contacts unless \`require_gxxxg\` is true |
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
inspect the \`selection_audit.json\` and supply a more appropriate selection
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

## Outputs

Per replicate:

- \`extensions/sdc4_mechanisms/timeseries.parquet\`: all frame-level scalar
  geometry, contact, hydration and membrane endpoints
- \`extensions/sdc4_mechanisms/lipid_density_maps.npz\`: per-leaflet 2D
  counts, density and fold maps, grid edges and frame count
- \`extensions/sdc4_mechanisms/selection_audit.json\`: exact matched
  motifs, GXXXG residue indices, HS termini and lipid counts
- \`extensions/sdc4_mechanisms/hs_conditioned_solvation.parquet\`:
  HS hydration and ion coordination conditioned on direct lipid contact
- \`extensions/sdc4_mechanisms/lipid_{density,enrichment}_maps.png\`
- \`extensions/sdc4_mechanisms/timings.json\`

Project-level (when \`shared_states.enabled\`):

- \`project_summary/shared_states/model.json\`
- \`state_centroids.parquet\`, \`frame_assignments.parquet\`
- \`replicate_populations.parquet\`, \`system_populations.parquet\`
- \`sampling_diagnostics.parquet\`, \`state_conditioned_hs_solvation.parquet\`
- \`state_populations.png\`, \`state_sampling_diagnostics.png\`

The scalar \`timeseries.parquet\` columns are also merged into MDscope's
existing \`summary/timeseries.parquet\` and analyzed by the existing
replicate-aware statistics framework. No new trajectory pass is added.


## Automatic covalent HS chain discovery

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
`REPORT.md` for PDF evidence-report generation and upload guidance.
