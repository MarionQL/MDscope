# Scientific analyses

**Start with [the full SDC4 YAML](../synd4_full_mechanisms.yaml)** for real selections and [Configuration](configuration.md) for syntax. This page defines **what measurements mean**, where their data live, and what they **cannot** establish.

## Choose an observable

| Question | Analysis | Main interpretation |
|---|---|---|
| Is protein conformation changing? | RMSD, RMSF, Rg, DSSP, SASA | Fit-defined structure, flexibility, secondary structure and accessible surface |
| Does the bilayer change? | Membrane APL/thickness, lipid order, leaflet normals, density | Packing/geometry proxies, **not** automatic fluidity or mechanical moduli |
| Are lipids recruited near SDC4? | Lipid enrichment, contact matrices, lipid H bonds | Composition, proximity, molecular contact and H-bond geometry |
| Does lipid mobility change? | Single-/multi-origin MSD | Descriptive lateral diffusion **slope**, fit-window dependent |
| Do waters/ions reorganize? | Solvation, RDF, hydration profiles | Geometric coordination and sampled persistence |
| Do HS/protein regions move differently? | SDC4 mechanisms, matched profiles, shared states | Distances, tilt, occlusion, contacts, conformational grouping |
| Are two conditions different? | Statistics, coupling and coverage | Replica-aware contrasts with explicit missing-data status |

**Units:** coordinates/contact cutoffs in Å; SASA in nm²; MSD in Å²; diffusion in nm²/ns. Frame time is in ps unless identified otherwise. Three independent simulations per condition constrain inference even with thousands of frames.

## Membrane and lipid observables

| Measurement | Definition | Important assumption |
|---|---|---|
| APL, thickness | Projected periodic area/lipid; separation of selected leaflet anchors | Approximately planar membrane, comparable anchor chemistry |
| Lipid `-S_CH` | Bonded C–H orientational order by sn-1/sn-2 site, species and initial leaflet | Requires bonded all-atom hydrogens; fixed initial leaflet assignment |
| Local bilayer normal | Plane fit to neighboring **same-leaflet** anchors; deviation from global axis (degrees) | Measures local tilt/undulation, **not** elastic modulus |
| Leaflet composition enrichment | Local species mole fraction / same-leaflet reference fraction | Defined populations; unavailable species must stay absent |
| Protein footprint enrichment | Lipid anchors in bands outside projected protein-disk union; species counts per accessible area | Geometric disk/slab/grid model, **not** binding affinity |
| Direct lipid proximity | Fraction of intact lipid molecules within 3D contact cutoff | Full bonded molecules required for multi-residue glycolipids |
| Lipid–lipid hydrogen bonds | MDAnalysis donor/acceptor geometry; bonded same-molecule contacts excluded | Donor/protonation chemistry and complete connectivity |
| Lipid contact matrices | Unique bonded molecule pairs in each leaflet | Random-mixing normalization, not size-corrected affinity |

### Geometry and chemistry

- Choose **one anchor per physical lipid**. Ceramide-based GM1a/BGAL may span multiple topology residues; do not count every sugar as a lipid.
- Normals, enrichment and lipid contacts use circular leaflet geometry and allow **in-plane** periodic shear; shear through the configured normal fails.
- SDC4 local maps/footprint enrichment require **orthorhombic** boxes. Initial APL/MSD cohorts, fixed lipid-order leaflets and dynamically reassigned leaflets are different populations.
- Empty/undefined denominators → NaN. A present species with measured zero contacts → zero contact fraction; composition enrichment can still be undefined.
- Global membrane-relative density centers each frame on the midpoint of circular phosphate-leaflet means. It follows **translation**, not tilt or undulation.

The two **audited release TPRs** use glycolipid **chain IDs**, not `resname GM1a`. Example anchors:

```yaml
species_selections:
  POPC: "resname POPC and name P"
  CHL1: "resname CHL1 and name O3"
  GM1a: "chainID GM1a and resname CER160 and name NF"
  BGAL: "chainID BGAL and resname CER160 and name NF"
  CER160: "chainID CER160 and resname CER160 and name NF"
allow_absent_species: true
```

| Release topology | GM1a | BGAL | Free CER160 |
|---|---:|---:|---:|
| Glycolipid-containing | 24 | 58 | 0 |
| No glycolipid headgroups | 0 | 0 | 82 |

These counts and identifiers are **release-specific**. Revalidate rebuilt TPRs; `resname BGAL` alone also matches some HS sugars. For bonded lipid H bonds, select complete molecules (phospholipids/cholesterol **or** glycosphingolipid chain IDs). Incomplete bond topology invalidates molecule-level exclusions.

### Lipid neighbors and enrichment

For species A/B on an instantaneous leaflet, the unique undirected pair count `C_AB` excludes self-pairs; AB and BA table rows represent the same contacts.

- Mean neighbors per molecule: `C_AB / N_A` (A≠B), `2C_AA / N_A` (A=B).
- Possible pairs: `N_A N_B` (A≠B) or `N_A(N_A−1)/2` (A=B).
- Random-mixing enrichment = observed contact share / possible-pair share. A complete graph gives **1**; this is **not** an affinity or steric-size correction.
- Absent population/reference → NaN; an available pair with zero contacts gives zero when the leaflet has other contacts.
- Output: `extensions/lipid_contacts/{matrix,timeseries}.parquet` and `selection_audit.json`. Heatmaps average **replicas equally**, mask missing species, and are descriptive.

For matched local SDC4 lipid specificity, `A_s` is the number of species-s anchors within a chosen local disk and `C_s` is the **subset of those same molecules** contacting the protein:

| Metric | Formula |
|---|---|
| Matched local contact fraction | `C_s / A_s` |
| Matched local composition fold | `(C_s / ΣC) / (A_s / ΣA)` |

No contacts anywhere → undefined fold; present species → zero contact fraction. Legacy **unmatched** local enrichment is kept separately; it cannot replace the corrected matched measure. See [Scientific contracts](developer/scientific-contracts.md).

Other lipid outputs: `extensions/bilayer_normal/`, `extensions/lipid_footprint/`, `interactions/lipid_hydrogen_bonds/`. Framewise ratios enter ordinary replicate-aware statistics. Check grid resolution, slab width, disk radius and anchor chemistry before interpreting small differences.

## Syndecan-4 (SDC4) mechanisms

`analyses.sdc4_mechanisms` joins the **shared source-coordinate pass**. `analyses.shared_states` fits **one** model across saved replicas during offline finalization. Validate topology selections before running.

| Endpoint | Reported quantities | Interpret with care |
|---|---|---|
| Ectodomain | Height, outward tilt, membrane-contact state and contact fraction | Terminal orientation assumes ordered CA selection; not full flexible-surface fitting |
| TM dimer | Crossing angle, handedness, GGIVG/GXXXG contacts and full interface contacts | Missing motif ≠ zero contact; PCA axes can become unavailable |
| NxIP/motifs | Chain-specific sequence matches; fraction of heavy atoms occluded by protein/HS/lipid; optional motif SASA | **Steric occlusion proxy**, not ligand-binding energy; protein-only SASA omits HS/lipid shielding |
| Individual HS chains | Attachment orientation, extension, Rg, anisotropy, lipid contact and membrane proximity | Bond topology; compact-chain minimum image; branches need appropriate endpoints |
| Hydration/ions | Unique oxygen/ion neighbors, per-residue states and conditioned means | Coordination, **not** electrostatic free energy |
| Local membrane | Thickness, quadratic curvature, local `-S_CH` | Comparable phosphate anchors, planar/small-slope approximation |
| Lipid density maps | TM-centered, per-leaflet maps, molecule density and fold | Maps do **not** divide by protein-excluded accessible area |
| Shared states | Replica-balanced standardized k-means; population, transitions and first/second-half diagnostics | Geometric clusters, not validated metastable states or thermodynamic free energies |

**Release topology:** each protomer has three **covalently attached** HS chains. Attachments: PROA **S39/S61/S63**, PROB **S309/S331/S333**. HS sugars share protein chain IDs, so `hs_auto` follows the **bonded glycan components**, requiring one protein attachment per chain. Rebuilt topologies need re-auditing. The older release found GGIVG in **PROA but not PROB**; absent optional motifs are unavailable, not measured zero.

- Local thickness/curvature use `surface_anchor_species` of comparable lipid headgroup heights. Quadratic curvature is a small-slope proxy in Å⁻¹.
- Protein occlusion excludes same-chain neighbors within ±2 sequence positions; partner fractions can overlap, and union fractions count each target atom once.
- State training balances input frames **per replica**. Missing-feature frames and temporal gaps change support and transitions; block-bootstrap of frames is not condition-level inference.
- Outputs: `extensions/sdc4_mechanisms/timeseries.parquet`, `selection_audit.json`, `lipid_density_maps.npz`, `hs_conditioned_solvation.parquet` and optional residue/motif/hydration tables. Project states: `project_summary/shared_states/`.
- A complete reference configuration is [`synd4_full_mechanisms.yaml`](../synd4_full_mechanisms.yaml). **Do not copy old-residue examples blindly.**

### Protein residue water and HS sugar hydration

| Option | Definition | Saved output |
|---|---|---|
| `residue_occlusion: true` + `residue_water_cutoff_A: 3.5` | Distinct water oxygens within cutoff of any selected protein N/O; protein/HS/lipid/no-partner conditional means | `extensions/sdc4_mechanisms/residue_occlusion.parquet` |
| `hs_sugar_water_cutoff_A: 3.5` | Distinct waters near **any oxygen per sugar**, numbered from attachment sugar outward along covalent HS chain | `extensions/sdc4_mechanisms/hs_sugar_hydration.parquet` |

- A water may count once for **multiple neighboring sugars**; per-sugar count is not a conserved chain total.
- HS profiles: `hs_water_by_position__*` (chain/position) and `hs_water_by_sugar_type__*` (type mean per frame, with exact positions documented).
- Residues without polar atoms have unavailable coordination. Conditional CIs preserve missing frames and require support across multiple temporal blocks.
- HS sugar hydration is a separate **optional trajectory pass**; saved cutoff and frame-grid compatibility are checked before reuse. Change the cutoff using fresh output or controlled regeneration.

## Single- and multi-origin diffusion

### Single-origin membrane MSD

`analyses.membrane.msd.enabled: true` uses the shared frame pass. The lateral 2D Einstein relation is `MSD = 4Dt`, so:

`D (nm²/ns) = 2.5 × slope (Å²/ps)`.

Per-replicate fits and R² are descriptive, not formal independent-replica uncertainty; output: `summary/diffusion_coefficients.parquet` and `project_summary/diffusion/`.

### Optional multi-origin MSD

```yaml
analyses:
  membrane:
    msd:
      enabled: true
      multi_origin: true
      multi_origin_max_lag_frames: 1000  # null = all sampled lags
      multi_origin_particle_batch: 64
      multi_origin_fit_diagnostics: true
```

- **Separate pass:** reconstruct lipid-center displacements by incremental minimum image (optionally remove common drift); stage float64 coordinates in an on-disk memmap.
- Approximate temporary disk bytes: `24 × sampled_frames × selected_lipids`; lipid batches bound FFT memory. Temporary coordinates are deleted after success/error.
- Average all overlapping origins; lag is in **sampled frames** and lag zero is excluded from fitting. Require regular sampling, ≥4 positive lag points and a documented fit window.
- Outputs: `summary/multi_origin_msd.parquet`, `summary/multi_origin_diffusion.parquet`; project `diffusion/multi_origin_{replicates,systems}.parquet`.
- Fit diagnostics read **saved** curves. Earlier/later overlapping windows report D, R², log–log slope and origin support: ≥10 origins or 10% of sampled frames (whichever larger). Heuristic flags: D spread >25% or exponent >0.20 from 1.
- Neither a high R² nor a screen pass establishes asymptotic diffusion. Interframe PBC aliasing, correlated origins, finite size, drift, wrapping and long-lag support remain checks.

## Solvation: water and ions

`analyses.solvation` shares the fused **raw-frame** pass. One sparse periodic search per solvent class/frame provides target counts and RDF histograms. Targets and ion species must be nonoverlapping selections; choose `group_by: residue`, `fragment` or bond-derived `hs_chain` appropriately.

```yaml
analyses:
  solvation:
    enabled: true
    preset: standard
    water_selection: "resname TIP3 and name OH2"
    ion_selections: {K: "resname POT", Cl: "resname CLA"}
    targets:
      protein:
        acidic:
          selection: "protein and ((resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2))"
          group_by: residue
    search_radius_A: 8.0
    rdf_bin_width_A: 0.2
    shell_cutoffs_A: {water: 3.5, ions: 4.0}
```

**This is an illustrative selector**; audit the actual CHARMM/TPR names, protonation, water hydrogens, HS sulfate oxygens and K/Cl chemistry.

| Preset / option | Adds | Cost / meaning |
|---|---|---|
| `standard` | Unique target-solvent counts, RDF, species timeseries | No solvent-identity event tables |
| `detailed` | Identity events and candidate geometric bridges | High cardinality for all water |
| `kinetics` | Identity tracking and censored sampled-frame episodes/survival | **Not** continuous physical residence |
| `options.track_target_classes` | Restrict expensive identities to selected classes | Counts and RDF still include all target classes |
| `options.water_hbond_bridges: true` | Two different water H each donate to a distinct acceptor (angle ≥150°, O···A ≤3.5 Å defaults) | Conservative subset; excludes target-donor→water bridges |
| `grid_frames` | Three ordered noncollinear atoms per selected site for moving-frame 3D water density | Maps depend on local coordinate frame and sampled sphere |

### RDF normalization, exposure and sampling

- **Historical `rdf.parquet:g_bulk`:** solvent count / full instantaneous periodic-cell volume. It is **not** corrected for protein, lipid or glycan excluded volume.
- Optional `analyses.solvation.aqueous_reference` uses an orthorhombic membrane-centered **excluded slab** and writes `g_aqueous_reference` separately. For box length `L`, half-width `h`: accessible reference *geometric* volume `V × (1 − 2h/L)`. It does **not** correct the inaccessible part of individual radial shells.
- Half-shell accessibility can yield `g(r)=0.5` even at the correct aqueous density. Density normalization alone does not prove hydration preference.
- Cutoffs are **specified**, not automatically inferred. `shell_candidates.parquet` suggests first minima; `shell_stability.parquet` checks sampled chunks (four resolvable chunks and ≤0.5 Å spread for stable candidates). Use the **same justified chemical-site cutoff** across conditions.
- `hydration_sasa.parquet` joins saved per-residue SASA at exact times; missing constituent residues → NaN, not partial SASA. This is contextual exposure, not a correction to RDF.
- Contact episodes and Kaplan–Meier survival represent **sampled observations**; left-censored episodes excluded, right-censored at risk. Frame/episode counts are **not independent replicas**.

### Saved solvation products

Under `extensions/solvation/`:

| File(s) | Contents |
|---|---|
| `counts.parquet`, `species_timeseries.parquet`, `target_summary.parquet` | Unique per-target counts; per-class mean and hydrated fraction (≥1 solvent); descriptive per-target summary |
| `rdf.parquet`, `shell_candidates.parquet`, `shell_stability.parquet` | Normalized RDF, candidate cutoff and chunk stability |
| `target_residues.parquet`, `target_atoms.parquet`, `metadata.json` | Topology/atom/selection identity and audit |
| Optional `events`, `bridges`, `episodes`, `survival` | Identity and sampled-contact histories |
| Optional `hydration_sasa`, `aqueous_reference`, `water_density_3d` | Exposure join, reference density, local 3D map |

Species class averages weight **each configured target equally per frame**; a single water can contribute to more than one target. Ion species `ions/K` and `ions/Cl` stay separate. Legacy `ion_selection` aggregates ions; use **either** that or `ion_selections`.

Existing temporal/replica bootstrap processes `species_timeseries` (`mean_n_solvent`, `hydrated_fraction`). RDF project mean/SD per bin is descriptive; shell-averaged replica contrasts require full support. No new trajectory pass is needed for statistical reuse of these outputs.

## Bridge evidence and potassium-specific audit

A geometric bridge requires **the same solvent entity** near both configured targets in the *same frame*. It does **not** imply a chemically qualified bridge or cooperativity.

```yaml
options:
  bridges: true
  bridge_target_pairs:
    - ["protein/carboxylate/7", "hs/sulfate/200"]
  bridge_spatial_ion: ions/K
```

Use actual `target_id` strings from `target_residues.parquet`; example IDs are placeholders. Both targets must exist and be tracked.

| Output / control | Interpretation |
|---|---|
| `bridge_pair_summary.parquet` | Every configured eligible pair/species, **including observed zeros** over all sampled frames |
| `bridge_pair_timeseries.parquet` | Same-entity count and presence per frame; valid input to replicate statistics |
| `bridge_pair_shift_reference.parquet` | Up to 100 nonzero shared circular time shifts, preserving each entity and occupancy |
| `bridge_pair_control_summary.parquet` | Observed−shift mean, reference percentiles and availability |
| `bridge_spatial_reference.parquet` | Eight fixed in-plane ion translations/frame preserving ion heights and relative arrangement |
| `bridge_residence*.parquet` | Continuous/intermittent sampled same-entity episodes and censoring |

- **Unknown** pair → error, not zero. Time-shift null unavailable for irregular or fewer than four observations.
- Shift percentiles are **not** confidence intervals, p-values or independent trials. The time reference breaks alignment, not crowding; spatial translations may place ions inside solute volume.
- Episode finalization streams frame-ordered contact identities and keeps active episode state; some survival/bridge/SASA joins still have non-bounded-memory stages.
- Molecular null tests and distinct K contact-vs-solvent-separated shell boundaries need independent chemistry and spatial evidence.

### Potassium bridge release audit

`analyses.potassium_bridge_audit` performs a **separate** per-replica pass, saving results under `extensions/potassium_bridge_audit/`. Run a release topology/trajectory check directly with:

```bash
python scripts/audit_potassium_bridges.py \
  --tpr real_test_data/carb-sdc4-step7_production.tpr \
  --xtc real_test_data/carb-sdc4-rep1_test.xtc \
  --stride 10 --output potassium-audit-carb-rep1
```

- `chemical_sites.parquet`: bond/charge-audited acidic patch and attached HS carboxylate/sulfate/sulfamate oxygens; unbonded +1 POT identity.
- `rdf.parquet` and `rdf_peak_candidates.parquet`: 0.1 Å oxygen–K radial bins, full-cell reference, first/second peak and intervening minimum candidates.
- `bridge_frames.parquet`, `bridge_summary.parquet`: same-K occupancy by patch, HS chain/class; contact cutoffs **3, 3.5, 4, 5 Å** and steric blocker-distance thresholds **none/1.8/2.1/2.4 Å**.
- Shared-water evidence requires one water within **3.5 Å** of K and donating an H bond to the target oxygen (≤3.5 Å donor–acceptor; ≥150°).
- Compare eligible-ion populations and observed-versus-shifted reference under identical steric rules. This is geometric sensitivity, **not** insertion energy or chemical binding free energy.
- Release files have only one trajectory per condition; use independent production replicas for condition inference.

## Practical validation and performance

1. Validate **both topologies** first: actual atom names, six HS attachments, intact lipid fragments, one anchor/molecule, target overlaps, water/ion identities and frame timing.
2. Run small standard solvation windows; inspect per-species counts, RDF shape, reference density, SASA and cutoff stability. Then activate limited kinetics/bridge sites.
3. Record total observer throughput, RSS and Parquet bytes on representative real frames. Optional benchmark: `scripts/benchmark_solvation_real.py`; batch/contact benchmarks live under `scripts/`.
4. Check all individual replicas, frame windows, direct measured zeros versus missing fields, blocked/censored episodes and small-n test resolution.
5. Use [Scientific audit](developer/scientific-audit.md) for **measured** historical microbenchmarks; they do not predict six-replica runtime.

## Adding an analysis

- Choose output contract first: scalar series, labeled time×position profile, interaction/event identities, or one estimate **per replicate**.
- Implement under the owning `src/mdscope/{observables,interactions,mechanisms}/` package; share geometry and sparse MDAnalysis kernels.
- Register outputs/tasks and **required companions** in `src/mdscope/registry.py`. Feed `summary/timeseries.parquet`, shared profile collectors or interaction tables.
- Cheap per-frame observables join the fused pass; expensive separate readers must be restartable. Do not add a custom statistical bootstrap or frame-pooled hypothesis test.
- Test a hand-calculated numerical case, invalid/missingness and required artifact behavior, then real TPR/XTC before production interpretation.

[Workflows](workflows.md) · [Scientific contracts](developer/scientific-contracts.md) · [Validation gaps](developer/validation-gaps.md)
