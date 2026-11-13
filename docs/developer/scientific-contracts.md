# Scientific and numerical contracts

These are regression invariants, **not** proof of sampling convergence or biological mechanism. Inspect [Scientific audit](scientific-audit.md) for unresolved limitations.

## Coordinate provenance and periodic geometry

| Role | Source | Allowed use |
|---|---|---|
| `source` | Primary XTC, **before** fitting | Membrane, contacts, solvation, local lipid/HS geometry |
| `msd_source` | Optional independent `msd_trajectory` | MSD only; topology/order, frames and time stamps must match primary |
| `fitted` | Derived protein-aligned NumPy arrays | Structural RMSD and other fit-defined metrics; **never diffusion** |

- `source` means **as loaded**, not necessarily unprocessed. Pre-rotated/recentered XTC can bias membrane observables even when MDscope itself never rotates the frames.
- Minimum-image incremental lipid steps assume intersample motion is unaliased; this is **not** guaranteed whole-molecule unwrapping.
- `pipeline/streaming.py` evaluates physical observers before protein fitting; `core/trajectory_sources.py` synchronizes optional MSD input.
- Coordinate positions/cutoffs use **Å**; MDAnalysis box dimensions are three lengths and three angles.

### Similar helpers with different scientific meanings

| Contract | Preserve |
|---|---|
| `core.geometry.minimum_image_displacements` | Shared triclinic displacement kernel, validated box and shape; caller dtype policy matters |
| `biophysical_common.lateral_distance` | **All-pairs** N×M lateral matrix, orthorhombic only |
| `mechanistic_common.planar_distances` / `core.lipid_geometry.lateral_distances` | Point-to-anchor lateral distances; different projection/shape and PBC policies |
| `AnchorTracker.positions` | Molecule-owner, PBC-aware center; not the same as simple arithmetic COG |
| `observables.msd._residue_centers` | Mass-weighted COM, with missing-mass COG fallback |
| `core.lipid_geometry.classify_leaflets` / `core.leaflet_partition.median_leaflet_mask` | Periodic circular grouping vs fixed median cohorts |
| Heavy-atom policies | Element-aware/name-based hydrogen inference differ for incomplete topologies |
| Motif matching | Chain identity and overlapping/nonoverlapping regex matches are separate contracts |
| Sparse contact wrappers | Selection-local atom-pair indices, optional exclusions, count vs pair output, precision |

**Do not deduplicate** pointwise vs pairwise, mass-weighted vs equal-weight, dynamic vs fixed leaflets, or different missing/invalid policies solely because functions look similar. Validate geometry against independent MDAnalysis/direct-sum references and actual TPRs.

## Measurements and estimands

| Measurement | Definition | Essential limitation |
|---|---|---|
| Single-origin lateral MSD | Mean squared displacement from integrated minimum-image steps (Å²) | Drift, unwrap aliasing and fit-window sensitivity |
| Multi-origin MSD | FFT average over available origins; separate on-disk second pass | Origins/lag bins are correlated; fewer origins at long lag |
| Diffusion | `D (nm²/ns) = 2.5 × lateral MSD slope (Å²/ps)` from `MSD = 4Dt` | Descriptive slope; high R² does **not** prove diffusion |
| Contacts | MDAnalysis capped search; indices local to selected AtomGroups, collapsed by **topology resindex**, not repeating `resid` | PBC requires current box; chemistry and selection are external inputs |
| H bonds | Donor–acceptor cutoff and D–H···A angle; lipid bond graph excludes same molecule | Charge/protonation and acceptor guess need validation |
| SASA | Probe-inflated geometric surface, nm²; environment burial = isolated minus occluded | Element radii, whole target, compatible frames and compact periodic target |
| Lipid order | `-S_CH` from bonded C–H vectors, by tail/site/leaflet | Initial cohort vs dynamic leaflet semantics |
| Local bilayer geometry | Projected lipid-anchor heights/tilt/quadratic small-slope curvature | Planar approximation; mixed headgroup heights are not interchangeable |
| Solvation counts | Distinct solvent entities near selected atoms; per target | Different polar atom numbers/site sizes confound raw counts |
| Solvation RDF | Pair density relative to **full-cell** selected solvent density | Excluded protein/HS/bilayer volume biases `g(r)`; no intrinsic binding interpretation |
| Shared states | Replica-balanced standardized k-means of geometric features | Geometric clusters, **not** validated metastable or causal states |

### Membrane geometry and populations

- Circular leaflet observers require the configured normal perpendicular to lateral cell vectors; **in-plane shear** may be supported, normal-direction shear fails.
- SDC4 height maps, projected footprint enrichment and certain local maps require **orthorhombic** boxes (and their configured axis limits). Curved/overhanging bilayers need another surface model.
- Initial median APL/MSD cohorts, initial order leaflets and dynamic enrichment/normals are **different populations**. Do not reinterpret one as another.
- One ceramide anchor per bonded GM1a/BGAL molecule; the glycan residues can share names with HS. Validate TPR chemical selections and missing species.
- PBC lipid footprint uses a leaflet-specific projected protein disk union; accessible area is grid quadrature, **not** molecular free area or affinity.

## Local lipid composition and residence

For species `s` in local region `R` and leaflet `L`:

`A_s` = available molecular anchors inside the prescribed disk; `C_s` = those **same molecules** making contact.

| Metric | Definition | Missingness |
|---|---|---|
| Matched contact fraction | `C_s / A_s` | Species absent → NaN |
| Matched composition fold | `(C_s / ΣC) / (A_s / ΣA)` | No contacts in any species → NaN |
| Original whole-leaflet fold | Contact fraction relative to whole leaflet species population | Different denominator; keep distinct |
| Legacy unmatched local fold | Historical comparison only | Do **not** substitute for corrected matched fold |

One molecule counts once even if it has multiple bonded glycan residues. A frame with no contacts has contact fraction **zero** for present species but undefined composition enrichment. Corrected summary columns: `mean_matched_local_contact_fraction`, `mean_matched_local_enrichment`. Old outputs without them remain missing, not fabricated from legacy columns. Contact folds are **not** affinity.

For a sampled contact starting at `i`, ending at `j`:
- **Observed span:** `t[j] - t[i]`.
- **Interval-imputed duration:** span + next observation interval; at trajectory end use the median interval. One-observation trajectories have unavailable imputation.
- Lipid `duration_ps` uses the **imputed** convention; solvation `duration_ps` uses **observed span**. `sampled_duration_ps`, `observed_span_ps` and `duration_convention` disambiguate.
- `occupancy_fraction`: molecule's contacting frames / total frames; `episode_contact_fraction`: episode contacting frames / total frames. Neither is time-weighted on irregular sampling.
- Intermittent gaps admit at most configured **absent frames**, whose elapsed time is part of the episode span. Left/right boundary censoring remains explicit.
- Kaplan–Meier survival excludes left-censored origins, retains right-censored at-risk observations and handles ties. Temporal undersampling means true physical lifetimes remain unresolved.

## Solvation/bridge evidence

| Evidence | Means | Does **not** mean |
|---|---|---|
| `rdf.parquet:g_bulk` | Historical full-cell normalized RDF | Solvent-excluded-volume corrected bulk RDF |
| Optional `g_aqueous_reference` | Density from a specified membrane-excluded orthorhombic slab | Shell-accessible-volume corrected RDF or measured true bulk |
| Geometric bridge | **Same solvent/ion identity** near two distinct targets in the same frame | H-bond-qualified or cooperative bridge |
| `water_hbond_bridges` | Two different H on a water each donate to a target acceptor | All possible donor–acceptor bridge orientations |
| Pair shift-reference | Up to 100 deterministic nonzero shared temporal offsets; relative occupancy | Confidence interval, p-value, spatial control or cooperativity |
| Spatial ion reference | In-plane translations retaining ion heights/relative arrangements | Sterically valid ion placement or free-energy null |

- Configured pair traces include **all** sampled frames and measured zeros. Unknown/untracked pair → error, not zero. Irregular or fewer than four frames → shift reference unavailable.
- Aqueous-region density uses `V_reference = V_cell (1 − 2h/L)` for axis length `L` and slab half-width `h`; unsupported rotations/triclinic or empty region are unavailable/fail explicitly. A hemisphere remains half-occupied even at correct reference density.
- A bridge shift breaks **time alignment**, not crowding. The steric K audit adds geometry/exclusion sensitivity but no insertion free energies.
- Radial bins, cutoff candidate minima, and sampled survival curves are not independent evidence of chemical binding.

## Replicate statistics and identity

- Condition-level independent units = **simulation replicas**. Paired inference matches replica labels; irregular/missing trajectory frames are not independent resamples.
- Scalar/profile/interactions: block within replica, hierarchy across replicas, exact/Monte Carlo permutations and BH FDR when configured. RDF bins are descriptive mean/SD; shell-integrated RDF and shared-lag diffusion use **replica-level** contrasts.
- Nonfinite/gapped/irregular series retain supported descriptive points but **not** unsupported temporal CIs or ACFs. Missing required replicas/coordinates cannot be silently omitted.
- A 3+3 unpaired two-sided exact permutation has 20 assignments, minimum p = 0.10. A bootstrap percentile CI excluding zero is **not** FDR-significant evidence.
- Compare common units, labels, molecular identities, chemical references, frame windows and normalization before interpreting a difference.

## Restart, completeness and verification

`registry.required_artifacts` describes companion files for configured analyses, not all possible outputs. Required products are checked before resume/finalization: primary columns and Parquet metadata, JSON objects, NPY shapes, NPZ CRC and relevant time/identity agreement. Disabled/optional products are not required.

- A *readable* artifact may still be scientifically stale. Compatible legacy results are reusable without hash-based gating; changed scientific definitions require deliberate regeneration.
- Missing/corrupt required companions must produce safe failure or reconstruction **without destroying valid outputs**.
- Checkpoint fingerprints use metadata, not a cryptographic content guarantee. Parquet footer integrity is not proof every data page or chemical selection is correct.
- Streaming contact-episode extraction is bounded; **other** finalization stages may still materialize tables. Do not infer production memory from toy tests.

## Reference tests

```bash
pytest -q tests/test_core_geometry_contracts.py tests/test_multi_origin_msd.py
pytest -q tests/test_publication_references.py tests/test_scientific_reference_oracles.py
pytest -q tests/test_mixed_observer_restart.py tests/test_planar_membrane_contracts.py
```

Independent numerical references cover PBC, mass/COG, sparse identity, exact permutations/FDR/bootstrap, area, SASA cap/sphere, MSD direct sums, residence and missingness. [Portable release validation](release-validation.md) adds saved-output oracles; [SDC4 fixtures](sdc4-fixtures.md) cover SDC4 geometry. Real GROMACS checks and six-production-replica convergence remain separate.
