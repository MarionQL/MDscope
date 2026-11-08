# Scientific and numerical contracts

These contracts protect numerical semantics during refactoring. Historical audits are consolidated here to retain important decisions without issue-numbered document filenames.

---

## Trajectory and coordinate contracts

MDscope uses three coordinate roles; their labels specify **provenance inside
MDscope**, not guarantees about how an input XTC was prepared.

| Role | Origin | Consumers | Requirements |
| --- | --- | --- | --- |
| `source` | `rep.trajectory`, before protein fitting | membrane observables, density, protein/HS–lipid contacts | Keep the membrane geometry and box physically meaningful; no protein-fit transforms |
| `msd_source` | Optional `rep.msd_trajectory`, independent Universe | membrane MSD only | Same topology/atom order, frame count, sampling and timestamps as primary; no rotation or per-frame recentering |
| `fitted` | Fresh NumPy arrays produced by rigid-body protein fitting | aligned structural RMSD and related metrics | Document fit selection and reference frame; never reuse for diffusion |

**Important:** `source` means *as loaded*, not necessarily raw, wrapped,
unwrapped, or centered. A previously preprocessed XTC may be unsuitable for
physical membrane observables. In particular, applying time-dependent protein
rotation or centering to the input trajectory can bias diffusion even when
MDscope itself does not modify that input.

MSD reconstructs displacements from minimum-image *incremental* steps. It
assumes each sampled per-lipid step is small enough that the minimum image is
unambiguous. This is not equivalent to guaranteed unwrapping at coarse
sampling intervals.

#### Execution boundary

In `analysis/pipeline/streaming.py`, fused membrane, extensions, densities and
membrane metrics are evaluated **before** fitted structural arrays are
computed. The separate MSD Universe is frame-synchronized in
`core/trajectory_sources.py`. Alignment constructs derived arrays; it does
not intentionally rotate the main Universe.

`core/coordinate_contracts.py` supplies explicit labels and a small
guard for programmatic handoffs. It does not attempt to infer PBC or
preprocessing history from coordinates. That would require additional input
metadata and validation and should be addressed separately.

#### Single-origin and multi-origin MSD

The standard single-origin estimator is evaluated incrementally in the fused
membrane streaming pass; its displacements are accumulated from minimum-image
steps. It reports framewise mean-squared displacement in Å² and a descriptive
per-replicate diffusion coefficient (nm²/ns).

The optional `analyses.membrane.msd.multi_origin: true` path **is now
implemented** in `observables/multi_origin_pipeline.py`. It makes a
separate, explicit trajectory traversal and stages reconstructed lipid-center
displacements in a bounded-RAM temporary disk-backed array, then calls the
pure NumPy FFT kernel in `observables/multi_origin_msd.py`. This deliberately
incurs additional I/O and uses the same minimum-image sampling assumption;
it does **not** reconstruct missing molecular unwrap history. See
[MULTI_ORIGIN_MSD.md](../analyses.md).

For publication-quality diffusion, validate the processing provenance of
input trajectories, fit-window sensitivity, sampling versus box displacement,
and independent reference trajectories. Neither an excellent linear-fit R²
nor a clean trajectory pass establishes the diffusive regime.

#### Contact and residue identity

`biophysical_common.atom_contacts` uses MDAnalysis capped-distance with
an optional periodic box; distances and cutoffs are in Å. It returns indices
**local to the two selected AtomGroups**, not global atom IDs. The pure NumPy
`collapse_atom_pairs_to_residues` function maps these indices to unique,
sorted **topology residue indices**; it must not collapse by the human-readable
`resid`, which can repeat between protein chains or lipid segments.
`residue_pairs` combines these two operations. This contract applies to
protein–membrane and membrane-only selections alike. A `None` box is
nonperiodic; callers must pass frame dimensions when PBC is required.

#### Remaining scientific checks

The standalone multi-origin adapter uses a second reader by design; do not
claim all optional observables are single-pass. Chemical selections and mass
fallback policies remain explicit; drift removal for both MSD estimators is
optional and must be recorded for each simulation.

---

## Geometry helper audit: shared operations versus distinct contracts

This audit compares the October 2026 implementations in
`biophysical_common.py`, `extended/common.py`,
`mechanistic_common.py`, and `lipid_geometry.py`.

### Consolidated now

The minimum-image transformation of *displacement vectors* is the same
MDAnalysis operation in both `extended.common.minimum_image` and
`mechanistic_common.planar_distances`. The implementations now call
`core.geometry.minimum_image_displacements`, preserving the extended
helper's explicit float32 coercion while keeping the mechanistic helper's
input precision. The shared function validates the box and coordinate shape.

### Deliberately separate

| Helpers | Why not merge without changing behavior |
| --- | --- |
| `biophysical_common.heavy`, `extended.common.heavy`, `mechanistic_common.heavy` | They differ in element-attribute availability, names starting with H, numbered hydrogen names, and behavior for missing topology elements. A single atom-selection policy needs separate chemistry/topology tests. |
| `biophysical_common.lateral_distance` | Calculates pairwise N×M distances; orthorhombic-only, configurable normal. |
| `extended.common.lateral_distance` | Calculates N point-to-one-center distances, fixed xy plane; float32 minimum image. |
| `mechanistic_common.planar_distances` | N point-to-one-anchor distances, configurable normal. |
| `lipid_geometry.lateral_distances` | Pointwise lateral offsets through its own PBC geometry and configurable normal. |
| `extended.common.centers` | Per-residue equally weighted coordinate averages, including documented multi-atom PBC restriction. |
| `lipid_geometry.AnchorTracker.positions` | Topology-fixed owner mapping, multi-atom PBC-aware averaging. |
| `biophysical_common.leaflet_groups`, `lipid_geometry.classify_leaflets`, `extended.common.load_lipid_anchors` | Different partition strategies, return types, and physical definitions; equivalence has not been established. |
| `box_of`, `valid_box` | Optional/nonperiodic fallback versus mandatory box with explicit error. |

### Follow-up safety work

1. Test the three hydrogen selection rules against real CHARMM/TPR topologies before standardizing.
2. Choose an explicit contract for triclinic lateral geometry and distinguish physical projection from deleting an array coordinate.
3. Audit the exact callers of each leaflet classification before unifying schemes.
4. Never change selection, return shape, dimensionality, precision, or failure policy solely to deduplicate code.

---

## numerical-contract — Scientific geometry consolidation audit

This audit distinguishes infrastructure duplication from intentional scientific
methods. **Never merge two operations solely because their names resemble each
other.** Cartesian positions and distances are Å; MDAnalysis dimensions are
three lengths and three angles.

| Implementations | Decision | Reason |
| --- | --- | --- |
| `core/lipid_geometry.planar_offsets` and `core/geometry.minimum_image_displacements` | Reuse shared kernel | Same MDAnalysis minimum-image displacement and dtype (float64), then omit the selected Cartesian axis. New tests compare to independent MDAnalysis reference, including triclinic PBC. |
| `extended/common.minimum_image` | Keep small adapter | Explicitly requests float32; changing dtype may affect historical results/performance. |
| `biophysical_common.lateral_distance` vs `mechanistic_common.planar_distances` vs `extended/common.lateral_distance` | Retain distinct contracts | Pairwise orthorhombic matrix versus per-anchor triclinic-compatible displacement versus XY-only float32; shape and accepted boxes differ. |
| `core/lipid_geometry.AnchorTracker.positions`, `extended/common.centers`, `observables/msd._residue_centers` | Retain distinct estimators | PBC-unwrapped equal-weight center versus wrapped COG per selected residue versus COM with COG fallback. The latter depends on masses. |
| `core/lipid_geometry.classify_leaflets`, `core/leaflet_partition.median_leaflet_mask`, `biophysical_common.leaflet_groups` | Retain | Circular two-cluster, Cartesian median, and circular-gap classification have different assumptions and output semantics. |
| `core/atom_selection.select_heavy_atoms` callers | Already consolidated, retain policies | Element-aware and name-based hydrogen detection differ for incomplete topologies. |
| `biophysical_common.atom_contacts`, `extended/common.pairs`, `mechanistic_common.contact_count` | Retain wrappers | Sparse output versus count, optional self-exclusion, and different float32/float64 semantics. All use capped-distance rather than dense NxN allocations. |
| `biophysical_common.sequence_motifs`, `mechanistic_common.sequence_matches`, `extended/common.sequence_hits` | Retain | Distinct nonoverlapping/overlapping match behavior, chain attribution and residue/object return contracts. |
| `mechanistic_common.orientation_vector` and `extended/common.axis_from_ca` | Retain | Minimum length 3 versus 4, CA selection policy, and invalid-value conventions differ. |

### Validation coverage

- `tests/test_numerical-contract_geometry_contracts.py`: independent MDAnalysis
  reference for minimum-image PBC displacement in rectangular and triclinic
  cells, circular centers crossing cell boundaries, invalid/empty geometry.
- `tests/test_numerical-contract_topology_contracts.py`: reordered atom-to-residue
  ownership, absent mass/element attributes, empty selections, selection policy.
- `tests/test_numerical-contract_glycolipid_motifs.py`: one ceramide anchor per physical
  glycolipid despite carbohydrate residues; wrapped versus unwrapped centers;
  overlapping species selection rejection; chain-specific motif identities
  and deliberate overlapping-match behavior.
- `tests/test_numerical-contract_contact_mass_contracts.py`: sparse atom/residue
  contact agreement with independent MDAnalysis calls; self-exclusion;
  COG fallback with missing mass and mass-weighted COM when masses exist.

All comparisons preserve currently documented output semantics. The
multi-residue fixture models distinct ceramide and sugar residues explicitly;
it does not claim that chemically connected residues share an MDAnalysis
`resindex`. Use a physical ceramide anchor for multi-residue species.

### Further consolidation policy

The remaining similarly named operations in the table above are **not**
demonstrated duplicates: their scientifically observable outputs differ.
No additional production-code merge is justified without selecting an
explicit common scientific contract, comparing real data, and obtaining
approval for any intentional numerical change. This issue's audit and
targeted mathematically equivalent kernel extraction are complete once
the CI/smoke suite verifies the added regression coverage.

---

## numerical-contract — numerical primitive audit

This table records **semantic contracts**, not just similarly named functions. The
goal is to remove unnecessary duplication without changing physical observables.

| Overlap | Decision | Rationale |
| --- | --- | --- |
| Minimum-image displacement in `core.geometry`, `core.lipid_geometry`, `mechanistic_common`, and `extended.common` | Consolidated `AnchorTracker.positions` and single-origin MSD steps onto `minimum_image_displacements`; retain narrow wrappers | Triclinic-aware MDAnalysis vectors; caller-specific dtype and output shape matter |
| `biophysical_common.lateral_distance` versus `core.lipid_geometry.planar_offsets` / `mechanistic_common.planar_distances` | **Retain separately** | The former rejects triclinic cells and computes an all-pairs matrix; latter methods compute one-to-one projected vectors after minimum imaging |
| `extended.common.centers`, `core.lipid_geometry.AnchorTracker` and `observables.msd._residue_centers` | **Retain separately** | Wrapped arithmetic COG for anchor selections, periodic unwrapped lipid anchor centers, and mass-weighted MSD COM with missing-mass fallback are intentionally different |
| `biophysical_common.leaflet_groups`, `core.lipid_geometry.classify_leaflets`, `extended.common.load_lipid_anchors` | **Retain separately** | Circular gap, periodic leaflet assignment and median split are not interchangeable |
| `biophysical_common.heavy`, `mechanistic_common.heavy`, `extended.common.heavy` | Shared implementation already exists in `core.atom_selection.select_heavy_atoms` | Distinct named policies protect name- vs element-based heavy-atom choices |
| Sparse contacts in `biophysical_common.atom_contacts`, `mechanistic_common.contact_count`, `extended.common.pairs` | **Retain adapters** | All use capped (sparse) contact searches but differ in coordinate precision, optional self-exclusion and return type (pairs vs count). No dense all-pairs replacement |
| `mechanistic_common.sequence_matches`, `biophysical_common.sequence_motifs`, `extended.common.sequence_hits` | **Retain separately** until output-contract migration is approved | Different regex overlap behavior, chain identification and output record types; unifying naively changes motif counts |
| Principal axes in `mechanistic_common.orientation_vector` and `extended.common.axis_from_ca` | **Retain separately** | Minimum atom requirements (3 vs 4 CA), atom selection, null return, and orientation normalizations differ |

### Current validation

Existing numerical-contract tests cover triclinic/orthorhombic PBC, mass fallbacks,
reordered residue indices, sparse contacts, multi-residue glycolipid anchors,
chain-local/overlapping motifs and empty selection behavior. These are regression
contracts, **not proof** of equivalence between the intentionally distinct
methods.

### Rules for subsequent consolidation

Before replacing a numerical primitive, document each caller, output shape,
units, periodic-box convention, dtypes and selection semantics. Compare
numerical outputs on both synthetic fixtures and opt-in real GROMACS systems.
Do not rewrite all-pairs calculations as pointwise projections (or vice versa),
or replace sparse search with dense matrices. Preserve CLI/YAML/API contracts.

### Completed code consolidation

`AnchorTracker.positions` and `observables.msd.msd_metrics` now use the
existing `core.geometry.minimum_image_displacements` adapter rather than
importing `MDAnalysis.lib.distances.minimize_vectors` separately. Both preserve
input displacement vectors, box dimensions, precision and result shapes;
the shared adapter adds explicit periodic box validation. Numerical regression
coverage is provided by existing numerical-contract PBC and MSD tests. No user-facing
migration is intended.

---

## Scientific numerical validation ledger

These deterministic tests are **independent numerical references**, not merely
round-trip execution checks. Run with:

```bash
pytest -q tests/test_publication_references.py
```

The independent reference suite covers direct, finite-sample-normalized FFT
autocorrelation; a zero-crossing automatic block-length case; exact independent
label permutations; paired sign-flip enumeration; Benjamini–Hochberg correction
with missing values and family separation; constant-data block and paired
hierarchical bootstrap controls; and direct, quadratic multi-origin lateral
MSD. Inputs are synthetic and contain no private simulation data.

Bootstrap uses NumPy Generator streams with a fixed seed; repeatability within
a fixed software environment is tested. Numerical results, rather than a
bitwise sequence of draws, should be compared across NumPy/Python versions
unless an RNG bit-generator and version are explicitly pinned.

Independent simulations (replicates), not autocorrelated frames, determine
between-condition test units. Paired tests preserve matching replicate labels.
Units: MSD coordinates are unwrapped Å; lateral MSD is Å²; ACF and corrected
p-values are dimensionless; block lengths are measured in sampled frames.

Outstanding before a publication-quality release: a redistributable YAML
end-to-end reference dataset with audited output tables and tolerances, mixed
observer interrupted-checkpoint restart comparison, profile-axis alignment
and missing-observation reference cases, representative real-GROMACS numerical
comparisons, benchmark resource budget, and final documented license/citation.
Passing these unit tests does not certify the remaining requirements.


## Matched local lipid populations and sampled residence (PR #117)

For species s, region R and instantaneous leaflet L, let A_s be the number
of molecular anchors within the configured xy disk and C_s the number of
those *same molecules* making at least one atom contact with R. Then
`matched_local_contact_fraction = C_s/A_s` and
`matched_local_composition_enrichment = (C_s/sum C)/(A_s/sum A)`.
Each molecule is counted once per species/region/frame; moiety subclasses
are not added into composition totals. Absent species give NaN. If there
are no contacts in any species, composition enrichment is undefined, while
contact fraction remains zero for available species. Means of enrichment
therefore condition on frames with contacts, and need their valid-frame
counts for interpretation. The disk is defined by anchor position, not by
headgroup reach or protein-excluded accessible area. Different molecular
sizes can still affect contact probabilities; enrichment is not affinity.

The whole-leaflet `contact_fraction` and `composition_enrichment` retain their
original populations. `local_composition_enrichment` and its project summary
`mean_local_enrichment` retain the legacy unmatched definition for comparison
only. Corrected project summaries are `mean_matched_local_contact_fraction`
and `mean_matched_local_enrichment`. Old persisted tables lacking corrected
fields produce missing corrected summaries, never legacy substitutions.
Figure inventory treats legacy local enrichment as supplementary. Corrected
figures and replica-averaged numerical summaries enter the PDF and bundle.

The shared numerical routine is `core/residence.py:sampled_episode_durations`.
For observed contact bounds i,j, span = t_j - t_i. The interval-imputed
estimate adds t_(j+1)-t_j, or the median observation interval at the final
sample. A one-observation trajectory has zero observed span and undefined
imputed duration. Lipid `duration_ps` retains its imputed convention;
`observed_span_ps` is explicit. Solvation `duration_ps` retains observed span;
`sampled_duration_ps` is explicit. The `duration_convention` column identifies
the imputed convention; it does not redefine solvation's legacy duration.

Lipid `occupancy_fraction` is the molecule's contact-frame fraction over the
whole observed window, repeated for each episode. `episode_contact_fraction`
is that episode's contact-frame count divided by the whole-window frame
count; it is **not** contact fraction within episode bounds. Neither is a
time-weighted occupancy under irregular sampling. Intermittent episodes
bridge at most `gap_frames` absent observations; gap time is included in the
span, but absent observations do not enter occupancy. Boundary censoring
remains explicit. These are sampled-contact estimators, not observed physical
on/off times. Existing solvation survival uses observed-span durations,
excludes unknown left origins and handles right censoring; interval censoring
between trajectory samples remains unresolved.

Time validation and terminal-interval calculation run once per episode table,
rather than once per event. Benchmark on Python 3.12, 2,000 samples and 200
solvents with 500 singleton episodes each (100,000 episodes): solvation episode
construction decreased from 4.714 s to 0.815 s (5.79x) with exact numerical
agreement. Python allocation peak increased from 47.30 to 48.58 MiB;
tracemalloc is not total process RSS. This is a synthetic postprocessing
benchmark, not a trajectory throughput prediction.

The passing release-backed artifact at commit `1110e89`, workflow
`41497488356`, contains four mechanism frames per topology. For `moiety=any`,
3/432 glycolipid-system rows and 4/384 control rows have finite enrichment
changes; maximum absolute differences are 0.3047619 and 0.1805556. The changed
anchor populations involve POPC and cholesterol respectively and alter other
species' contact shares too. These four-frame samples do not establish
production effect sizes, convergence or biological significance. The expanded
real smoke test writes `matched_local_impact.json` alongside mechanism tables.

Restart validation compares corrected specificity and both lipid residence
tables after interruption and clean execution. Use a fresh output directory when
recomputing historical results with the corrected estimator. An old output may
remain structurally readable even when it predates a scientific correction.

## Streamed solvent episodes and descriptive bridge references

The persisted contact-event writer orders observations by frame. The episode
extractor streams those records and uses the full observation-frame map; an
absent event is a sampled noncontact, not an omitted time point. It closes
contacts after missed observations, preserves duplicates without double-counting,
and uses the shared interval rule for terminal imputation. A single observation
has zero observed span and unavailable imputed duration. Episode output order
is completion order; scientific identity is in solvent/target/entity/frame keys.
The legacy in-memory `contact_episodes` reference remains useful for arbitrary
unordered inputs and independent parity tests. Episode writes are individually
atomic; this is not a multi-product transaction or a global memory guarantee.

Survival risk counts are computed from grouped duration multiplicities rather
than rescanning every episode for every distinct duration. Left-censored origins
are still excluded; right censoring and ties retain the product-limit convention.

Configured bridge pairs receive complete sampled-frame traces and deterministic
nonzero shift references. Same-entity identity and all-frame denominators are
preserved. At most 100 evenly spaced nonzero offsets are evaluated; small windows
use all offsets. Fewer than four samples or irregular sampling suppresses the
reference. Shift percentiles are descriptive, not Monte Carlo error bounds,
bootstrap confidence intervals, spatial controls or cooperativity p-values.
The existing replica statistics consume occupancy/count traces without treating
individual ions or frames as independent condition-level replicas. Raw bridge
counts are checked against target-contact reconstruction using integer counts.
Definitions and caveats are in [analyses](../analyses.md).

## Required companion products and offline completion (2026-11-08)

The existing registry's `analysis_output_plan` extends fused output planning
to standalone DSSP/SASA scheduling. `required_artifacts` lists required
companions for lipid matrices/enrichment, membrane-order mappings/index,
H-bond counts/audits, residue-distance axes/mappings/profiles, solvation
RDF/coordination/membership metadata, and SDC4 maps/audits/conditioned tables.
Configuration-dependent SASA atoms/burial, motif residence, solvation
kinetics/grid/reference/pair summaries are required only when requested.
Sparse bridge/event tables that are legitimately absent when no events occur
are not required. Disabled analysis families are not required.

`run_replicate` checks existing products before any reader/fallback can
replace them, and validates finalized products before summary completion.
`finalize_project` checks the same plan plus enabled multi-origin MSD tables
before archiving prior project statistics. Missing/corrupt companions produce
an explicit non-destructive failure with restoration/separate-output guidance;
missing primary products retain the established restart/fallback behavior.
No code/config/content digest determines reuse. Compatible historical markers
and data remain eligible regardless of source version.

Checks include Parquet footer/schema, required primary columns, valid JSON
objects, memory-mapped NPY readability, streamed ZIP CRC integrity for NPZ,
indexed membrane-order files, and residue-distance matrix/time/residue shape
and time-axis agreement. These checks do not inspect every Parquet data page,
prove scientific selection correctness, or infer historical analysis settings.
Existing streaming committed-chunk/replica-identity checks are preserved.
Core scalar summaries still use their established time/checkpoint contract;
this is not an exhaustive new schema/manifest framework.

Aqueous-region references and explicit pair occupancy use the definitions and
limits in [analyses](../analyses.md). Both are optional; the default full-cell
RDF, sampled residence estimators, and statistical RNG behavior are unchanged.
The report evidence overview counts finite scalar contrasts by the existing
category and q-based support policy. It neither ranks effects by interval
sign nor treats the category counts as independent discoveries.

## Periodic projected protein footprint and lipid area references

`core/lipid_geometry.py:accessible_area_grid` owns the existing SDC4 slab/disk
area estimator; mechanism and general lipid-enrichment callers share it. The
nearest-disk query is periodic in both lateral coordinates. Previously a protein
atom and a grid point on opposite sides of the cell could appear far apart in
a nonperiodic KD-tree despite being periodic neighbors. A hand-counted two-point
reference fails against that original implementation (1 Å² instead of 0 Å²).
The correction may change previously underestimated occlusion at periodic edges;
ordinary compact local grids away from this boundary case keep their results.

The optional general area adapter evaluates full-plane and annular equal-weight
quadratures using this same kernel. Annular midpoint sampling is uniform in
r² and angle. Its density ratio is distinct from composition enrichment, with
raw counts and both area denominators persisted. Duplicate obstacles cannot
subtract area twice; atoms outside the normal slab cannot occlude a leaflet.
Absent species and a fully occluded shell remain unavailable ratios. These
operational planar areas are not solvent-accessible surfaces, molecular free
areas or thermodynamic binding constants. Numerical resolution sensitivity is
separate from trajectory/statistical uncertainty.

## Cross-observer planar geometry contract

`core/lipid_geometry.planar_box` guards Cartesian-height leaflet assignment.
The two box angles involving the configured normal must be 90 degrees;
in-plane shear is permitted. Normals, enrichment and molecular lipid contacts
validate each current frame, not just initialization. Per-leaflet membrane
order validates the initial frame used to assign its fixed cohorts. General
3D minimum-image distances, periodic centers and MSD retain triclinic support;
this guard does not redefine those kernels. SDC4 height maps and projected
protein-excluded areas remain orthorhombic-only.

The configured axis is a laboratory reference. Rotation about it preserves
scalar distances, enrichment and normal-deviation angles. Rotating a tilted
membrane relative to that fixed axis changes its deviation by definition;
there is no general rotation-invariance claim for height fields. Local tilt
and weak curvature are tested with analytic planes and quadratic surfaces;
vesicles/overhangs require a different surface model.

`tests/test_planar_membrane_contracts.py` checks shared orthorhombic and
in-plane-skew trajectories under translation/independent lattice wrapping,
a hand-counted square contact graph, zero flat-plane tilt, leaflet exchange
and reversed reference labels, known tilted/quadratic surfaces and local
thickness, sparse-neighbor NaNs, initial-cohort MSD against known displacements,
and unsupported boxes at initialization and after a valid frame. Existing
coordinate/dual-MSD tests protect source versus fitted provenance and timestamp
synchronization. Legacy membrane APL/thickness/MSD cohorts still use the
initial median partition; their centered/resolved first-frame requirement is
unchanged and is not an instantaneous circular-leaflet guarantee.
