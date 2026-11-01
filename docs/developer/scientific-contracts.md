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
