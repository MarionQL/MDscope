# Issue #31 — numerical primitive audit

This table records **semantic contracts**, not just similarly named functions. The
goal is to remove unnecessary duplication without changing physical observables.

| Overlap | Decision | Rationale |
| --- | --- | --- |
| Minimum-image displacement in `core.geometry`, `core.lipid_geometry`, `mechanistic_common`, and `extended.common` | Already centralized in `minimum_image_displacements`; retain narrow wrappers | Triclinic-aware MDAnalysis vectors; caller-specific dtype and output shape matter |
| `biophysical_common.lateral_distance` versus `core.lipid_geometry.planar_offsets` / `mechanistic_common.planar_distances` | **Retain separately** | The former rejects triclinic cells and computes an all-pairs matrix; latter methods compute one-to-one projected vectors after minimum imaging |
| `extended.common.centers`, `core.lipid_geometry.AnchorTracker` and `observables.msd._residue_centers` | **Retain separately** | Wrapped arithmetic COG for anchor selections, periodic unwrapped lipid anchor centers, and mass-weighted MSD COM with missing-mass fallback are intentionally different |
| `biophysical_common.leaflet_groups`, `core.lipid_geometry.classify_leaflets`, `extended.common.load_lipid_anchors` | **Retain separately** | Circular gap, periodic leaflet assignment and median split are not interchangeable |
| `biophysical_common.heavy`, `mechanistic_common.heavy`, `extended.common.heavy` | Shared implementation already exists in `core.atom_selection.select_heavy_atoms` | Distinct named policies protect name- vs element-based heavy-atom choices |
| Sparse contacts in `biophysical_common.atom_contacts`, `mechanistic_common.contact_count`, `extended.common.pairs` | **Retain adapters** | All use capped (sparse) contact searches but differ in coordinate precision, optional self-exclusion and return type (pairs vs count). No dense all-pairs replacement |
| `mechanistic_common.sequence_matches`, `biophysical_common.sequence_motifs`, `extended.common.sequence_hits` | **Retain separately** until output-contract migration is approved | Different regex overlap behavior, chain identification and output record types; unifying naively changes motif counts |
| Principal axes in `mechanistic_common.orientation_vector` and `extended.common.axis_from_ca` | **Retain separately** | Minimum atom requirements (3 vs 4 CA), atom selection, null return, and orientation normalizations differ |

## Current validation

Existing issue-31 tests cover triclinic/orthorhombic PBC, mass fallbacks,
reordered residue indices, sparse contacts, multi-residue glycolipid anchors,
chain-local/overlapping motifs and empty selection behavior. These are regression
contracts, **not proof** of equivalence between the intentionally distinct
methods.

## Rules for subsequent consolidation

Before replacing a numerical primitive, document each caller, output shape,
units, periodic-box convention, dtypes and selection semantics. Compare
numerical outputs on both synthetic fixtures and opt-in real GROMACS systems.
Do not rewrite all-pairs calculations as pointwise projections (or vice versa),
or replace sparse search with dense matrices. Preserve CLI/YAML/API contracts.
