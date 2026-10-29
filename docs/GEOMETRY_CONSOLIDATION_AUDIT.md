# Issue #31 — Scientific geometry consolidation audit

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

## Further work before closing issue

Add topology-backed tests for absent masses/elements, reordered residue IDs,
selections spanning multi-residue glycolipids, empty groups and motif
boundaries. Compare existing real-trajectory outputs before modifying other
numerical implementations. No implicit change to COM/COG, axis conventions,
leaflet estimators, or contact precision is authorized by this PR.
