# Geometry helper audit: shared operations versus distinct contracts

This audit compares the October 2026 implementations in
`biophysical_common.py`, `extended/common.py`,
`mechanistic_common.py`, and `lipid_geometry.py`.

## Consolidated now

The minimum-image transformation of *displacement vectors* is the same
MDAnalysis operation in both `extended.common.minimum_image` and
`mechanistic_common.planar_distances`. The implementations now call
`core.geometry.minimum_image_displacements`, preserving the extended
helper's explicit float32 coercion while keeping the mechanistic helper's
input precision. The shared function validates the box and coordinate shape.

## Deliberately separate

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

## Follow-up safety work

1. Test the three hydrogen selection rules against real CHARMM/TPR topologies before standardizing.
2. Choose an explicit contract for triclinic lateral geometry and distinguish physical projection from deleting an array coordinate.
3. Audit the exact callers of each leaflet classification before unifying schemes.
4. Never change selection, return shape, dimensionality, precision, or failure policy solely to deduplicate code.
