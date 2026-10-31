# Trajectory and coordinate contracts

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

### Execution boundary

In `analysis/pipeline/streaming.py`, fused membrane, extensions, densities and
membrane metrics are evaluated **before** fitted structural arrays are
computed. The separate MSD Universe is frame-synchronized in
`core/trajectory_sources.py`. Alignment constructs derived arrays; it does
not intentionally rotate the main Universe.

`core/coordinate_contracts.py` supplies explicit labels and a small
guard for programmatic handoffs. It does not attempt to infer PBC or
preprocessing history from coordinates. That would require additional input
metadata and validation and should be addressed separately.

### Single-origin and multi-origin MSD

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
[MULTI_ORIGIN_MSD.md](MULTI_ORIGIN_MSD.md).

For publication-quality diffusion, validate the processing provenance of
input trajectories, fit-window sensitivity, sampling versus box displacement,
and independent reference trajectories. Neither an excellent linear-fit R²
nor a clean trajectory pass establishes the diffusive regime.

### Contact and residue identity

`biophysical_common.atom_contacts` uses MDAnalysis capped-distance with
an optional periodic box; distances and cutoffs are in Å. It returns indices
**local to the two selected AtomGroups**, not global atom IDs. The pure NumPy
`collapse_atom_pairs_to_residues` function maps these indices to unique,
sorted **topology residue indices**; it must not collapse by the human-readable
`resid`, which can repeat between protein chains or lipid segments.
`residue_pairs` combines these two operations. This contract applies to
protein–membrane and membrane-only selections alike. A `None` box is
nonperiodic; callers must pass frame dimensions when PBC is required.

### Remaining scientific checks

The standalone multi-origin adapter uses a second reader by design; do not
claim all optional observables are single-pass. Chemical selections and mass
fallback policies remain explicit; drift removal for both MSD estimators is
optional and must be recorded for each simulation.
