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

In `analysis/streaming.py`, fused membrane, extensions, densities and
membrane metrics are evaluated **before** fitted structural arrays are
computed. The separate MSD Universe is frame-synchronized in
`core/trajectory_sources.py`. Alignment constructs derived arrays; it does
not intentionally rotate the main Universe.

`core/coordinate_contracts.py` supplies explicit labels and a small
guard for programmatic handoffs. It does not attempt to infer PBC or
preprocessing history from coordinates. That would require additional input
metadata and validation and should be addressed separately.

### Future verification

For publication-quality diffusion, validate the processing provenance of
input trajectories, fit-window sensitivity, sampling versus box displacement,
and independent reference trajectories. Do not treat the `multi_origin_msd`
numerical kernel as wired into production until its unwrapped input adapter
and storage policy are implemented.
