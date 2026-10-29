# Opt-in multi-time-origin lipid diffusion

Enable in `analyses.membrane.msd`:

```yaml
msd:
  enabled: true
  multi_origin: true
  multi_origin_max_lag_frames: 1000  # null = all sampled lags
  multi_origin_particle_batch: 64
  slope_start_fraction: 0.2
  slope_end_fraction: 0.8
```

Runs a **separate trajectory pass** per replicate after the standard streaming
stage. This is deliberately opt-in and preserves the existing single-origin
timeseries and diffusion estimates. The optional `msd_trajectory` is used if
specified; otherwise the main trajectory is used. Neither can be assumed
unmodified merely from filenames.

For each selected lipid, the reader accumulates minimum-image frame-to-frame
center displacements (and optionally subtracts common drift), then writes
these into a temporary on-disk float64 NumPy memmap. The FFT MSD is computed
over small batches of lipids and averaged across lipids and all available time
origins. Disk usage is approximately
`24 * sampled_frames * selected_lipids` bytes, plus temporary FFT memory
bounded by `multi_origin_particle_batch` (the final tables are much smaller).
Temporary coordinates are deleted when processing finishes or errors.

The maximum lag is in **sampled frames** after trajectory stride. The fitting
window is expressed as fractions of maximum retained lag time, excluding lag
zero, and must contain at least four lag points. Sampling must be regular.
Output units are Å² for MSD and nm²/ns for D, using the 2D Einstein relation
`D = slope/4`. Reported D values are descriptive slopes, not uncertainty
estimates; independent replicates remain the inference units.

Replicate outputs:
- `summary/multi_origin_msd.parquet`: lag time, MSD, number of time origins
- `summary/multi_origin_diffusion.parquet`: per-replicate fit by all/upper/lower lipids

Project outputs:
- `project_summary/diffusion/multi_origin_replicates.parquet`
- `project_summary/diffusion/multi_origin_systems.parquet`

**Limits:** Minimum-image reconstruction assumes sampled displacements are less
than half the relevant box dimension in an unambiguous image sense. The FFT
kernel assumes stable lipid identity, a fixed analysis frame, and that the
appropriate raw trajectory has been supplied. Long-lag MSD estimates have
fewer time origins and correlated lag values; high R² is not evidence that the
diffusive regime was reached. The second pass increases trajectory I/O but
does not alter standard streaming checkpoints.
