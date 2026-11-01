# Solvation hot-path audit and benchmark (2026-11-01)

This reviews the existing merged PR #78 implementation, rather than
introducing an additional analysis framework.

## Code footprint

PR #78: **1,608 added lines in 11 files** (1 deletion), including 298 lines
of documentation, 196 lines of benchmark scripts, and 355 lines of tests.
The numerical/postprocessing implementation was 505 + 203 = 708 lines;
the remaining lines were configuration and fused-pipeline wiring. Removing
documentation or independent tests would reduce line count without eliminating
runtime overhead. No new abstraction layer is justified.

## Locally measured synthetic benchmark

Measured with MDAnalysis 2.10.0 and NumPy, 150 Å cubic periodic box,
2,400 selected atoms grouped into 800 targets, 80,000 water oxygens,
8 Å RDF search and 3.5 Å coordination shell. Uniform random coordinates;
**not** a biomolecular trajectory or a validation of physical hydration.

| Operation | Baseline | Optimized | Notes |
| --- | ---: | ---: | --- |
| Periodic search, 8 Å | 85.97 ms | unchanged | Median over 5 invocations; 122,166 pairs |
| Deduplicate + count, 3.5 Å | 1.61 ms | 1.08 ms | Median over 150 iterations; 10,276 shell pairs |
| Neighbor search backend auto vs nsgrid | 123.1 ms | 98.9 ms | Different short benchmark, 4 repetitions; **not a verified production speedup** |
| Neighbor search backend pkdtree | 277.2 ms | — | Same short benchmark, 4 repetitions |

The optimized coordination counts match the original Python set algorithm
exactly, including repeated site–water pairs and unoccupied target IDs.
**The full workload gain will be small** because the periodic neighbor
search dominates this synthetic scenario. The code change is valuable because
it avoids per-pair Python objects entirely in the standard preset. It should
not be sold as a large end-to-end speedup.

## Code changes

- `unique_target_solvent_pairs` performs lexicographic sort, adjacent-pair
  deduplication and per-target bincount, retaining exact target/water identity
  when detailed analysis needs it.
- Standard mode **does not create** the large Python `dict[int, set[int]]`
  that detailed/kinetics mode requires for solvent identity output.
- Optional `neighbor_method` defaults to `auto` to preserve existing
  behavior. Alternatives are `nsgrid`, `pkdtree`, `bruteforce`.
- The existing *real* trajectory benchmark takes `--neighbor-method`;
  compare options before setting any HPC-wide default.

## HPC benchmark to perform before tuning further

On Tinkercliffs, with validated solvation target selections, run the same
representative trajectory window with each neighbor backend:

```bash
python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --neighbor-method auto \
  --output solvation-auto.json

python scripts/benchmark_solvation_real.py synd4_solvation.yaml \
  --system synd4_glycolipids --rep rep1 --frames 100 \
  --preset standard --neighbor-method nsgrid \
  --output solvation-nsgrid.json
```

Repeat with `--preset kinetics` and an explicit
`options.track_target_classes` for interesting protein/HS sites, then
repeat in the nonglycolipid condition. Compare **total observer throughput**,
maximum RSS and bytes of output, not just the native PBC search speed.

## Remaining performance risks

1. Full-system water identity tracking and all-target bridge enumeration
   can still produce enormous tables; use the built-in class filter.
2. `finalize()` concatenates all identity-event chunks in memory before
   writing outputs, which may become limiting on very long trajectories.
   A subsequent bounded-memory streaming finalizer is worthwhile **only
   after measuring the affected realistic kinetics preset**.
3. The radial search radius determines a large fraction of candidate-pair
   work. Preserve a radius sufficient for RDF structure/3D targets;
   do not shorten it purely to improve benchmark speed.
4. The current benchmark does not profile other fused MDscope observers
   (SASA, DSSP, protein/lipid mechanisms) running concurrently.
5. Profile a representative multi-thousand-frame run for queue-time
   extrapolation and space usage before comparing complete 2 µs replicates.

**No production TPR/XTC was accessible in this environment; performance
numbers are synthetic and do not establish real MDscope throughput.**

## Shared PBC implementation

The review also found duplicate calls to MDAnalysis
`minimize_vectors` inside the solvation H-bond kernel and local 3D
water-density alignment. Both now use MDscope's existing
`core.geometry.minimum_image_displacements`, which preserves a single
triclinic minimum-image validation contract. A dedicated triclinic wrapped
H-bond regression was added. The sparse neighbor search still calls
MDAnalysis `capped_distance` directly because its atom-pair + distance
output is not equivalent to the existing `biophysical_common.atom_contacts`
(pair-only API); forcing one API onto the other would add complexity and
unnecessary distance computation for existing callers.

RDF histogram values at *exact* bin boundaries remain sensitive to
floating-point rounding between native distance backends. This is a
histogram numerical-convention issue, not a distinct PBC convention;
ordinary nonboundary backend parity is tested separately.
