# MDscope architecture

See [REPOSITORY_RULES.md](REPOSITORY_RULES.md) for cleanup and refactoring
policy. This document records the current canonical module layout and dependency boundaries.

## Current canonical locations

- `analysis/core/`: alignment, trajectory and provenance helpers, periodic
  geometry, atom selection, checkpoints, persistence.
- `analysis/observables/`: RMSD, RMSF, radius of gyration, MSD, bilayer
  normal, lipid enrichment, membrane metrics/order, density, DSSP, SASA,
  structural profiles and residue–membrane distance.
- `analysis/mechanisms/`: SDC4 conformational mechanisms, HS topology,
  mechanistic geometry, local environment and shared conformational states.
- `analysis/statistics/`: replicate-aware inference, resampling, time-series, profile and interaction summaries.
- `analysis/reporting/`: report PDF/ZIP generation, membrane and SDC4 figures,
  coupling-screen reporting and plots.
- `analysis/cli.py`, `analysis/config.py`, `analysis/registry.py`: current
  CLI, configuration and registry entry points.

- `analysis/interactions/`: hydrogen bonds, salt bridges, ligand contacts,
  lipid hydrogen bonds and shared interaction tables.
- `analysis/pipeline/`: streaming, fused trajectory processing, orchestration,
  finalization and execution progress.

The former flat runner, streaming, progress, profiles and timeseries modules
have been removed; import implementations directly from their owning packages.

The repository intentionally does **not** preserve legacy Python import
aliases after a replacement is migrated. `mdscope` CLI commands and YAML
configuration stay unchanged unless a separate change is approved.

## Dependency direction

Core numerical kernels and observables must not import reporting. Pipeline
coordinates readers and observers; statistics operates on persisted numerical
results; reporting consumes finalized tables and images without rerunning
trajectories. Keep framewise computations fused and sparse where applicable.

## Quality gates

Use the existing Python 3.11/3.13 tests, wheel and CLI checks, and real
GROMACS smoke workflow. Structural refactors must preserve units, PBC
handling, atom and residue identity, output schemas and replicate indexing.
No new regression tests are required merely for relocating modules.
Verify the package from outside the source checkout as well as in-tree, and
check the real-GROMACS trajectory workflow before treating a refactor as
validated. GitHub check names correspond to `.github/workflows/mdscope-checks.yml`
and `.github/workflows/real-trajectory-smoke.yml`.

A passing smoke test establishes coverage for its tested configurations,
not numerical equivalence for every supported topology or analysis option.

## Issue #33 — observable and interaction ownership audit (2026-10-31)

The former `analysis/timeseries.py`, `analysis/profiles.py` and
`analysis/interactions.py` facades were intentionally removed by earlier
architecture PRs; do not recreate compatibility wrappers for unsupported
imports. Numerical observables live under `observables/`, interaction
serializers and identity tables under `interactions/`, source/frame handling
under `core/trajectory_sources.py`, and the fused frame consumers under
`pipeline/`. `core/coordinate_contracts.py` documents provenance roles,
not absolute guarantees of input unwrapping. Pure kernels should be extracted
only when a real caller or numerical reference benefits; avoid a second
abstract observable framework.

Current tests cover direct FFT multi-origin MSD reference sums, independent
PBC geometry, topology-residue contacts, mixed protein–membrane and
membrane-only selections, explicit coordinate roles, and fused frame
lifecycle. Broader numerical, topology, and trajectory provenance validation
belongs to issue #36 rather than repeated architecture reorganization.
