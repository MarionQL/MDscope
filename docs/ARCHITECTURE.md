# MDscope architecture

See [REPOSITORY_RULES.md](REPOSITORY_RULES.md) for cleanup and refactoring
policy. This document records the **current** layout and the planned package
boundaries; do not assume that every planned package already exists.

## Current canonical locations

- `analysis/core/`: alignment, trajectory and provenance helpers, periodic
  geometry, atom selection, checkpoints, persistence.
- `analysis/observables/`: RMSD, RMSF, radius of gyration, MSD, bilayer
  normal, lipid enrichment, membrane metrics/order, density, DSSP, SASA,
  structural profiles and residue–membrane distance.
- `analysis/mechanisms/`: SDC4 conformational mechanisms, HS topology,
  mechanistic geometry, local environment and shared conformational states.
- `analysis/statistics/`: replicate-aware inference, resampling and summaries.
- `analysis/reporting/`: report PDF/ZIP generation, membrane and SDC4 figures,
  coupling-screen reporting and plots.
- `analysis/cli.py`, `analysis/config.py`, `analysis/registry.py`: current
  CLI, configuration and registry entry points.

## Remaining migration

- `analysis/interactions/`: hydrogen bonds, salt bridges, ligand contacts,
  lipid hydrogen bonds and shared interaction tables (now grouped).
- `analysis/pipeline/`: fused trajectory execution is grouped here; streaming,
  runner and finalization remain at top level pending a separate migration.

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
