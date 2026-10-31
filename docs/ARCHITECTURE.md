# MDscope architecture

See [REPOSITORY_RULES.md](REPOSITORY_RULES.md) for cleanup and refactoring
policy. This document records the **current** layout and the planned package
boundaries; do not assume that every planned package already exists.

## Current canonical locations

- `analysis/core/`: alignment, trajectory and provenance helpers, periodic
  geometry, atom selection, checkpoints, persistence.
- `analysis/observables/`: RMSD, RMSF, radius of gyration, MSD, bilayer
  normal and lipid enrichment.
- `analysis/reporting/`: report PDF/ZIP generation, membrane and SDC4 figures,
  and coupling-screen reporting.
- `analysis/cli.py`, `analysis/config.py`, `analysis/registry.py`: current
  CLI, configuration and registry entry points.

## Remaining migration

- `analysis/interactions/`: hydrogen bonds, salt bridges, ligand contacts,
  lipid hydrogen bonds and shared interaction tables (now grouped).
- `analysis/mechanisms/`: organize SDC4 mechanisms and system-specific
  topology/mechanistic calculations.
- `analysis/statistics/`: move statistical postprocessing without changing
  estimators or aggregation.
- `analysis/pipeline/`: move streaming, fused passes, runner and finalization
  after updating all dependent imports and HPC entry points.
- Other physical observables and plotting utilities still have flat modules;
  move them by coherent family and remove the originals.

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
