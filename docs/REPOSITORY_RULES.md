# MDscope repository rules

This document is the single source of truth for structural refactoring.

1. **No backward compatibility for superseded internal modules.** Migrate consumers and delete old paths in the same PR.
2. **Minimize code and complexity.** Remove redundant functions, wrappers, unused configuration and unnecessary classes.
3. **Organize by responsibility.** Keep geometry, selections and trajectory infrastructure under `core/`; scientific metrics under `observables/`; molecular contacts under `interactions/`; system-specific mechanisms under `mechanisms/`; inference and comparisons under `statistics/`; figures and exports under `reporting/`; streaming and orchestration under `pipeline/`.
4. **One documentation source of truth.** Delete or rewrite stale guides rather than retaining contradictory instructions.
5. **One implementation per operation.** Distinct scientific semantics must stay distinct; do not consolidate merely because names resemble one another.
6. **Preserve scientific correctness.** A structural move must not silently alter algorithms, units, mass weighting, PBC handling, atom selection, frames or outputs.
7. **Tests follow current architecture.** Update current tests, remove tests for deleted APIs and retain valuable numerical/integration coverage. Avoid adding tests just to check wrappers.
8. **No speculative abstractions.** Introduce helpers/classes only when they remove real duplication.
9. **Minimize trajectory passes, memory and I/O.** Preserve fused processing, streaming and sparse distance calculations.
10. **Separate execution from analysis.** HPC orchestration, trajectory processing and scientific statistics must have clear boundaries.
11. **Require meaningful CI.** Run existing Python 3.11/3.13 checks and real-GROMACS workflows for changed paths.
12. **Prefer net simplification.** Explain what a cleanup PR removes, and avoid mass formatting.

## Target layout

- `analysis/core/`: geometry, coordinate and topology contracts, selections, trajectory inputs, persistence and checkpoints
- `analysis/observables/`: RMSD, RMSF, MSD, radius of gyration, membrane and protein observables
- `analysis/interactions/`: hydrogen bonds, salt bridges and molecular contacts
- `analysis/mechanisms/`: SDC4, glycan/lipid mechanisms and topology-specific analyses
- `analysis/statistics/`: resampling, statistical summaries and system comparisons
- `analysis/reporting/`: plots, tables and PDF/report exports
- `analysis/pipeline/`: trajectory streaming, fused observers and run/finalize orchestration
- `analysis/cli.py`, `analysis/config.py`, `analysis/registry.py`: minimal top-level entry points

Move one coherent group at a time, migrate all imports, and confirm CI before proceeding. These are ownership boundaries, not a requirement to invent additional wrappers.
