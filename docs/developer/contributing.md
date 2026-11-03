# Contributing and repository rules

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

## Package ownership

- `src/mdscope/core/`: geometry, coordinate and topology contracts, selections, trajectory inputs, persistence and checkpoints
- `src/mdscope/observables/`: RMSD, RMSF, MSD, radius of gyration, membrane and protein observables
- `src/mdscope/interactions/`: hydrogen bonds, salt bridges and molecular contacts
- `src/mdscope/mechanisms/`: SDC4, glycan/lipid mechanisms and topology-specific analyses
- `src/mdscope/statistics/`: resampling, statistical summaries and system comparisons
- `src/mdscope/reporting/`: plots, tables and PDF/report exports
- `src/mdscope/pipeline/`: trajectory streaming, fused observers and run/finalize orchestration
- `src/mdscope/cli.py`, `src/mdscope/config.py`, `src/mdscope/registry.py`: minimal top-level entry points

Move one coherent group at a time, migrate all imports, and confirm CI before proceeding. These are ownership boundaries, not a requirement to invent additional wrappers.

## AI-assisted implementation and review

AI-generated changes are reviewed against the **behavioral and scientific contract**, not judged by whether the code runs or the agent reports success. Before editing, identify the intended observable, units and populations, input/outputs and provenance, mutable/persisted state, failure outcomes, and affected call sites. Keep the review proportional to the change: a typo does not need an adversarial test suite; an estimator, restart, or reporting change does.

For a substantive bug fix, first reproduce the failure with a small fixture or a previously failing case. Add a test that checks the required outcome **through the actual affected entry point**, not a tautology about the test fixture. Where practical, show that the test fails on the original behavior and passes with the correction (or use a targeted mutation/fault injection). Prefer analytical or independent reference calculations for scientific claims; comparing two paths that reuse the same implementation cannot establish their shared correctness. Avoid redundant tests that only duplicate implementation details.

At integration boundaries, check identity, shape, units, coordinate frame/PBC, time sampling, absent-versus-zero semantics, and intended missing-data behavior. For persisted products, distinguish existence, integrity, freshness, scientific completeness, and provenance. Decide explicitly whether a missing or corrupted required product must be rebuilt or cause a clear, non-destructive failure. Do not classify optional or deliberately disabled analyses as failures.

Before claiming completion, report the exact tests/checks executed, the behavior they establish, known assumptions and limitations, and any real-trajectory or production-scale validation that remains unperformed. A green CI run, matching hashes, successful PDF export, or a plausible plot is not an independent scientific validation. See [scientific contracts](scientific-contracts.md), [scientific audit](scientific-audit.md), and [validation effectiveness/gaps](validation-gaps.md) rather than restating their domain-specific caveats here.
