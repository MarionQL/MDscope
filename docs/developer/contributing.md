# Contributing

**Rule:** one implementation per operation; preserve scientifically different contracts even when code looks similar. Prefer net simplification.

## Where code belongs

| Module | Owns |
|---|---|
| `core/` | Geometry, atom selections, trajectory sources, persistence, checkpoints |
| `observables/` | RMSD/RMSF/MSD, membrane/lipid/structural observables |
| `interactions/` | Molecular contacts and H bonds |
| `mechanisms/` | SDC4 and other system-specific analyses |
| `statistics/` | Replicate-level inference and comparisons |
| `reporting/` | Findings, scientific figures and exports |
| `pipeline/` | Shared streaming, orchestration, finalization |
| `cli.py`, `config.py`, `registry.py` | CLI, YAML and task discovery |

## Change checklist

- **Scope:** inspect callers, relevant contracts, tests and docs. Change one coherent responsibility at a time; avoid speculative abstractions.
- **Compatibility:** migrate consumers and remove superseded *internal* modules in the same PR; preserve public CLI/YAML/output contracts unless explicitly changed.
- **Scientific invariants:** retain units, populations, coordinate frame, PBC, mass weighting, selection chemistry, time sampling, missing-vs-zero semantics and results.
- **Performance:** minimize trajectory traversals, memory, I/O and dense atom-pair materialization. Preserve fused and sparse processing.
- **Testing:** use independent numerical oracles for mathematical claims; check integration through the actual entry point. A successful run or two agreeing implementations may share the same bug.
- **Failure safety:** test invalid inputs and checkpoint corruption without destroying valid existing outputs. Verify artifact identity, integrity, completeness and provenance—not just existence.
- **Validation:** run relevant targeted tests, Python 3.11/3.13 CI, and real GROMACS checks when affected. Document what ran and what remains unverified.
- **Documentation:** one authoritative location per concept. Prefer short tables and bullets; remove stale guides and avoid mass formatting.

See [Scientific contracts](scientific-contracts.md), [Scientific audit](scientific-audit.md) and [Validation gaps](validation-gaps.md) for the detailed acceptance criteria.

**AI-assisted changes:** review behavior and scientific evidence, not the origin of the code. For substantive fixes, demonstrate the prior failure or a targeted fault-injection/mutation detection where practical. A green CI run is not evidence of physical convergence.
