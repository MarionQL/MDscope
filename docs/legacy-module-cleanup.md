# MDscope cleanup policy and completed migration

MDscope permits breaking changes to historical Python import paths. Once a
canonical implementation replaces a compatibility module, migrate all internal
consumers and delete the old wrapper. No deprecation period is required.

## Completed in this PR

Deleted the six superseded re-export modules:
`analysis/alignment.py`, `analysis/lipid_geometry.py`,
`analysis/msd.py`, `analysis/bilayer_normal.py`,
`analysis/lipid_enrichment.py`, and `analysis/profiling.py`.

Production callers now import the respective canonical modules in
`analysis/core/` or `analysis/observables/`. Synthetic and real-trajectory
tests were migrated, and `tests/test_module_imports.py` checks the canonical
API rather than deleted compatibility exports. Obsolete shim-identity tests
were removed while substantive numerical tests were retained under core tests.

## Rules for future cleanup

1. One implementation per scientific operation; delete superseded wrappers.
2. Organize code by responsibility, not by its creation date.
3. Change imports, tests, launchers and packaging references together with moves.
4. Remove obsolete documentation rather than retaining contradictory guides.
5. Avoid speculative frameworks, extra indirection and parallel pipelines.
6. Preserve scientific numerical definitions (PBC, COM/COG, units, selections).
7. Prefer sparse/vectorized computation; preserve fused trajectory pass counts.
8. Retain independent numerical and end-to-end tests, not tests of deleted APIs.
9. Require Python 3.11/3.13 and relevant real GROMACS checks for structural PRs.
10. Measure simplification by fewer concepts and less duplication, not just LOC.
11. Keep generated caches and large trajectories out of committed source.
12. Delete only after repository consumers are migrated in the same branch.

The cleanup intentionally breaks historical Python module import paths.
No numerical outputs, CLI commands, or YAML schema are intentionally changed.
