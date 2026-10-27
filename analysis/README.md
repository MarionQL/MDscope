# Analysis source architecture

The public API remains `MDscope.analysis.*` (also importable as
`mdscope.analysis.*` in editable/symlink installs). Modules retain their
legacy entrypoints so project YAML, SLURM invocation and external scripts do
not break during incremental organization.

- `core/`: execution infrastructure and optional per-stage profiler.
- `observables/`: physics/trajectory algorithms (MSD, lipid enrichment, local bilayer normals).
- `reporting/`: replicate/project-level postprocessing (coupling screen).
- `extended/`: shared geometry for the extended membrane observers.
- Root `cli.py`, `runner.py`, `streaming.py`: command entrypoint, orchestration,
  one-pass frame iteration and checkpoints.
- Remaining root modules: original public compatibility surface; they can
  migrate incrementally only with scientific numerical regression tests and
  import compatibility checks.

Raw-coordinate membrane observers execute before protein alignment in
`streaming.py`. Never rotate or center the whole trajectory inside the
shared frame loop: MSD depends on physical lipid displacements. Aligned
protein feature coordinates are derived in separate arrays.

New code should import concrete implementations from their grouped packages.
Legacy imports (for example `analysis.msd`) are intentionally supported
through lightweight compatibility wrappers. `tests/test_module_imports.py`
asserts that both paths export the same implementation objects.
