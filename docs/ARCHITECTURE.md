# MDscope architecture and compatibility policy

## Current package boundaries

`mdscope.analysis.cli`: command-line user interface.
`mdscope.analysis.runner`: replicate/project orchestration and output merging.
`mdscope.analysis.streaming`: shared frame pass, checkpoint state and data capture.
`mdscope.analysis.core`: trajectory provenance, rigid-body alignment, reusable lipid/PBC geometry and runtime infrastructure. The canonical implementations are `core/alignment.py` and `core/lipid_geometry.py`; legacy flat imports remain compatibility re-exports.
`mdscope.analysis.observables`: independently testable numerical observables.
`mdscope.analysis.extended`: supporting extended observables.
`mdscope.analysis.registry`: specialized analysis registration and output contracts.
`mdscope.analysis.statistics`: statistical postprocessing.
`mdscope.analysis.reporting` and `report`/`plotting`: presentation output.

Not all legacy implementations have been moved into those boundaries yet. Flat
`mdscope.analysis.*` imports intentionally remain supported. Avoid moving
modules merely for symmetry; use small moves with regression coverage.

## Dependency direction

- Core numerical kernels should depend on NumPy and their explicit inputs,
  not CLI configuration, run directories, or report generation.
- Trajectory adapters translate MDAnalysis AtomGroups/frames into the numeric
  inputs and preserve units, atom identity, timestamps and periodic coordinates.
- Pipeline orchestration calls adapters and kernels; it owns checkpointing and
  output paths rather than individual scientific algorithms.
- Statistics operates on validated outputs, preserving frame/replicate identity.
- Reporting never changes numerical measurements.

## Publication-oriented quality gates

1. Independently test each important numerical kernel on analytic and
   reference-implementation cases, including PBC and edge cases.
2. Test imports, CLI entry points and wheel installation in a clean environment.
3. Run representative real-trajectory workflows and verify output schemas.
4. Check checkpoint restart behavior and source/config provenance.
5. Explicitly document units, selection policies, time origins, estimators,
   statistical uncertainty and known limitations.

## Packaging

A distribution called `mdscope` installs the import namespace
`mdscope.analysis` and console command `mdscope`. The physical repository
directory is still named `MDscope`; existing HPC commands such as
`python -m MDscope.analysis.cli` continue to work when the repository's
parent directory is on `PYTHONPATH`. Prefer the installed lower-case command
for new deployments, and retain older submission scripts until their
environments have been tested.

Only the standard runtime dependencies are mandatory. Report generation and
some specialized scientific analyses may require optional dependencies.
Do not change user-facing YAML keys, public imports or output schemas as part
of purely structural refactors.
