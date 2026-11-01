# Architecture and packaging

This guide consolidates the authoritative technical details formerly scattered across topic and implementation-history files. Sections retain their scientific caveats and validation limitations; consult the source code for the executable contract.


---

### MDscope architecture

See [REPOSITORY_RULES.md](../developer/contributing.md) for cleanup and refactoring
policy. This document records the current canonical module layout and dependency boundaries.

### Current canonical locations

- `src/mdscope/core/`: alignment, trajectory and provenance helpers, periodic
  geometry, atom selection, checkpoints, persistence.
- `src/mdscope/observables/`: RMSD, RMSF, radius of gyration, MSD, bilayer
  normal, lipid enrichment, membrane metrics/order, density, DSSP, SASA,
  structural profiles and residue–membrane distance.
- `src/mdscope/mechanisms/`: SDC4 conformational mechanisms, HS topology,
  mechanistic geometry, local environment and shared conformational states.
- `src/mdscope/statistics/`: replicate-aware inference, resampling, time-series, profile and interaction summaries.
- `src/mdscope/reporting/`: report PDF/ZIP generation, membrane and SDC4 figures,
  coupling-screen reporting and plots.
- `src/mdscope/cli.py`, `src/mdscope/config.py`, `src/mdscope/registry.py`: current
  CLI, configuration and registry entry points.

- `src/mdscope/interactions/`: hydrogen bonds, salt bridges, ligand contacts,
  lipid hydrogen bonds and shared interaction tables.
- `src/mdscope/pipeline/`: streaming, fused trajectory processing, orchestration,
  finalization and execution progress.

The former flat runner, streaming, progress, profiles and timeseries modules
have been removed; import implementations directly from their owning packages.

The repository intentionally does **not** preserve legacy Python import
aliases after a replacement is migrated. `mdscope` CLI commands and YAML
configuration stay unchanged unless a separate change is approved.

### Dependency direction

Core numerical kernels and observables must not import reporting. Pipeline
coordinates readers and observers; statistics operates on persisted numerical
results; reporting consumes finalized tables and images without rerunning
trajectories. Keep framewise computations fused and sparse where applicable.

### Quality gates

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

### Issue #33 — observable and interaction ownership audit (2026-10-31)

The former `src/mdscope/timeseries.py`, `src/mdscope/profiles.py` and
`src/mdscope/interactions.py` facades were intentionally removed by earlier
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

---

### Packaging and source layout

MDscope uses a conventional `src/mdscope/` package layout.
The public imports remain `mdscope.*` and the CLI remains `mdscope`.

Install from a checkout with `python -m pip install -e ".[dev]"`.
For HPC where dependencies already exist, `python -m pip install --no-deps -e .`
registers the checkout in the active environment; do not rely on a global
`PYTHONPATH` pointing to a named parent directory. Run `mdscope --help`,
`python -m mdscope.cli validate <config.yaml>`, or the existing
SLURM scripts after installing.

Build with `python -m build` and validate using `python -m twine check dist/*`.
Check both wheel and source distributions in clean virtual environments from
outside the repository. All seven functional subpackages are included through
setuptools package discovery, not a manually curated list.

The package relocation did not change scientific kernels, YAML schemas, output
paths, units, or trajectory passes. The supported imports are `mdscope.*`;
`mdscope.analysis.*` is not an alias. HPC scripts use the installed lowercase
package consistently.
