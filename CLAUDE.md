# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Read [AGENTS.md](AGENTS.md) first.** It is the authoritative checklist (scope, merge criteria, which guide owns which topic). This file only adds commands and the big-picture map; do not duplicate the docs it links.

## Commands

```bash
python -m pip install -e '.[test]'        # add ,sasa for FreeSASA (no Windows wheel: use Linux/WSL)
ruff check src tests scripts              # CI-blocking: E9 + full pyflakes (F)
pytest -q tests                           # default excludes -m benchmark
pytest -q tests/test_solvation.py::test_hs_chain_grouping_reuses_covalent_topology_not_chain_ids
pytest -q -s -m benchmark tests/test_performance_benchmarks.py
git ls-files -z '*.sh' | xargs -0 shellcheck   # CI-blocking on the 3.11 job
python scripts/architecture_map.py --output /tmp/mdscope-architecture   # AST module/call graph
```

- CLI tests (`tests/mdscope_test_support.py::invoke_cli`) spawn the `mdscope` executable from `PATH`. When testing a worktree without installing it, put the venv's scripts dir on `PATH` and set `PYTHONPATH=<worktree>/src` (inherited by the subprocess) so the right source is exercised.
- Intentional re-exports must use `import name as name` (Ruff F401 is blocking).
- CI (`.github/workflows/mdscope-checks.yml`) runs Python 3.11/3.12/3.13 plus a 3.11 minimum-dependency job. Figure review (`figure-review.yml`) is the only `workflow_dispatch` job; dispatch it for renderer/figure changes.

### Real GROMACS data (local only)

Trajectories are never committed or uploaded. The four `test-data-v1` TPR/XTC files live in git-ignored `examples/data/`, verified against `tests/data/real_test_data.sha256`.

```bash
scripts/run_real_data_checks.sh smoke     # also: hbond solvation performance finalization box-sweep box-clearance all
python scripts/stride_real_data.py --stride 5                 # fast copy in examples/data/stride5
MDSCOPE_REAL_DATA=examples/data/stride5 scripts/run_real_data_checks.sh smoke   # iteration only, not merge evidence
```

`MDSCOPE_REAL_DATA` / `MDSCOPE_REAL_OUTPUT` override the data and output dirs (output defaults to git-ignored `real_test_results/`). Which checks a change needs: [check policy](docs/developer/check-policy.md).

## Architecture

`mdscope` CLI (`cli.py`) → YAML resolved and validated by `config/` / `config_schema.py` → topology preflight (`pipeline/preflight.py`) → `pipeline/runner.py` per replica → **one shared streaming pass** over the trajectory (`pipeline/streaming.py`) feeding fused observers (`pipeline/fused_membrane.py`, `pipeline/fused_extensions.py`, SDC4 mechanisms, core observables) → chunked, checkpointed per-replica products → replicate summary + artifact validation → offline `statistics/` → findings and optional figures (`reporting/`).

Contracts that span several files:

- **Observers, not passes.** New per-frame measurements attach as observers to the shared pass. Adding a separate trajectory read is a design smell; the few exceptions (multi-origin MSD, potassium-bridge audit, DSSP variants) are scheduled in `pipeline/replica_tasks.py`.
- **Coordinate frames.** Membrane observers see source coordinates *before* protein fitting; never rotate/recenter shared frames. PBC/unit semantics: [scientific contracts](docs/developer/scientific-contracts.md).
- **Persisted-output contract.** `contracts/outputs.py` (canonical per-analysis paths, static companions, schemas and producer constants), `output_paths.py` (existing path exports), `registry.py` (`analysis_output_plan`, configuration-dependent `required_artifacts`) and `planning.py` (introspection) are the single source for what a run must produce. Restart/resume and `statistics`/`finalize` read saved products only and rely on these; a new output file must be registered there, and changing a required companion affects resumability of existing outputs.
- **Writes are atomic** via `core/persistence.py` / `core/checkpoints.py`; scientific kernels should not write files directly.
- **Statistics are replica-aware.** Replicas, not frames, are the independent units (`statistics/`).
- **YAML generalization.** User-selectable scientific choices go in the shared validated YAML (defaults in `config.py`, schema in `config_schema.py`, surfaced by `planning.py` and Studio), while the maintained `synd4_*.yaml` SDC4 configurations must load and reproduce unchanged. See the [YAML contract](docs/configuration.md#yaml-contract-for-generalization).
- **Package layering.** `core` → `observables`/`interactions` → `mechanisms` → `statistics`/`reporting` → `pipeline` → `cli`/`studio`; lower packages never import higher ones (`tests/test_import_layering.py`). Shared residue-identity and progress helpers live in `core/sequence.py` and `core/progress.py`.
- **No internal aliases.** Old internal module paths are not kept; CLI/YAML/output compatibility is the public contract.
- Similar-looking helpers are sometimes deliberately distinct (e.g. `biophysical_common.heavy` vs `mechanistic_common.heavy`); check [architecture](docs/developer/architecture.md) before merging "duplicates".

## Repo conventions

- GitHub Flow, squash merges, issue-titled branches (`feat/...-NNN`, `fix/...-NNN`); issue/PR conventions in [contributing](docs/developer/contributing.md).
- PRs record executed evidence: commit, commands, outcome, skipped checks and limits.
- `*.sh` and `*.sha256` must stay LF (`.gitattributes`; Slurm rejects CRLF). On Windows, write such files with `newline="\n"`.
