# ADR-0001: Target package architecture

**Status:** accepted (2026-11-22) · **Tracking:** [#269](https://github.com/MarionQL/MDscope/issues/269)

## Context

The package is grouped by *kind* of code (`observables/`, `interactions/`, `mechanisms/`), with 15 modules at the package root. An import graph of `src/` (graphify code graph, 2026-11-22) showed:

- **Lower layers importing higher ones.** At least 12 module-level edges went the wrong way: `core` imported `mechanisms`/`observables`, and observables imported `pipeline`. #256 removed the subpackage cases and added a layering test. Remaining: `config` ↔ `core`, `config` → `reporting`, `registry` → `mechanisms`/`pipeline`/`reporting`, `statistics` → `mechanisms`.
- **An analysis's output contract is spread over about 8 places:** `DEFAULTS`, `config_schema`, `FUSED_OUTPUTS`, `required_artifacts` companions, `artifact_schema` branches, `analysis_selection`, planner figure candidates and statistics families.
- **Lightweight contracts live at the package root** (`statistics_*.py`, `permutation_resolution.py`, `state_feature_spec.py`). This is only because `statistics/__init__.py` imports every estimator up front.
- **SDC4-specific code is mixed with generic analyses in `mechanisms/`.** #169 showed the alternative: generic glycan discovery plus an SDC4 preset that reproduces the previous results exactly.
- **The graph's communities form along scientific domains** (solvation postprocessing, lipid contact matrices, HS topology), not along the current folders.

## Decision

Organize by dependency layer, and within the analysis layer by scientific domain:

```
src/mdscope/
├── contracts/     outputs.py (per-analysis descriptors) · statistics.py · features.py
├── config/        defaults.py · validate/ (one module per analysis) · runs.py
├── core/          geometry/PBC, selections, sequence, residence, persistence, checkpoints, progress
├── analyses/      structure/ · membrane/ · interactions/ · solvation/ · states/
│                  (each: streaming observer, offline postprocess, config validator)
├── systems/       sdc4/  (mechanisms, K⁺ audit, HS presets, SDC4 figures and exports)
├── engine/        streaming pass, fused execution, runner, preflight, restart/recovery
├── statistics/    estimators and dispatch only; __init__ imports nothing heavy
├── reporting/     figures, PDF tables, findings
└── interfaces/    cli.py · planning (read-only inspection) · studio/
```

1. **Dependency rule.** A package imports only from packages listed above it. `contracts/` imports nothing from `mdscope`, so planning, Studio and the CLI stay import-light. `tests/test_import_layering.py` enforces the rule; it gains each layer as that layer is completed, and its exception list only shrinks.
2. **One descriptor per analysis.** Each analysis declares its primary output, companions, column/shape schema, named-run support, statistics families and figure IDs in `contracts/outputs.py`. The engine, registry, planner, Studio and statistics read that descriptor instead of keeping their own lists. Adding an analysis becomes one `analyses/<domain>/` module plus one descriptor.
3. **Systems sit on top of generic analyses.** System-specific selections, presets and figures live in `systems/<name>/`. They configure generic analyses, and add mechanisms only where the science is genuinely system-specific. A new system does not touch `engine/`.
4. **`engine/` is separate from `interfaces/`.** The streaming engine is the restart- and performance-critical code. CLI, planning and Studio are thin callers, so a Studio change cannot quietly load trajectory machinery.

## Consequences

- Adding or changing an analysis is local: one domain folder plus one descriptor. Generalizing beyond SDC4 adds a `systems/` entry rather than editing shared code.
- Every step changes internal import paths. Internal aliases are not kept, so each move updates callers, tests, docs and `CLAUDE.md` together, and conflicts with open branches touching the moved modules.
- Unchanged public contracts: CLI commands, YAML keys and defaults, persisted output paths and schemas, and config hashes.

## Migration

Incremental, one issue per PR, in the order listed in [#269](https://github.com/MarionQL/MDscope/issues/269):

1. #256: move shared helpers down a layer; add the layering test.
2. #260: per-analysis output descriptors in `contracts/outputs.py`.
3. #257: statistics contracts into `contracts/`; `statistics/__init__.py` becomes import-light.
4. #259: `config/` package; remove the `config` ↔ `core` cycle.
5. #258: `systems/sdc4/`.
6. Domain regrouping into `analyses/`, `engine/` and `interfaces/`. This step is purely structural, so schedule it when few PRs are open.

**Every step must show:**
- `mdscope plan --format json` is identical before and after for the maintained `synd4_*.yaml` configs and the test fixtures.
- Resolved configs and config hashes are unchanged.
- The full suite and CI pass.
- The real-data `smoke` check passes: the stride-5 copy is enough for pure moves; use the full files when an execution path changes.

## Alternatives considered

- **Keep grouping by kind of code.** Folder names are easy to guess, but one scientific feature is spread over four or five packages, and the root keeps accumulating lightweight modules.
- **Move everything in one PR.** It is faster on paper, but it conflicts with every open branch, cannot be reviewed meaningfully, and mixes moves with the descriptor redesign.
- **Plugin entry points for analyses.** Not needed while all analyses ship in this package. The descriptor registry can be exposed as entry points later without moving code again.
