# ADR-0002: Evidence-first architecture and cohesive modules

**Status:** accepted (2026-11-22); implemented through [#269](https://github.com/MarionQL/MDscope/issues/269) and [#279](https://github.com/MarionQL/MDscope/issues/279). Supersedes ADR-0001's mandatory `analyses/` · `engine/` · `interfaces/` migration (removed; see git history). The current layout is described in [Architecture](../architecture.md).

## Context

- Repeated ownership required many edits to add one runnable, resumable, selectable and plottable analysis.
- Independent scheduling paths reread trajectories; folder moves and forwarding layers did not reduce either cost.
- Reject mandatory domain-folder migration: moving files alone does not reduce edit count or trajectory passes.

## Decision

Keep role-based packages and fix *ownership* instead:

1. **One owner per fact.** Each scheduling decision, output identity, family list or capability set is defined once and imported elsewhere. A test pins every rule that has regressed.
2. **Correct dependency direction.** Lower layers never import higher ones; `contracts/` is dependency-free; system-specific science sits above generic science and below `pipeline/`. `tests/test_import_layering.py` enforces it.
3. **One execution decision per replica.** The runner derives pending work from registered outputs and replays compatible missing work in a single source-frame pass.
4. **Move or split modules only on evidence**: a demonstrated independent reason to change, a forbidden dependency, or real duplication. Size alone is never a reason.

## Hard gates for any architectural change

- Preserve CLI/YAML keys and defaults, resolved config hashes, `mdscope plan --format json`, persisted paths and schemas. Internal import paths are not public and get no permanent adapters.
- Preserve coordinate frames (membrane observers see source coordinates before protein fitting), PBC, units, identities, temporal order, replica independence and missing-versus-zero meaning.
- Preserve restart safety: atomic publication, failure before mutation, no-op reruns that reread nothing, and no additional full trajectory traversal without measured justification.
- Keep planning, validation and the CLI import-light.
- Verify with config/plan parity on the maintained configurations, the Python 3.11–3.13 matrix, a built wheel/sdist, and real-data checks when execution changes ([Check policy](../check-policy.md)).

## Design test before adding a module or abstraction

Before adding a module or abstraction, identify:

- Its owner and independent reason to change.
- The duplication/dependency it removes and navigation cost it saves.
- Any adapters it would require.

Without specific benefits, keep the cohesive owner. All analyses ship together; plugin entry points, per-analysis base classes and extra registries have no demonstrated need.

## Outcome (measured, 2026-11-22)

| Gate | Before | After |
|---|---|---|
| Import cycles (architecture map) | 1 (`figure_inventory` ↔ `table_reports`) | 0 |
| Reverse dependencies | `registry`/`systems.sdc4` → `pipeline`, `core.alignment` → `config`, `contracts` → `output_paths` | none |
| Family and capability lists | statistics/figure families and group-capable analyses copied in up to 9 places; one copy rejected `statistics.families.lipid_contacts` | one owner each in `contracts/` or `config/defaults.py` |
| Replica replay | four family-specific fallback paths | one `run_fused_only` replay; specialized tasks declare `produces` |
| Real-data recovery, three products deleted | 4 trajectory opens, 126 decoded frames, 63.1 s | 1 open, 44 frames, 22.7 s; outputs identical |
| Full-run wall time | — | unchanged within noise |
| YAML / plan / config parity | — | identical for all maintained configurations |

Reviewed and deliberately **not** moved, for lack of a seam with an independent reason to change: `observables.multi_origin_pipeline`, the SDC4 streaming observer, `reporting/plotting.py`, `mechanisms/state_outcomes.py`.

## Consequences

- Adding an analysis still touches several owners (YAML, kernel, scheduling, output contract, metadata) because those change for different reasons, but no step copies a fact owned elsewhere ([Adding an analysis](../architecture.md#adding-an-analysis)).
- A domain folder for one area remains possible later, if a measured change-map shows it would reduce the edit surface.
