# Contributing

## Branches and issues

- Use GitHub Flow: branch from current `main`, open a focused PR, review and squash-merge after relevant checks. Delete merged branches.
- Example branches: `feat/statistics-metric-selection-179`, `fix/resume-cache-123`. Resolve conflicts by checking both behaviors and rerunning affected tests.
- Issue title: `[P0-P3][Area] Imperative, testable outcome`.
- Priority: P0 = correctness/data-loss blocker; P1 = important defect/required capability; P2 = routine improvement; P3 = polish.
- One issue per observable outcome: goal, current behavior, acceptance criteria, scope/non-goals, dependencies and reproducer where applicable.
- Reuse established labels. Issues/PRs own status; use Projects only with an actively maintained board and milestones for actual releases.
- Use `Closes #N` only when every acceptance criterion is met; otherwise use `Refs #N` and name the remaining work.

## Cohesion and scope

- Complete the authorized milestone; avoid unrelated cleanup and speculative features.
- Keep related functions together. Split for different ownership, dependencies, invariants or change lifecycles, not line counts or lint scores.
- Before a substantial refactor, consider consolidation. Assess import depth, indirection, cycles, navigation and modules touched per feature.
- Add an abstraction or dependency only for demonstrated correctness or maintenance value. Inline trivial single-use wrappers; retain helpers that isolate complex logic or external I/O.
- Remove superseded internal paths after migrating callers. Keep adapters only for demonstrated public compatibility or a tracked staged migration.
- Reuse algorithms when semantics agree; preserve distinct estimators and independent test oracles.
- Validate at boundaries. Catch expected failures with actionable context; avoid swallowed errors, silent scientific fallbacks and redundant guards.
- Keep scientific choices in the [validated YAML](../configuration.md#yaml-contract-for-generalization), with SDC4 compatibility and non-SDC4 verification.
- Document each concept once. Keep scientific caveats and non-obvious rationale; remove narration, trivial docstring boilerplate and stale guides.

Package ownership and adding analyses: [Architecture](architecture.md).

## Scientific review

- Define behavior, inputs/outputs, units, populations, provenance and failure outcomes.
- Check assumptions, mathematics, state/data flow, integration, verification gaps and recovery.
- Preserve coordinate frames, PBC, mass weighting, selection chemistry, time alignment, replica independence and missing-versus-zero meaning.
- Reproduce substantive bugs and verify the affected entry point; prefer independent analytical references. Never weaken assertions to pass CI.
- Check artifact identity, integrity, completeness and provenance; file existence alone is insufficient.
- Preserve fused/sparse processing; measure trajectory passes, memory, I/O and dense pair allocations when changing hot paths.
- Review behavior and scientific evidence regardless of code origin. Use targeted fault injection where practical; green CI does not establish convergence.

Details: [Scientific contracts](scientific-contracts.md) · [Scientific audit](scientific-audit.md) · [Check policy](check-policy.md).

## Figures

- Same observable and units: shared y limits after drawing series and uncertainty bands. Missing text-only panels inherit limits without determining them.
- Center signed contrast/autocorrelation axes symmetrically where applicable. Different observables/units keep separate scales.
- Same-quantity heatmaps share color limits and a colorbar. Use `reporting.shared_scales.synchronize_y_axes` for comparable panels.
- Dense categorical figures: aggregate each replica independently, keep region order, use one external legend, short display names and full source-table identities.
- Put repeated categorical labels on the final row and use a shared y title.
- Keep zeros in numerical outputs. Record omitted zero-only series; if all region series are zero, omit the figure and archive stale PNGs.
- Test high-cardinality cases through the production renderer.

### Visual review

```bash
python scripts/plot_preview.py --output /tmp/mdscope-preview
```

- Inspect PNGs in `/tmp/mdscope-preview/project_summary/figures` for labels, legends, axes, units and clipping.
- Previews use production renderers and deterministic saved fixtures; no trajectories or inference. Each image is labeled **SYNTHETIC EXAMPLE — NOT MD RESULTS**.
- Extend these fixtures for new renderers; do not build substitute plotting functions.
- Dispatch `MDscope figure review` for plotting PRs. The default comparison is the merge base with `main`; `base_ref` selects another revision.

| Path in `mdscope-figure-review` artifact | Purpose |
|---|---|
| `index.html` | Updated-code gallery organized by biological condition |
| `review/index.html` | Code-revision comparison |
| `review/selection/index.html` | One-plot comparison |
| `current/table_previews/table_inventory.csv` | Synthetic Parquet inputs; adjacent `*.columns.csv` and `*.sample.csv` give schema and ≤5 rows |
| `review/previous`, `current`, `review/previous-selected`, `review/updated-selected` | Code-review inputs, not biological conditions |

Visual review does not certify scientific validity.
