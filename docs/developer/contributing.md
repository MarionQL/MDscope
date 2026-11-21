# Contributing

**Rule:** reuse one authoritative implementation per scientific operation where possible; keep distinct estimators separate even when their code looks similar.

## Branches and pull requests

Use **GitHub Flow**: `main` is the integration branch; each issue or focused change normally uses a short-lived branch from `main` (for example, `feat/statistics-metric-selection-179` or `fix/resume-cache-123`). Open a PR, run relevant checks, review the diff, then squash-merge. Do not commit unfinished features directly to `main`, create long-lived `develop` branches, or merge failing changes. Update a PR branch against `main` when needed; resolve conflicts by verifying both behaviors and rerunning affected tests. Delete merged branches.

## Issue tracking

- Title format: `[P0-P3][Area] Imperative, testable outcome` (example: `[P1][Statistics] Select inference per metric`). Use P0 for blocking correctness/data loss, P1 for important defects or required capabilities, P2 for routine improvements, and P3 for low-priority polish.
- One issue = one observable outcome. Include **goal**, **current behavior**, **acceptance criteria**, **scope/non-goals**, and any dependencies or reproducible examples. Avoid speculative implementation details.
- Use existing `bug`/`enhancement` labels and established area labels; do not create near-duplicate labels. Link PRs with `Closes #N` only when all acceptance criteria are met; otherwise use `Refs #N` and state remaining work.
- Use GitHub Projects for prioritization only when a shared board is actively maintained; issues/PRs remain authoritative. Milestones are for actual releases, not arbitrary work buckets.

## Avoid codebase bloat

- **Dead paths:** remove demonstrably unused code within scope after checking callers and public compatibility; do not expand a bug fix into unrelated cleanup.
- **Error handling:** validate at boundaries and propagate scientific failures; catch expected exceptions when adding actionable context or preserving recoverable state. Avoid silent fallbacks.
- **Duplication:** reuse existing utilities where semantics agree; allow small adapters and independent test oracles when they verify different contracts.
- **Dependencies:** prefer existing capabilities when suitable; add a maintained dependency if it materially improves correctness or maintainability.
- **Generalization:** implement user-selectable scientific choices in the shared validated YAML/CLI, preserve SDC4 configurations and explicit overrides, and verify non-SDC4 plus applicable SDC4 parity. See [YAML contract](../configuration.md#yaml-contract-for-generalization).
- **Scope and docs:** do not add unrequested features, configuration options, wrappers, or lengthy commentary. Update the authoritative short guide and delete obsolete wording.

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
| `cli.py`, `config.py` | CLI and YAML configuration |
| `output_paths.py`, `registry.py` | Shared output path constants, persisted artifact schemas and required companions |
| `pipeline/replica_tasks.py` | Specialized executable task dispatch and saved scalar inputs |

## Comparable subplot scales

Within one figure, panels showing the same y-axis observable and units
must have identical y-axis limits. Derive the common range from measured
panels after drawing all series and uncertainty bands; text-only missing
panels do not influence the limits but inherit them. Center signed
autocorrelation/contrast axes symmetrically around zero when applicable.
Do not synchronize different observables or units. Multi-panel heatmaps
of one quantity must share color limits and colorbar; axis coordinates
remain in their documented units. Use
`reporting.shared_scales.synchronize_y_axes` rather than separate
per-renderer scaling formulas.

## Dense multi-replica categorical figures

For saved-data small multiples (e.g., SDC4 lipid-region specificity),
pre-aggregate each replicate independently, keep the same region ordering
and y-axis scale across comparable panels, and avoid subplot-local
legends. Reserve a figure-level legend outside the data area and display
short region names while retaining full names in the source tables.
Repeated axes should show categorical tick labels on the final row only,
with a shared y-axis title instead of three overlapping copies.

Do not draw every all-zero lipid/leaflet series on top of one another;
keep those zero values in the CSV/Parquet outputs and record how many
zero-only series were omitted from a mixed-value figure. If every
series in a region figure is zero, omit the figure and archive a stale
PNG instead of leaving an uninformative blank plot in active results.
Test representative high-cardinality figures with the production renderer.

## Visual review for plotting changes

```bash
python scripts/plot_preview.py --output /tmp/mdscope-preview
# Open PNGs under /tmp/mdscope-preview/project_summary/figures in VS Code
```

Previews call the **production** figure renderers using deterministic fake saved tables and arrays (no TPR/XTC or statistical inference). They include time series, profiles, condition contrasts, interaction occupancy/lifetimes, RDFs, survival curves, lipid enrichment/contacts, diffusion, and SDC4 maps/states. Every PNG carries a prominent **SYNTHETIC EXAMPLE — NOT MD RESULTS** header outside the plotting area. Plotting PRs run the `MDscope figure review` GitHub Action; its `mdscope-figure-review` artifact opens at **`index.html`**, a condition-first preview of the *updated code* (not side-by-side duplicates). The separate **`review/index.html`** compares code revisions; **`review/selection/index.html`** isolates a one-plot request. `current/table_previews/table_inventory.csv` lists every synthetic Parquet input, with matching `*.columns.csv` for schema and `*.sample.csv` with at most five rows. `review/previous`, `current`, `review/previous-selected` and `review/updated-selected` are code-review inputs, **not biological conditions**. Inspect labels, legends, axes, units and clipping. The gallery is not an automatic scientific-validity gate; when adding a new renderer, extend the saved-data fixture rather than build a substitute plotting function.

## Change checklist

- **Scope:** inspect callers, relevant contracts, tests and docs. Complete the authorized milestone across its affected responsibilities; keep ownership clear and avoid speculative abstractions. Small commits are review units, not delivery limits.
- **Compatibility:** migrate known callers before removing an internal module; retain a tested adapter when external use or staged migration requires it. Preserve public CLI/YAML/output contracts unless explicitly migrated.
- **Scientific invariants:** retain units, populations, coordinate frame, PBC, mass weighting, selection chemistry, time sampling, missing-vs-zero semantics and results.
- **Performance:** minimize trajectory traversals, memory, I/O and dense atom-pair materialization. Preserve fused and sparse processing.
- **Testing:** use independent numerical oracles for mathematical claims; check integration through the actual entry point. A successful run or two agreeing implementations may share the same bug.
- **Failure safety:** test invalid inputs and checkpoint corruption without destroying valid existing outputs. Verify artifact identity, integrity, completeness and provenance—not just existence.
- **Validation:** run relevant targeted tests and applicable Python 3.11–3.13 CI; run real GROMACS checks for affected behavior. Distinguish passed, failed, skipped and unavailable checks.
- **Documentation:** one authoritative location per concept. Prefer short tables and bullets; remove stale guides and avoid mass formatting.

See [Scientific contracts](scientific-contracts.md), [Scientific audit](scientific-audit.md) and [Validation gaps](validation-gaps.md) for the detailed acceptance criteria.

**AI-assisted changes:** review behavior and scientific evidence, not the origin of the code. For substantive fixes, demonstrate the prior failure or a targeted fault-injection/mutation detection where practical. A green CI run is not evidence of physical convergence.
