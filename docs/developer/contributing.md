# Contributing

**Rule:** one implementation per operation; preserve scientifically different contracts even when code looks similar. Prefer net simplification.

## Branches and pull requests

Use **GitHub Flow**: `main` is the integration branch; each issue or focused change uses a short-lived branch from `main` (for example, `feat/statistics-metric-selection-179` or `fix/resume-cache-123`). Open a PR, run relevant checks, review the diff, then squash-merge. Do not commit incomplete work directly to `main`, create long-lived `develop` branches, or merge failing changes. Update a PR branch against `main` when needed; resolve conflicts by verifying both behaviors and rerunning affected tests. Delete merged branches.

## Issue tracking

- Title format: `[P0-P3][Area] Imperative, testable outcome` (example: `[P1][Statistics] Select inference per metric`). Use P0 for blocking correctness/data loss, P1 for important defects or required capabilities, P2 for routine improvements, and P3 for low-priority polish.
- One issue = one observable outcome. Include **goal**, **current behavior**, **acceptance criteria**, **scope/non-goals**, and any dependencies or reproducible examples. Avoid speculative implementation details.
- Use existing `bug`/`enhancement` labels and established area labels; do not create near-duplicate labels. Link PRs with `Closes #N` only when all acceptance criteria are met; otherwise use `Refs #N` and state remaining work.
- Use GitHub Projects for prioritization only when a shared board is actively maintained; issues/PRs remain authoritative. Milestones are for actual releases, not arbitrary work buckets.

## Avoid codebase bloat

- **Dead paths:** remove superseded functions, unused exports, stale branches and their tests when changing behavior.
- **Defensive noise:** validate at boundaries; do not add redundant guards or broad `except` blocks that suppress real errors.
- **Duplication:** reuse existing utilities rather than copying similar implementations; preserve genuinely different scientific estimators.
- **Dependencies:** prefer the standard library and installed packages; justify new packages and avoid overlapping functionality.
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

## Visual review for plotting changes

```bash
python scripts/plot_preview.py --output /tmp/mdscope-preview
# Open PNGs under /tmp/mdscope-preview/project_summary/figures in VS Code
```

Previews call the **production** figure renderers using deterministic fake saved tables and arrays (no TPR/XTC or statistical inference). They include time series, profiles, condition contrasts, interaction occupancy/lifetimes, RDFs, survival curves, lipid enrichment/contacts, diffusion, and SDC4 maps/states. Every PNG carries a prominent **SYNTHETIC EXAMPLE — NOT MD RESULTS** header outside the plotting area. Plotting PRs run the `MDscope figure review` GitHub Action; its `mdscope-figure-review` artifact opens at **`index.html`**, a condition-first preview of the *updated code* (not side-by-side duplicates). The separate **`review/index.html`** compares code revisions; **`review/selection/index.html`** isolates a one-plot request. `current/table_previews/table_inventory.csv` lists every synthetic Parquet input, with matching `*.columns.csv` for schema and `*.sample.csv` with at most five rows. `review/previous`, `current`, `review/previous-selected` and `review/updated-selected` are code-review inputs, **not biological conditions**. Inspect labels, legends, axes, units and clipping. The gallery is not an automatic scientific-validity gate; when adding a new renderer, extend the saved-data fixture rather than build a substitute plotting function.

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
