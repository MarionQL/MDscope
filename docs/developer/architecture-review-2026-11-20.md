# MDscope architecture review — 2026-11-20

**Source snapshot:** `935a16c30679b9e93eb1983c46e9798b6bbec0b4`, whose production source matches `main` at `dc0615311b875d9c7361182cdd18fbb5c6f1589a`. This replaces the earlier review with fresh source, dependency, history and targeted execution evidence. Scope: identify bounded improvements; no runtime refactor or scientific calculation.

## Assessment

Keep the existing package ownership. Shared source-coordinate streaming, selective observer recovery, registry-owned required products and saved-data replica statistics provide useful depth. Studio delegates numerical execution and topology validation, and its draft/export model is separated from Streamlit widgets. Preserve distinct scientific estimators.

The generated AST inventory contains **105 modules, 317 import edges, 953 statically resolved call edges, and zero multi-module import SCCs**. This argues against a cycle-driven package rewrite, but does not establish correct dynamic calls or lifecycle behavior. Execution/recovery graphs are conceptual maps; import graphs derive from current source.

In the latest 40 source-touching commits, the most frequently changed files are `config.py` (15), `pipeline/runner.py` (14), `reporting/plotting.py` (9), and `registry.py`, `statistics/profiles.py`, `statistics/time_series.py` (8 each). Counts identify exploration targets; they do not measure change coupling or justify splitting files by length. No glossary or ADR files were found in this checkout.

## Ranked candidates

| Order | Improvement | Strength | Effort / migration risk |
|---|---|---|---|
| 1 | Complete report publication, identity and provenance | Strong | Medium / medium |
| 2 | Give statistics-manifest interpretation one owner | Strong | Medium / low–medium |
| 3 | Unify saved-result validation facts and readiness policies | Worth exploring | Medium–high / high |
| 4 | Move SDC4 numerical summary export out of rendering | Strong | Small–medium / medium |

Effort is relative scope, not a calendar estimate. Rank reflects verified gaps, caller leverage and migration cost. These are candidate milestones, not authorization to implement all four.

### 1. Complete report publication, identity and provenance

**Evidence:** `reporting/table_reports.py:26–46` groups pages and publishes `statistics_tables.pdf`; production plotting uses that session at `reporting/plotting.py:989`. However, `reporting/figure_inventory.py:75–86` discovers/snapshots only PNGs, `:133` records provenance from that discovery, and `:179` verifies only PNG content. A direct probe generated a genuine PDF through the production report session: the PDF existed with a `%PDF-` header, while `_figure_paths` returned `[]`.

An empty successful session left a previous PDF unchanged. This proves retention, not that every retained PDF is stale: desired behavior for zero selected pages needs an explicit product policy. The fixed pending filename at `table_reports.py:31` is a conditional concurrency risk; support for concurrent same-project plotting has not been established.

**Deepen:** extend the existing report publication module to own both PNG/PDF product identity, publication and source/page membership; inventory consumes those facts. Give report sessions unique temporary files and resolve stale-report handling. Avoid a generic artifact framework.

**Before → after:** publication writes PDFs while inventory separately scans PNGs → one reporting product owner supplies coherent inventory and provenance.

**Deletion test:** callers no longer reconstruct PDF contents from filenames or know pending paths and archival rules. Locality improves when selection/publication/provenance changes land together; leverage improves when CLI and Studio consume the same published-product facts.

**Acceptance:** through production plotting, inventory each PDF once with source/page membership, validate PDF separately from PNG, preserve old outputs on failed rendering, and verify populated-to-empty reruns against the agreed archival policy. Extend `tests/test_table_pdf_reports.py` and inventory-consumer coverage. Test interruption/concurrency against the writer model actually supported. No trajectory reads or estimator changes.

### 2. Give statistics-manifest interpretation one owner

**Evidence:** lifecycle writing is at `statistics/__init__.py:49–75,112–118`; independent readers are at `reporting/figure_inventory.py:194`, `planning.py:536,707`, `pipeline/runner.py:773`, and `reporting/report.py:157`. The figure reader defaults absent status to complete at `figure_inventory.py:206`; runner `:792–797` and findings generation `:170` require explicit completion.

A direct call to production `statistics_status` returned `valid` for `{"methods": {}}`, `incomplete` for an explicit running manifest, and `invalid_manifest` for a JSON list. This is a verified policy mismatch, not evidence of corrupt numerical results. Some legacy compatibility is deliberate.

**Deepen:** put lightweight manifest parsing, shape checks, completion and family facts under one read-only owner near existing statistics contracts. Consumers choose explicit strict-publication or permitted-legacy policies and retain their own error/presentation context. Preserve output schemas and estimands; leave bootstrap implementations separate.

**Before → after:** runner, planner and renderers separately infer completion → one manifest interpretation with explicit consumer policies.

**Deletion test:** remove duplicate JSON/default/completion decisions from callers. Merely wrapping `json.loads` would be shallow. Locality means malformed/legacy behavior is specified once; leverage comes from every saved-statistics consumer sharing the same facts.

**Acceptance:** one behavior matrix for missing, malformed, nonobject, running, complete, missing-status, disabled-family and legacy manifests, plus planning/plotting/findings entry-point coverage. Retain input-change/incomplete-manifest coverage in `tests/test_report_integrity.py` and failure preservation. No new hash gating or unplanned rejection of compatible historical outputs.

### 3. Unify saved-result validation facts and readiness policies

**Evidence:** `pipeline/runner.py:82` validates registered artifact structure; `:536` adds offline finalization checks; `:618` assembles Studio's saved-result inventory. `studio/model.py:309` delegates inventory, but `:315,345` reconstructs coverage/readiness from inventory rows.

Guarantees differ: inventory `runner.py:679` checks summary readability; offline preflight `:555–563` also checks nonempty, finite, unique, monotonic time, and `:577–586` compares checkpoint identity/count. Thus `verified_existing` and ready for offline finalization are distinct guarantees. Source inspection establishes the difference; an actual production misclassification was not reproduced. Some policy differences may be intentional.

**Deepen:** concentrate persisted-product inspection facts and contextual readiness policies in a pipeline-owned module, reusing `registry.py` and existing scientific validators. Studio presents those facts; finalization enforces the relevant policy. Keep topology checking separate from saved-data validation.

**Before → after:** inventory and finalization separately assemble readiness → shared inspection facts with explicit reuse/finalization policies.

**Deletion test:** eliminate repeated coverage-key construction and validation decisions from consumers. If callers still reconstruct readiness rules, the extraction is shallow. Locality puts diagnostics and validation facts together; leverage lets Studio and CLI reuse them without importing unrelated trajectory orchestration.

**Acceptance:** inventory and actual finalization entry points agree on facts for valid, malformed-time, identity-mismatch, missing-companion and partial-coverage fixtures. Preserve stage-marker/config/schema regressions in `tests/test_studio.py:1320–1361`, historical reuse, no additional trajectory passes and no mutation on failed preflight. Run relevant recovery/entry-point tests and applicable real-data checks before landing runtime changes.

### 4. Move SDC4 numerical summary export out of rendering

**Evidence:** `reporting/sdc4_condition_figures.py:60–90` aggregates and writes `extensions/sdc4_mechanisms/lipid_region_replicate_summary.csv` for each replica; `:121` calls it while drawing region figures. A direct probe of the production helper wrote that replica CSV from two saved-style rows and produced mean contact fraction `0.3`. This was a helper-level probe, not a full plot-only project run.

**Deepen:** give saved-data summary/export ownership of aggregation and persistence. The renderer consumes the result, or emits a clearly designated reporting export under `project_summary/`. Check consumers before choosing the destination; preserve columns and independent replica aggregation. Historical outputs need an explicit compatibility decision.

**Before → after:** region rendering computes and rewrites replica summaries → one summary owner; rendering presents its result.

**Deletion test:** remove numerical CSV publication from rendering rather than adding a pass-through function. Locality puts aggregation and independent expected values together; leverage lets other consumers reuse it without invoking Matplotlib.

**Acceptance:** snapshot replica numerical artifacts, run production plotting, and prove they remain unchanged. Verify destination, schema and values against the established contract; test missing moieties and empty data. Retain selection and shared-scale behavior. No estimator or trajectory changes.

## Corrections to the previous review

- **A1 retired as a missing-feature claim:** `planning.py:255–281` already emits a complete implicit hydration descriptor labeled required by SDC4; `tests/test_analysis_plan.py:197–212` checks it. Keep this behavior.
- **A2 + A4 consolidated** into candidate 1, with PDF omission and empty-session retention directly probed. Concurrency remains conditional.
- **A3 retained** as candidate 4 with direct helper evidence.
- **A5 replaced** by concrete manifest/readiness seams. File length alone is not a recommendation.
- **A6 is documentation hygiene:** date historical delivery notes rather than treating them as an architectural refactor.

## Verification and limitations

Executed against the cited source snapshot:

```bash
python scripts/architecture_map.py --output /tmp/mdscope-architecture-review-20261120
python scripts/pipeline_graphs.py --output /tmp/mdscope-pipeline-review-20261120 --inventory /tmp/mdscope-architecture-review-20261120/architecture.json
git log -40 --format= --name-only -- src/mdscope
```

Both graph scripts succeeded; Graphviz generated execution, recovery and import SVGs. An independent fresh-context source exploration checked seams and existing tests. Direct production-function probes used temporary files only: `table_report_session`/`add_table_page`, `_figure_paths`, `statistics_status`, and `_write_region_replica_summaries`. No production project was modified.

| Probe | Observed |
|---|---|
| Genuine PDF publication | `%PDF-` header, PDF exists |
| Discovery after PDF-only publication | Empty list |
| Empty session after existing PDF | Previous bytes remain unchanged |
| Manifest with methods but no status | Figure reader: `valid`; explicit-completion predicate: false |
| Running manifest | Figure reader: `incomplete` |
| Complete manifest | Figure reader: `valid` |
| JSON-list manifest | Figure reader: `invalid_manifest` |
| SDC4 renderer summary helper | Replica CSV written; mean contact fraction `0.3` |

No pytest suite, GROMACS/trajectory run, benchmark or Studio browser session was executed. This shell lacks pytest, pyarrow and MDAnalysis; probes used installed NumPy/pandas/Matplotlib and lightweight modules. They do not establish full entry-point behavior, convergence or speed. Follow [check policy](check-policy.md) when implementing a selected candidate.

**Recommended first milestone:** candidate 1 closes a reproduced gap in a recently changed reporting path, needs no scientific-kernel changes and can be verified with saved-data rendering fixtures. Candidate 2 is a separate follow-up; defer candidate 3 until its policy differences are resolved.

## Issue tracking

Follow-ups: [#215](https://github.com/MarionQL/MDscope/issues/215) (report publication), [#216](https://github.com/MarionQL/MDscope/issues/216) (manifest interpretation), [#202](https://github.com/MarionQL/MDscope/issues/202) (existing saved-result readiness issue), and [#217](https://github.com/MarionQL/MDscope/issues/217) (rendering/export ownership). See the [backlog review](backlog-review-2026-11-20.md) for current priorities, superseded claims and YAML generalization requirements; [#218](https://github.com/MarionQL/MDscope/issues/218) tracks explicit residue-response YAML control.
