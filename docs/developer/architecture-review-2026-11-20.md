# MDscope architecture design review — 2026-11-20

**Scope:** read-only inspection of current `main` source and developer documents; no
trajectories, CI benchmarks, generated dependency graph, or browser UI was run.
The user-installed **Improve Codebase Architecture** skill was visible in a
screenshot but its full instructions were not exposed to this agent, so this
is an evidence-led architecture review, **not a claim to have run that skill**.
Risk levels are proposed priorities, not verified production failures.

## Overall assessment

The existing package boundaries are well chosen. Preserve the separation of
`core` (geometry and persistence), `observables` / `interactions` /
`mechanisms` (scientific measurements), `pipeline` (scheduling/recovery),
`statistics` (replica-level inference), `reporting` (saved-data display),
and `studio` (config authoring). See
[`architecture.md`](architecture.md) and
[`pipeline-graphs.md`](pipeline-graphs.md).

Strong architectural decisions to retain:

- Shared source-coordinate streaming and selective fused-observer execution
  avoid re-reading trajectories for every metric. Membrane observables
  inspect source frames before alignment.
- `registry.analysis_output_plan` and `required_artifacts` anchor
  persisted-product validation and the read-only planner.
- Statistics reuse saved replica products; the code distinguishes replicas
  from trajectory frames as units of independent inference.
- Studio's model holds draft and export reconciliation separate from
  Streamlit widgets; `planning.authoring_catalog` reuses backend metadata.
- The repo already has architecture graph generators, real-data check policy,
  saved-output reference tests, and synthetic production-renderer previews.
  Build on these instead of introducing a new framework or a second registry.

## Findings and bounded recommendations

### A1 · P1 · Expose derived execution requirements in planning and Studio

**Observed code:** `pipeline/analysis_selection.py::extensions_enabled_names`
enables the hydration output when either `analyses.hydration.enabled` is true
**or** `analyses.sdc4_mechanisms.hs_sugar_water_cutoff_A` is not null.
`pipeline/runner.py::_validate_replica_products` in turn validates required
configured outputs before statistics. Thus an explicitly disabled general
hydration toggle does **not** imply that hydration artifacts are optional.

**User-visible consequence:** the documented SDC4 finalization failure
required `extensions/hydration/counts.parquet` for all six replicas despite
`hydration.enabled: false`. The behavior may be scientifically necessary
for HS sugar hydration, but the configuration dependency is surprising.

**Recommendation:** expose a single backend-derived **effective requirement**
in the planner/Studio alongside the stated toggle. Explain which enabled
feature requires hydration and exactly which artifacts must exist. Keep
producer and registry requirements authoritative; do not weaken validation
or silently rewrite `hydration.enabled`.

**Acceptance evidence:** a fixture with general hydration disabled and the
HS sugar cutoff enabled produces the same expected hydration requirement in
planner, UI metadata, and finalization preflight. A configuration with both
disabled does not require it. Verify no extra trajectory passes are scheduled.

### A2 · P1/P2 · Reconcile figure inventory with PDF report outputs

**Observed code:** `reporting/table_reports.py` writes
`statistics_tables.pdf` in each relevant figure directory. In contrast,
`reporting/figure_inventory.py::_figure_paths` recursively discovers only
`*.png`; `figure_snapshot`, `record_generated_figures` and
`figure_inventory` therefore cannot list or record those PDFs.

**Consequence:** the figure/inventory/reporting surfaces can diverge:
table PDFs exist on disk but are not included in the same output audit as
PNG figures. They also lack a dedicated PDF integrity/provenance branch
in `figure_status`, which currently uses Pillow PNG verification.

**Recommendation:** extend the artifact model to support distinct
`image/png` and `application/pdf` render products. Use one authoritative
figure/report inventory and explicit source dependencies. Do not route PDFs
through Pillow; validate that a PDF is readable and preserve provenance.
Keep PNG-only consumers supported until migrated.

**Acceptance evidence:** plot a saved-data fixture with
`plots.table_reports_pdf: true`; inventory lists each PDF once with the
right family and provenance. Re-running plotting updates the PDF record and
does not mistakenly advertise archived table PNGs. No trajectory reads.

### A3 · P2 · Remove numerical CSV side effects from presentation

**Observed code:** `reporting/sdc4_condition_figures.py::_region_figures`
calls `_write_region_replica_summaries`, which writes
`extensions/sdc4_mechanisms/lipid_region_replicate_summary.csv` beneath
each replica. By contrast, `docs/studio.md` and
`docs/developer/pipeline-graphs.md` describe saved-data figure rendering
as presentation rather than a scientific producer.

**Risk:** a figure-only regeneration can change replica outputs, including
file timestamps and possibly downstream provenance or user expectations.
The CSV may be useful, but its authoring stage is unclear.

**Recommendation:** produce the summary in the replica finalization/derived
saved-data stage, or explicitly designate it as a reporting export under
`project_summary/`. Preserve its column contract and existing consumers.
Keep plotting idempotent with respect to per-replica numerical artifacts.

**Acceptance evidence:** snapshot replica extension files, rerun only
plotting, and confirm no numerical files changed; confirm the summary
still exists at the documented destination and matches the old values.

### A4 · P2 · Harden report publication for concurrent jobs and absent pages

**Observed code:** `reporting/table_reports.py` uses the fixed temporary
name `.statistics_tables.pending.pdf` inside each folder, then
`os.replace` and archives old table PNGs on success.

**Risk to investigate:** two plotting processes targeting the same project
can collide on the pending filename. A run that generates no pages in a
folder does not replace an older report for that folder. Whether this
occurs in the supported execution model needs verification.

**Recommendation:** use a unique temporary path plus a per-output-directory
writer ownership rule; define the expected behavior for a now-empty or
unselected report (archive/unregister the previous report versus retain it
with a stale status). Avoid automatic deletion without an archive path.

**Acceptance evidence:** interrupted and concurrent mock writers never
publish partial PDFs or overwrite each other's temporary files. Switching
from a populated report to a configuration with zero selected table pages
does not advertise the old PDF as current.

### A5 · P2 · Reduce orchestration/UI change surface without splitting by size

**Observed code:** `config.py` (~1,050 lines) owns defaults,
parsing/validation, labels, and a few persistence reexports;
`pipeline/runner.py` (~900 lines) owns replica orchestration,
output validation and offline finalization; `planning.py` (~880 lines)
covers output planning, saved-statistics inspection and figure catalogs;
`studio/app.py` (~955 lines) owns numerous Streamlit UI sections.

**Interpretation:** size alone is not an architecture defect. The strongest
specific opportunities are to keep UI presentation functions from acquiring
backend semantics; to provide small typed, read-only result models at the
planner/Studio boundary; and to isolate finalize/preflight **policies**
from runner control flow when a concrete change requires it.

**Recommendation:** use `scripts/architecture_map.py` to inspect import
edges/SCCs and use actual changed-file coupling over several PRs to choose
a candidate seam. Extract one seam at a time with parity tests; never build a
second config parser, analysis registry, or artifact path map.

**Acceptance evidence:** one existing CLI/Studio plan fixture yields
identical IDs, eligibility and output contracts after extraction; actual
producer data and restart behavior remain unchanged.

### A6 · P3 · Refresh release-status documentation

**Observed documentation:** `docs/developer/studio-plan.md` still lists
PR #200 as a remaining review/merge step, though PR #200 is already merged.
Its notes are useful as historical acceptance context but should not read
as current release status.

**Recommendation:** mark the section explicitly historical, date its
snapshot, and point to the merged PR and live issue tracker for current
status. This is documentation hygiene, not an executable architecture change.

## Suggested implementation order

1. **A1** (planner + CLI/Studio effective prerequisites): directly addresses
   costly and confusing ARC validation/recomputation decisions.
2. **A2 + A4** (report identity/provenance and atomic lifecycle): one coherent
   reporting-artifact milestone, with no changes to statistical estimators.
3. **A3** (reporting purity): safeguard numeric outputs against plot-only runs.
4. **A5** (architecture mapping and one evidence-backed seam): only if real
   dependencies or change frequency justify extraction.
5. **A6** (doc refresh): small, independent maintenance PR.

Do not undertake a package-wide rewrite. Preserve user-facing YAML, saved
Parquet/NPZ schemas, replica identities, scientific numerical definitions
and existing result directories.

## Verification before changing runtime code

Generate fresh dependency evidence:

```bash
python scripts/architecture_map.py --output /tmp/mdscope-architecture
python scripts/pipeline_graphs.py --output /tmp/mdscope-graphs
```

Then run the affected targeted tests, automatic Python 3.11/3.13 checks,
and applicable manual figure-review or real-data parity workflows as defined
in [`check-policy.md`](check-policy.md). Neither these scripts nor the
related full skill workflow were run for this design review.
