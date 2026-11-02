# Coordinated sprint checkpoint — 2026-11-02

Base: `main` at `6c3927b` (merged scientific-audit PRs #117/#118).
All 27 open issues and their descriptions were reviewed; none had discussion
comments. Recent merged PRs #90/#94/#96 and #110–#118 were checked to avoid
reimplementing completed grouping, chemistry, motif and residence work.
No issue is closed and no sprint PR is merged; owner review is required.

## Grouped work

1. **Statistics and validation** — `sprint/statistics-validation-controls`,
   PR #119. Independent family enable/draw/CI/test controls, pre-read kind
   filtering, compact method manifest, scalar missing-time safeguards,
   prespecified scalar BH families, exact nonconstant block/permutation
   references, real periodic SASA reference, packaging-test consolidation and
   opt-in timing runs. Default regular-series RNG behavior is preserved;
   invalid temporal intervals are intentionally unavailable rather than
   computed after compressing missing frames.
2. **Molecular observables** — `sprint/molecular-contact-observables`, PR #120.
   Bond/anchor-validated batched lipid matrices, finite-size normalization,
   standard statistics integration and replica heatmaps; offline same-ion
   bridge episodes with censoring; explicit chemical atom subsets after
   sequence mapping; HS shape/attachment coordinates in the existing pass.
   The real test selects `same fragment as (validated anchor selections)`;
   broad residue-name queries can also select protein-attached HS linkers and
   correctly fail the strict lipid population check.
3. **Report integrity** — `sprint/report-integrity-controls`, PR #121. Independent
   plot/PDF controls, pre-read kind exclusions, human-readable labels/units,
   q-based interpretation labels, figure config+code signatures, image hashes,
   source freshness, complete rejection inventory, sampling/source provenance,
   no separator-only pages, and offline adversarial report tests. No new
   analysis/plot registry or second metadata framework was introduced.

Each branch is based independently on the same main checkpoint. The final
combined source tree passed **299 tests**, with two optional timing benchmarks
deselected (see final PR validation). The paired-invalid-replica contrast and render time/selection hash regressions
are included in the final check. Configuration/doc
EOF additions require retaining both `statistics_family_config` and
`presentation_enabled`; remaining code hunks merged automatically. This is a
local integration check, not a merge into main. Preserve every branch's numerical
changes and check CI at the exact reviewed heads before owner merges.

## Issue acceptance accounting

| Issue | Sprint disposition and remaining acceptance |
|---|---|
| #75 | Implementation complete: overlapping import/layout suites consolidated; public config exports intentionally retained because callers use them; atomic rollback/hash and wheel/sdist/outside-checkout gates retained. Owner review/CI determines closure. |
| #76 | Implementation complete: timing tests opt-in, numerical references retained; nondefault cadence is checked for monotonicity, bounded output and completion. Exact 1% assertions cover the documented default cadence. Owner review determines closure. |
| #92 | Partial: independent families, draw/CI/test overrides, method manifest and stale-output guards implemented. A separate point-estimates-only mode and per-family hierarchy controls are not claimed. |
| #89 | Partial: existing #90/#94 grouped lookup retained, early include/exclude reads and progress improved. Strict six-replica production baseline/RNG benchmark remains required before closure. |
| #69 | Partial: real periodic sphere/occluder and target-only residue references supplement existing sphere/cap and burial contracts. Additional heterogeneous lipid geometry/burial sensitivity remains. |
| #70 | Partial: exact nonconstant block distribution, 3-vs-3 p-resolution and invalid-time regression added. Broader nonconstant hierarchy/coverage study remains. |
| #74 | Partial: invalid family/draw/CI/test controls fail early. The full selection, multi-run, pairing and missing-file negative-case matrix remains. |
| #99 | Partial: scalar time validity and report missing/stale/disabled states strengthened. Full six-replica semantic audit remains. |
| #105 | Partial: scalar prespecified BH families and q-based display, finite-replica oracle. Matched windows, convergence and prespecified endpoint decisions remain. |
| #91 | Partial: new lipid matrices, normalization, batching, uncertainty via standard scalar pathway, heatmaps, replay and real-topology validation. Production cutoff/null/area sensitivity and final publication contrasts remain. |
| #103 | Partial: chemical motif site subset and same-ion continuous/intermittent bridge episodes from saved data. Validated K+/HS chemistry, eligible never-observed pair denominators, spatial nulls, direct-vs-solvent-separated interpretation and replica evidence remain. |
| #104 | Partial: HS Rg/anisotropy/vertical shape and attachment coordinates plus linear-chain oracle. Validated microstates, state-conditioned survival and matched windows remain. |
| #71 | Partial: HS analytic geometry supplements existing TM/PBC/area/enrichment oracles. A complete major-observable fixture manifest and all failure cases remain. |
| #72 | Partial: molecular PBC/atom-multiplicity/batch invariance; unsupported planar-box geometry fails explicitly. Wider tilted/triclinic cross-observer policies remain. |
| #73 | Partial: new lipid matrix clean/chunk/replay equivalence. Broader mixed-observer fault injection and finalized-statistic equivalence remain. |
| #108 | Partial: finite-population matrix reference, unavailable-vs-zero semantics, explicit box-derived APL labels. Six-replica composition inventory, scarce-species uncertainty and excluded-area sensitivity remain. |
| #93 | Substantial implementation: independent family/section controls across existing entry points, pre-read identity exclusion, strict stale PNG rejection and non-destructive preservation. Final combined regression passed 299 tests; owner acceptance of the full control matrix determines closure. |
| #106 | Partial: readable labels/units, q-based evidence labels, source capsules and nine-page synthetic PDF review. Full production figure-family readability/summary-dashboard review remains. |
| #107 | Partial: strict image/source/config/code verification, timestamps and actual replica windows. Family-level source metadata is not exact per-datapoint provenance; original legacy analysis hashes remain explicitly unknown. |
| #109 | Partial: expected/present/disabled/missing/invalid figure inventory and adversarial clean offline export. Six-replica report rerender and human production-page review remain. |
| #34 | Partial: existing inventory layer extended and SDC4/shared-state plot responsibilities separated. Renderer/bundle and reusable CLI-validation separation remain where justified. |
| #36 | Partial: new independent numerical, molecular, restart and report references; no execution-only validity claim. Complete redistributable scientific fixture/release checklist remains. |
| #100 | Shared support improved (inference/report controls); no claim to finish motif-specific remaining science or production comparisons. |
| #101 | Shared temporal/interpretation/report safeguards improved; state robustness, matched windows and biological evidence remain. |
| #102 | Shared lipid comparison/statistics/reporting support improved; no new affinity interpretation; production chemical/cutoff sensitivity remains. |
| #98 | No new mapping change: exact mapping and real references already on main. Six-replica recomputation remains blocked on ARC data access. |
| #77 | Partial: parent-process statement/branch gap mapping and six hand-reviewed scientific mutant detections are in `validation-gaps.md`; direct interaction parity and a missing burial-sign reference were added. Historical suites remain until unique behaviors are mapped; this is not a global mutation score. |

## Validation and limits

- PR #119 local: 285 passed, two optional benchmark tests deselected. Independent
  regular-data references, invalid-time behavior and family-isolation tests pass.
  A new paired contrast regression failed before the fix: one invalid replica
  pair could otherwise be silently omitted by a NaN-aware reducer.
- PR #120 local: 269 passed after adding a no-read disabled-matrix-plot regression. Benchmark: 1,024 molecules / 12,288 atoms / three
  frames; batch 256 = 0.750 s / 170.5 MiB RSS; batch 4096 = 0.619 s / 174.7 MiB,
  identical matrix digest. These are synthetic batch tradeoffs, not whole-run
  production speedups. Reproduce with `scripts/benchmark_lipid_contacts.py`.
- Report batch local: 271 tests passed after inline-plot exclusion validation;
  focused report/integrity references pass. Full CI must be checked at the
  final branch head. Synthetic scalar/solvation/coverage pages were rendered
  with Poppler and visually inspected; this is not the historical production PDF.
- PR #120 final-head Python, real-GROMACS smoke and release-solvation benchmark
  CI passed. PR #119 previously passed Python and real-trajectory checks; the final
  paired-contrast correction and report time/selection hash correction require
  reruns at their new heads. CI conclusions
  belong to exact branch heads; check current Actions runs.
  A real-topology test failure exposed broad-selection contamination by HS
  linkers; the selection was corrected, not the population assertion weakened.
- Original six production XTC/TPR paths are unavailable here. No scientific
  issue requiring their convergence, contrasts or report rerender is closed.

## Next work in scientific order

1. Review current CI logs, then the integrated tree; retain strict selection,
   zero/NaN, source-freshness and RNG contracts. Do not merge automatically.
2. On ARC, use the fresh-output commands in `scientific-audit.md`. Compare
   corrected/raw measures within replicas before interpreting condition effects.
3. Finish eligible-pair ion-bridge occupancy/nulls, validated HS microstates,
   and matched-window/threshold/stride sensitivity for #100–#105.
4. Address long-event finalization RAM and physically measured bulk-water
   references. Current streaming updates do not imply bounded finalization RAM.
5. Finish production report/page review and finer source lineage; never bless
   unchanged legacy PNGs or replace unavailable values with zero.
6. Complete broad failure-mode, mutation/coverage and mixed-observer validation.
