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

## Offline reporting completion — 2026-11-03

The reporting batch completes the remaining #34 separation: PDF rendering,
archive export, reusable table inventory/eligibility and lightweight labels
are separate responsibilities; CLI/programmatic structural validation was
completed in #124. Stable CLI commands, PDF/ZIP names, figure paths and prior
bundle member names are retained. Original Parquet members and table/content
manifests are additive. Unchanged numerical estimators and fused passes remain
outside reporting.

Validation: **423 local tests passed**, two optional benchmarks deselected;
wheel/sdist, twine, compileall and adversarial offline export checks passed.
A seven-page synthetic evidence PDF was rendered with system Poppler at
150 dpi and visually inspected: sampling/warnings, dense scalar tables with
long acidic-patch labels and separate units, configured 90% intervals,
replicate tables, methods, atlas figures/captions and sharing instructions.
This is a presentation fixture, not physiological evidence. Original report
pages and the remaining unique scientific plot families still require review.
No binary snapshots or production trajectories were added.

#106/#107/#109 remain **partial**: larger table font, explicit estimate/interval/
p/q lines, human labels and named contrast direction; common PDF/ZIP eligibility;
new statistical input/config freshness lineage; incomplete-run exclusion;
duplicate-image rejection and displayed-content audit. Export failures retain
the previous PDF/ZIP/manifests; disabling bundling archives a previous ZIP.
Legacy statistical inputs remain explicitly unverified rather than claiming
original parameter provenance. Source fingerprints are metadata checks, not
content hashes. The final individual file replacements are not an atomic
multi-file transaction. Production-page QA, exact per-datapoint lineage,
summary dashboard and broader plotting review remain open.

## Consolidated completion follow-up — 2026-11-08

This section supersedes earlier issue-disposition/freshness statements above
where they conflict with merged PRs #124/#125/#132/#133 or the owner's
structure-compatible legacy-output policy. Historical signatures are not
reuse eligibility gates.

At the owner's direction, #91, #99, #100, #101, #104, #105, #107 and #109
were closed with explicit outstanding six-production-replica checks recorded
in each issue. Closure does not establish scientific convergence, endpoint
prespecification for a future study, biological mechanisms, or production
report acceptance.

One follow-up branch consolidates:
- #128: required companion/standalone output checks and preservation on failure;
- #130: optional measured geometric aqueous-region reference, distinct from
  full-cell density and shell accessible-volume correction;
- #103: optional explicit eligible-pair sampled bridge occupancy, including
  never-observed pairs; existing identity/censoring semantics remain;
- #106: compact q-based evidence overview using existing scalar statistics;
- #72/#36: simultaneous source-space geometry/contact observer invariance,
  independent density/excluded-shell references and failure injection.

Remaining code/method work must not be closed merely because it has some
passing tests: #102/#108 need validated protein-excluded-area/composition
sensitivity; #103 still needs chemical/spatial-null interpretation and real
ion/site validation; #130 needs real-trajectory slab/plateau sensitivity and
an accessible-shell strategy if justified; #72 needs broader tilted/triclinic
cross-observer policies; #36 retains release-wide validation/license/citation
work. #131's adaptive endpoint-precision stopping is unimplemented: fixed
per-family draws remain the tested default. It needs a justified resource
benchmark and coherent inferential-family integration rather than extra
configuration without demonstrated benefit. These limitations remain open.

Validation for this follow-up:
- Local Python 3.12 suite with real FreeSASA: 454 passed; two optional
  benchmark tests deselected. Focused completion/invariance/reference checks
  passed again after the final NPZ-member validation change.
- Two new entry-point regressions were run against the original code and
  failed: standalone DSSP omission reached finalization, and a completed run
  reused a distance matrix after its time axis was removed. Both now fail
  explicitly before mutation or trajectory reads. Further regressions cover
  damaged JSON/NPZ, valid-Parquet wrong columns, and mismatched axes.
- A separate subprocess run against an untouched main snapshot and the
  modified tree produced 34 exactly matching finalized Parquet tables on the
  mixed-observer fixture, including statistics. No default estimator change
  is claimed. This is synthetic numerical parity, not production throughput.
- Wheel/sdist builds and twine checks passed; compile/whitespace/shell checks
  passed. A 12-page synthetic report was generated, and its new evidence
  overview/dense scalar page was visually inspected with system Poppler at
  150 dpi. No binary fixtures were committed.
- Python 3.11/3.13 and real-GROMACS checks remain GitHub CI gates. GROMACS and
  the six full production trajectories are unavailable in this workspace;
  no real-trajectory normalization or production report claim is made.
