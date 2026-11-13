# Validation effectiveness checkpoint — 2026-11-03

Measured on the consolidated inference/validation branch, Python 3.12.14,
coverage.py 7.17.0: **409 tests passed**; two optional timing tests deselected.
Coverage traces this
pytest process, **not subprocess CLI executions**. Source is `src/mdscope`;
large real-GROMACS scripts, generated files and optional timing workloads are
not part of this measurement. These fractions are a gap map, not a scientific
validity score or a percentage gate. The current mixed-observer, report,
nonconstant statistical and SDC4 reference suites are included.

| Scientific module | Statements covered / total | Branches covered / total |
|---|---:|---:|
| `observables/multi_origin_msd.py` | 69/81 | 22/34 |
| `core/geometry.py` | 11/11 | 4/4 |
| `core/lipid_geometry.py` | 101/117 | 23/40 |
| `observables/sasa.py` | 143/179 | 56/82 |
| `mechanisms/sdc4_mechanisms.py` | 514/667 | 154/266 |
| `statistics/common.py` | 194/209 | 52/64 |
| `statistics/time_series.py` | 134/152 | 52/66 |
| `statistics/profiles.py` | 149/191 | 49/70 |
| `statistics/interactions.py` | 156/163 | 62/72 |
| `config.py` | 347/522 | 211/386 |
| `core/validation.py` | 54/119 | 29/84 |
| `reporting/report.py` | 435/595 | 178/282 |

Reproduce after installing the optional developer tool (`pip install coverage`):

```bash
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage run --branch --source=src/mdscope -m pytest -q
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage report
```

**Test organization update (2026-11-13):** The issue #31 geometry and
cross-module semantic test modules were consolidated into
`test_core_geometry_contracts.py`. Distinct triclinic, periodic-center,
pairwise-versus-pointwise, and invalid-geometry checks remain. The duplicate
multi-origin FFT/fit integration fixture was removed; the direct-sum oracle,
end-to-end fit reference, and explicit invalid-window test remain in
`test_multi_origin_msd.py`. This is a structural change, not a new coverage
measurement; the 409-test figure above is historical. Standard PR CI still
runs Python 3.11/3.13 synthetic checks and validates wheel/sdist once on 3.11;
optional profiling is available through manual workflow dispatch. The separate
real-trajectory checks remain unchanged.

## Hand-reviewed scientific mutants

Six isolated temporary **copies** of the current package were tested in fresh
pytest subprocesses; the working source was not mutated. Each error below caused the
listed numerical/identity assertion to fail; collection/import errors were
not counted as detections. Six selected detections are not a global mutation
score. All six selected mutations were detected by assertions; there were no
survivors in this set to classify as equivalent or missing coverage. Untested
mutations elsewhere remain unknown.

| Deliberate error | Detecting independent reference |
|---|---|
| Remove particle count from FFT MSD normalization | `test_multi_origin_msd.py::test_fft_matches_independent_direct_sum` |
| Omit minimum-image displacement | `test_core_geometry_contracts.py::test_minimum_image_matches_independent_mdanalysis_reference` |
| Classify every lipid in the lower leaflet | `test_physical_validation.py::test_anchor_groups_vectorized_mapping` |
| Reverse isolated minus environment SASA | `test_scientific_reference_oracles.py::test_sasa_burial_subtraction_and_exact_reference_identity` |
| Reverse replica order in only one paired bootstrap population | `test_publication_references.py::test_nonconstant_paired_hierarchy_preserves_labels_and_draws` |
| Replace unavailable enrichment by zero | `test_scientific_reference_oracles.py::test_composition_reference_two_species_and_missing_denominator` |

The earlier burial-sign probe motivated the quantitative subtraction/fraction
oracle and wrong-time/wrong-residue joins retained here. Fake-FreeSASA streaming
parity and the real sphere/cap/periodic/heterogeneous-radius references serve
different contracts; neither replaces the other.

## Current default-result parity and new controls

Against `main` at `f7af0e1`, a six-replica fixture with unequal lengths,
nonconstant scalar/profiles/RMSF data, reordered profile columns, absent
interactions and censored lifetimes produced **42 exactly equal Parquet
tables** (values, schema, dtypes, row order). The next 100 integers from each
of the scalar, profile and interaction generators matched exactly; a separate
common-kernel RNG probe also matched. The original all-finite profile
reduction order was retained to avoid tiny floating-point drift.

The missing/irregular-profile correction intentionally changes intervals for
invalid data: unavailable positions cannot silently drop an invalid replica
from a system CI, and valid columns retain their temporal intervals. Point-only
and no-hierarchy modes have explicit zero-draw output metadata; hypothesis
tests are independent controls. Disabled kernels are checked before execution.
Previous statistics/method manifests are archived and partial runs are not
eligible for reporting. These are not wall-clock or production speedup claims.

## Interaction baseline parity

The existing grouping optimization was compared against historical source
`7fc100d` (before PR #90), using six replicas, unequal lengths, absent
interactions, censored/uncompleted lifetimes and both paired/unpaired designs.
All three output Parquet tables matched exactly in values, schemas, dtypes and
row ordering; the next 100 RNG integers matched exactly. The small fixture
uses three identity classes and 100 bootstrap draws. Timings were approximately
0.09 s per variant and establish **no material production speedup**. The
545,719-identity / 5,000-draw production runtime/RSS comparison remains open.
All-NaN lifetime bootstrap draws retain NaN and their existing warnings; they
were not replaced by zeros to obtain parity.

## Prioritized remaining gaps

1. Parent-process coverage is incomplete for CLI-driven interaction/finalizer
   paths; retain direct six-replica interaction references and add subprocess
   tracing before treating untraced lines as untested behavior.
2. Structural preflight now has shared CLI/programmatic references for empty
   selections, alignment size, overlapping anchors, SASA references, YAML/run
   IDs, paired labels and failure-before-output. The SDC4/HS parts of its local
   coverage map are limited; release-TPR validation is a separate real-data CI
   gate. This is not exhaustive validation of every possible topology/selector.
3. Long-event finalization peak memory and the full six-replica interaction
   runtime/RSS benchmark (#89) remain open. Passing small-fixture parity does
   not establish production scalability.
4. The named major-observable SDC4 positive/boundary manifest is in
   `sdc4-fixtures.md`. Cross-observer planar policies now have independent tilted/curved, supported
   triclinic and explicit unsupported-box references (see scientific contracts).
   Six-production-replica recomputation/convergence (#98–#109) remains unverified;
   branch coverage and analytic fixtures do not establish their sampling.
5. Historical issue-prefixed topology, contact/mass and glycolipid fixtures
   retain distinct numerical contracts. The overlapping issue #31 geometry
   and semantic checks were integrated into the geometry test module without
   discarding their distinct edge cases. No suite was dropped on the basis
   of a coverage percentage. New tests are grouped by scientific responsibility.

## Targeted AI-assisted failure-mode review (2026-11-03)

Use this as a **review checklist**, not a declaration of defects or a new coverage score. In particular, distinguish numerical tests from scientific assumptions and negative tests from merely successful execution. Relevant mechanisms already exist; strengthen them only where a demonstrated gap affects a real contract.

| Failure mode | Existing evidence and boundary | Useful adversarial check |
| --- | --- | --- |
| Incorrect assumptions | [Scientific audit](scientific-audit.md) records RDF reference-density, fit-window, leaflet and sampling limitations; explicit documentation does not resolve the underlying physical assumptions. | Compare an implementation to an independently specified toy system and challenge the physical interpretation, not just units or schema. |
| Scientific/mathematical errors | `test_scientific_reference_oracles.py`, `test_publication_references.py` and numerical-contract tests cover selected estimators. | Deliberately change a critical sign/denominator/unit/coordinate convention in an isolated test run and show the relevant reference catches it. |
| State and data flow | Streaming checkpoint identities, artifact inventory and fault-injection restart tests provide strong protection for committed chunks. `pipeline/runner.py:RunState.complete` checks explicit completion status; completed-stage reuse is not uniformly structure-validated. | After completion, inject a missing/corrupt required derived product or a noncomplete marker and verify safe recomputation or explicit failure without loss of valid results. |
| Integration | `_merge_timeseries` rejects conflicting metric columns and duplicate keys; mixed-observer restart parity tests exist. The outer time merge allows absent observations where legitimate. | Where sources promise identical sampling, remove one observation and require a coverage mismatch diagnostic; keep legitimate partial sampling distinct. |
| Incomplete verification | Independent numerical references and selected hand-reviewed mutant detections are documented above. These do not establish a diffusive regime, adequate six-replica sampling, or correctness of every supported topology. | Verify that a newly added test fails under the target regression; add production-data sensitivity or external comparisons only when relevant and available. |
| Failure handling | Input preflight, checkpoint preservation, and stale statistical/report rejection are tested. `finalize_project` checks persisted time axes and checkpoint counts, not completeness of every enabled observable. | Remove an enabled, required extension output while retaining a plausible base timeseries; assert missingness is detected and clearly reported, rather than silently presented as a complete result. |

These are **source-inspection observations and candidate regressions**, not evidence that currently persisted production results are wrong. Re-check at the implementation commit before writing an issue or claiming a fix. Prefer correcting the existing stage/manifest contracts over adding another state or validation framework. Do not silently broaden required-artifact rules to optional/disabled analyses. Preserve production outputs while exercising failures on temporary fixtures.
