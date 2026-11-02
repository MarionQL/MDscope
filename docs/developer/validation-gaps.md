# Validation effectiveness checkpoint — 2026-11-02

Measured on the statistics/validation branch, Python 3.12, coverage.py 7.17.0:
283 tests passed; two optional timing tests deselected. Coverage traces this
pytest process, **not subprocess CLI executions**. Source is `src/mdscope`;
large real-GROMACS scripts, generated files and optional timing workloads are
not part of this measurement. These fractions are a gap map, not a scientific
validity score or a percentage gate. Molecular/report branch changes require
their own integrated measurements after review.

| Scientific module | Statements covered / total | Branches covered / total |
|---|---:|---:|
| `observables/multi_origin_msd.py` | 69/81 | 22/34 |
| `core/geometry.py` | 11/11 | 4/4 |
| `core/lipid_geometry.py` | 101/117 | 23/40 |
| `observables/sasa.py` | 134/175 | 47/76 |
| `mechanisms/sdc4_mechanisms.py` | 484/649 | 137/258 |
| `statistics/common.py` | 175/193 | 43/58 |
| `statistics/interactions.py` | 156/164 | 63/74 |
| `config.py` | 180/455 | 84/328 |
| `reporting/report.py` | 367/537 | 124/244 |

Reproduce after installing the optional developer tool (`pip install coverage`):

```bash
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage run --branch --source=src/mdscope -m pytest -q
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage report
```

## Hand-reviewed scientific mutants

Temporary source changes were run in fresh pytest subprocesses, then original
source and compiled caches were restored. Each mutation below caused the
listed numerical/identity assertion to fail; collection/import errors were
not counted as detections. Six selected detections are not a global mutation
score. No equivalent or surviving mutation was inferred from this small set.

| Deliberate error | Detecting independent reference |
|---|---|
| Remove particle count from FFT MSD normalization | `test_multi_origin_msd.py::test_fft_matches_independent_direct_sum` |
| Omit minimum-image displacement | `test_core_geometry_contracts.py::test_minimum_image_matches_independent_mdanalysis_reference` |
| Classify every lipid in the lower leaflet | `test_physical_validation.py::test_anchor_groups_vectorized_mapping` |
| Reverse isolated minus environment SASA | `test_scientific_reference_oracles.py::test_sasa_burial_subtraction_and_exact_reference_identity` |
| Reverse replica order in only one paired bootstrap population | `test_publication_references.py::test_paired_hierarchical_bootstrap_constant_reference` |
| Replace unavailable enrichment by zero | `test_scientific_reference_oracles.py::test_composition_reference_two_species_and_missing_denominator` |

The burial-sign probe exposed missing direct coverage before this batch.
A quantitative subtraction/fraction oracle and wrong-time/wrong-residue joins
were added. Do not substitute fake-FreeSASA streaming parity for absolute
geometry: real sphere/cap/periodic occluder tests serve a separate purpose.

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
2. Many configuration branches are topology/run-specific: extend negative
   selections, run IDs, incompatible references and failure-before-output
   fixtures. Cheap invalid stride/chunk/time/CI/family/window cases are now
   parameterized, but that is not the full matrix in #74.
3. SASA periodic/heterogeneous radii, temporal missingness in profile statistics,
   and long-event finalization memory need further physical/scaling references.
4. Complete SDC4 major-observable failure/boundary fixtures and the six-replica
   production comparisons; branch coverage does not establish their sampling.
5. Keep historical test suites until their distinct behavior is mapped.
   Packaging/layout/import overlap was consolidated, numerical oracles retained,
   and timing tests moved out of the default gate. No mass renaming was done.
