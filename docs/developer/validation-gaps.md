# Validation gaps

**Historical snapshot (2026-11-03):** Python 3.12.14 / coverage.py 7.17.0; **409 tests passed**, two benchmarks excluded. Covers only the measured in-process pytest run—**not** subprocess CLI execution, release trajectories or the current test count. Coverage is not a scientific validity score.

| Module | Statements | Branches |
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

```bash
python -m pip install coverage
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage run --branch --source=src/mdscope -m pytest -q
COVERAGE_FILE=/tmp/mdscope-coverage python -m coverage report
```

## Hand-tested mutations

Six **isolated code-copy** changes were caught by numerical assertions, not import failures. This is **not** a global mutation score.

| Deliberate error | Test |
|---|---|
| Omit particle normalization in FFT MSD | `test_multi_origin_msd.py::test_fft_matches_independent_direct_sum` |
| Omit minimum-image step | `test_core_geometry_contracts.py::test_minimum_image_matches_independent_mdanalysis_reference` |
| All lipids assigned lower leaflet | `test_physical_validation.py::test_anchor_groups_vectorized_mapping` |
| Reverse isolated/environment SASA subtraction | `test_scientific_reference_oracles.py::test_sasa_burial_subtraction_and_exact_reference_identity` |
| Reverse paired bootstrap labels | `test_publication_references.py::test_nonconstant_paired_hierarchy_preserves_labels_and_draws` |
| Replace unavailable enrichment with zero | `test_scientific_reference_oracles.py::test_composition_reference_two_species_and_missing_denominator` |

Real FreeSASA sphere/cap geometry and fake-FreeSASA streaming parity verify **different** contracts.

## Historical parity

| Reference | Finding | Limit |
|---|---|---|
| Six-replica fixture vs `f7af0e1` | **42 identical Parquet tables** and next 100 RNG integers per family | Invalid/missing profile corrections intentionally differ; no performance claim |
| Interaction grouping vs `7fc100d` | Three identical tables and next 100 RNG integers; ~0.09 s/toy variant | 545,719 identities / 5,000 draws still unbenchmarked |
| 2026-11-13 suite refactor | Geometry/semantic oracles consolidated; duplicate multi-origin fixture removed but unique failure check retained | Historic 409 count is stale; no new coverage measurement |

## Highest-value work

| Gap | Next check |
|---|---|
| Subprocess CLI code not covered in parent process | Trace spawned executions before inferring missing verification |
| SDC4/HS chemistry and selection preflight | Validate real TPRs; toy topologies are not exhaustive |
| Long-event and large interaction memory | Measure peak RSS, output size, wall time on real six-replica work |
| Convergence/sensitivity | Inspect independent replicas, time windows, cutoff, RDF and diffusion behavior |
| Completed-stage artifact failure | Fault-inject missing/corrupt companion; require safe failure or rebuild |
| Irregular/missing time grids | Preserve partial valid observations, distinguish NaN from measured zero |
| Scientific interpretation | Challenge physical assumptions with independent reference models |

## Six failure modes during review

1. **Assumptions:** plausible numerical output, incorrect physical model.
2. **Mathematics:** sign, normalization, units, PBC and population errors.
3. **State/data flow:** checkpoints, committed parts, provenance.
4. **Integration:** identity/time/shape consistency, missingness.
5. **Verification:** tests must detect the intended regression.
6. **Recovery:** invalid/corrupt work must not destroy valid outputs.

Retain topology, contact/mass and crash-regression tests despite old issue-number filenames. Keep reporting tests while reporting code exists. [Scientific contracts](scientific-contracts.md) · [Release checklist](release-validation.md).
