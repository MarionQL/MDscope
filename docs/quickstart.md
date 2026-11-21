# Synthetic quickstart

> **Synthetic fixture, not MD evidence.** These trajectories are generated
> coordinates of a 45-atom toy model with a built-in offset. They exercise the
> pipeline; they say nothing about any protein, membrane or biological system.

Runs install → two conditions × three replicas → `findings.md` in seconds,
without trajectories or the real-data release.

```bash
git clone https://github.com/MarionQL/MDscope.git && cd MDscope
python -m pip install -e '.[test]'
python tests/validation_fixture.py quickstart ../mdscope-quickstart
mdscope validate ../mdscope-quickstart/quickstart.yaml
mdscope run ../mdscope-quickstart/quickstart.yaml
python -c "import pandas as pd; print(pd.read_parquet('../mdscope-quickstart/results/project_summary/statistics/findings.parquet').drop(columns='source_table').to_string())"
```

- **A checkout is required.** The generator reuses the analytic test geometry
  (`tests/test_sdc4_geometry_references.py`) and is deliberately not packaged,
  so the wheel ships no test fixtures. `[test]` provides `pytest` for that import.
- Run from the checkout root. Output goes beside the checkout, in a new or
  empty directory; the YAML uses paths relative to it.
- `tests/test_synthetic_quickstart.py` executes this block (after its clone and
  install lines) in the automatic Python checks.

## Design

Both conditions share `synthetic.pdb`. Each replica has 20 frames in which
chain B is translated along x, away from chain A, by the condition offset plus
seeded jitter in [0, 0.2) Å:

| System | Offset | Expected |
|---|---|---|
| `synthetic_reference` | 0 Å | CA-chain separation ≈ 3.1 Å |
| `synthetic_displaced` | 1 Å | separation ≈ 4.1 Å |

Only alignment, RMSD, RMSF and radius of gyration run. Project name
`synthetic_quickstart_fixture` and the system names mark every output as synthetic.

## Expected outputs

| Path under `../mdscope-quickstart/` | Content |
|---|---|
| `results/project_summary/statistics/findings.md` | Coverage per family and the comparisons whose interval excludes zero |
| `results/project_summary/statistics/findings.parquet` | Those comparisons, one row per estimator |
| `results/project_summary/statistics/coverage/contrast_coverage.parquet` | Every registered comparison, including `includes_zero` and unavailable rows |
| `results/project_summary/statistics/time_series/system_comparisons.parquet` | Full time-series contrasts with permutation p/q values |
| `results/synthetic_displaced/rep1/summary/timeseries.parquet` | Per-frame `time_ps`, `rmsd_A`, `radius_of_gyration_A` for one replica |

## Interpret the Parquet tables

`findings.parquet` columns: `family`, `observable` (JSON identity, here
`metric`), `system_1`/`system_2`, `estimator` (`mean`/`median`), `status`,
`difference_system2_minus_system1` (source units), `ci_low`/`ci_high`,
`ci_level`, `p_value`, `q_value`, `source_table`.

| Row | Expected | Why |
|---|---|---|
| `radius_of_gyration_A`, both estimators | `above_zero`, ≈ +0.058 Å | Analytic: `sqrt(var_z + (4.1/2)²) − sqrt(var_z + (3.1/2)²)`, with `var_z` the variance of the 24 CA z coordinates |
| `rmsd_A` (coverage only) | `includes_zero` | RMSD is relative to each replica's first frame, which already carries the offset |
| RMSF profiles (coverage only) | `includes_zero` | Per-atom fluctuations come from the same jitter in both conditions |
| `p_value` 0.1, `q_value` 0.2 | Floor, not weak data | Three versus three replicas allow only 20 permutation assignments; FDR then adjusts |

An interval excluding zero is a pointwise statement about these replicas. For
real projects follow the [interpretation checks](../README.md#interpret-cautiously),
the [result paths](workflows.md#find-and-interpret-results) and, for numerical
verification against independent references, the
[portable validation project](developer/release-validation.md#portable-numerical-project).
