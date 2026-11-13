# Release validation

## Portable numerical project

From an installed checkout (`python -m pip install -e '.[test]'`):

```bash
export MDSCOPE_STATS_WORKERS=1
python tests/validation_fixture.py create /tmp/mdscope-validation
mdscope validate /tmp/mdscope-validation/validation.yaml
mdscope run /tmp/mdscope-validation/validation.yaml
python tests/validation_fixture.py verify /tmp/mdscope-validation
```

- The destination must be **new**; the generator refuses to overwrite inputs.
- Keep generated PDB, DCD, YAML and `reference.json` together (YAML uses relative paths).
- Verification uses **saved numerical outputs**, without reopening trajectories.
- Eight synthetic frames; analyze **0,2,4,6**. The 45-atom model is not an equilibrated or chemically complete MD system.
- Independent values: `tests/data/validation_reference.json` and [SDC4 fixtures](sdc4-fixtures.md).

| Product | Expected independent result | Tolerance |
|---|---|---|
| Scalar time series | 0/2/4/6 ps, bilayer 32 Å, TM angle 0°, ectodomain 12 Å, HS 6 Å | 1e−5 |
| Bilayer normals | Flat plane, 0°; all 18 anchors supported | 1e−6 |
| Density | 9 anchors/occupied 5 Å bin, 9/50,000 Å⁻³, total 18 | 1e−12 |
| Local enrichment | 7/9 lipid anchors/leaflet in 15 Å; one-species fold 1 | 1e−12 |
| RDF pairs | 6 CA–water pairs/frame within 6 Å → 24 total | 1e−12 counts |
| Chain A residue distances | 20,16,12,8,4,4,8,12,16,12,8,4 Å | 1e−5 Å |
| Chain B distances | `sqrt(h²+(3+0.03×frame)²)` for A distances `h` | 1e−5 Å |

Density centers: −32.5/+32.5 Å (uncentered box-midpoint reference). RDF bins **at exact edges** can differ by backend rounding; check total counts. DCD float32/time arithmetic motivates tolerances.

The fixture covers fused structural/membrane/SDC4/solvation/distance observers and 100-draw statistics, **not** native FreeSASA/DSSP or between-replica uncertainty.

## Release checklist

- [ ] Record commit, Python/NumPy/MDAnalysis/FreeSASA versions, seed, worker allocation, commands, elapsed time, peak RSS, skips.
- [ ] Run `pytest -q tests` and the portable oracle; keep independent mathematical/PBC/chemical/failure references.
- [ ] Build wheel/sdist and install/test **outside checkout** for Python 3.11 and 3.13 ([Architecture](architecture.md)).
- [ ] Run native FreeSASA sphere/cap/periodic tests; document optional skips.
- [ ] Check crash/restart parity, companion artifact identity, and non-destructive input failures.
- [ ] Review coordinate provenance, units, topology identities, leaflet definitions, sampling/censoring and missing-vs-zero.
- [ ] Run release-backed GROMACS smoke/physics **when data/quotas permit**, then separately all six production replicas and sensitivity/convergence.
- [ ] Cite [`CITATION.cff`](../../CITATION.cff) with actual version/commit and comply with the repository [MIT License](../../LICENSE) (retain copyright and license notice).

**Local investigation budget:** ~1 minute / 512 MiB for the tiny reference project with `MDSCOPE_STATS_WORKERS=1`. These are guidance, **not** timed CI assertions or production throughput predictions. Optional profiler: `pytest -q -s -m benchmark tests/test_performance_benchmarks.py`.

Record resource-blocked checks as **unverified**, never passed. A short reference run cannot support manuscript-level claims.
