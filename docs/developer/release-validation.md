# Numerical reference project and release checklist

## Portable synthetic project

From a checkout with `python -m pip install -e '.[test]'`, run:

```bash
export MDSCOPE_STATS_WORKERS=1
python tests/validation_fixture.py create /tmp/mdscope-validation
mdscope validate /tmp/mdscope-validation/validation.yaml
mdscope run /tmp/mdscope-validation/validation.yaml
python tests/validation_fixture.py verify /tmp/mdscope-validation
```

Use a fresh destination: creation refuses to replace inputs, configuration,
reference values or results. The generated PDB, DCD, YAML and `reference.json`
can be moved together; input and output paths in YAML are relative to that
file. Verification reads saved numerical products without reopening trajectories.
It requires neither a release download nor GitHub Actions. The constructor is
shared with the mixed-observer restart tests, rather than maintaining a second
simulation fixture or a binary trajectory archive.

The 45-atom geometry is specified in [SDC4 fixtures](sdc4-fixtures.md): two
12-CA proteins, 18 lipid anchors, two HS sites and one water oxygen. Its eight
frames translate chain B laterally by 0.03 Å per frame; frames 0, 2, 4 and 6
are analyzed. It contains self-generated synthetic coordinates and no private
simulation/force-field assets. It is not an equilibrated or chemically complete
molecular system. PDB/DCD input exercises the ordinary CLI; real GROMACS
selection/topology validation remains separate.

`tests/data/validation_reference.json` specifies seven products. Expected values
come from the geometry below, not from copying an MDscope run. Unlisted columns
remain available and are not required to retain an exact software-version schema.
Compatible historical results are not rejected because of a source/config hash.

| Product | Independent reference | Absolute tolerance |
| --- | --- | --- |
| Replicate scalar table | Four times 0, 2, 4, 6 ps; leaflet separation 32 Å; parallel TM axes 0°; ectodomain heights 12 Å; HS extension 6 Å | 1e-5 in stated units |
| Bilayer normal table | Two flat planes, zero mean deviation, all 18 anchors valid with the specified 30 Å neighborhood | 1e-6 |
| Lipid density profile | Nine anchors in each occupied 5 Å bin; each slab is 100×100×5 Å³, giving 9/50,000 Å⁻³ and integrated count 18 | 1e-12 |
| Lipid enrichment counts | Seven of nine POPC anchors per leaflet in the 15 Å protein-centered disk; one-species composition fold 1 | 1e-12 |
| Solvation RDF table | Six CA–water pairs per sampled frame within 6 Å, total 24 across four frames | 1e-12 for total count |
| Chain A residue-distance matrix | Distances to planes at ±16 Å: 20,16,12,8,4,4,8,12,16,12,8,4 Å at every frame | 1e-5 Å |
| Chain B residue-distance matrix | For each height distance h above, sqrt(h² + (3 + 0.03×frame)²) Å | 1e-5 Å |

The distances use the x=y=0 lipid anchor as the nearest point in these frames.
The two density bins have centers −32.5 and +32.5 Å because uncentered density
uses the box midpoint as origin. The RDF check uses a total count: individual
bins at exact distance boundaries may differ by backend rounding. Coordinate
and time tolerances allow float32 DCD/distance arithmetic; count references do
not permit missing or additional contacts. A successful file read alone fails
to establish any of these numerical values.

The generated YAML enables the existing fused structural, membrane, SDC4,
solvation and residue-distance observers, plus scalar statistics. It disables
plots, native FreeSASA/DSSP and large optional profile bootstraps. The
100-draw, seed-2026 statistics exercise finalization; one synthetic replicate
cannot establish independent-replicate uncertainty or scientific convergence.
Numerical statistical references live in `test_publication_references.py`,
`test_interaction_statistics_reference.py`, `test_statistics_modes.py` and
`test_parallel_bootstrap.py`; they are not replaced by this end-to-end fixture.

## Release evidence ledger

Record the commit, Python/NumPy/MDAnalysis/native-library versions, seed,
worker allocation, command, elapsed time, peak RSS and any skipped checks.
A commit is provenance, not a compatibility lock on saved outputs. Fixed seeds
are reproducible in a fixed environment; cross-version validation compares
numerical estimates within their stated tolerances, not unpinned RNG bitstreams.

- Run the ordinary synthetic suite and the reference commands above. Existing
  restart tests inject failures during a frame, after chunk writes and before
  checkpoint commit, then compare finalized mixed-observer products.
- Build wheel/sdist, check metadata and install each outside the checkout using
  the commands in [architecture](architecture.md). Existing workflows specify
  Python 3.11/3.13; a local 3.12 run does not certify that matrix.
- Run native FreeSASA analytical references when that optional library is
  available; record skips. Preserve analytical MSD/PBC/area/normal/contact,
  direct ACF/permutation/FDR/bootstrap and profile/missing-data references.
- Review coordinate provenance, units, selection identities, leaflet cohorts,
  sampling/censoring, replicate-level comparisons, and missingness. Scalar and
  profile temporal intervals are suppressed for unsupported gaps/irregularity;
  descriptive finite-observation means do not eliminate missingness bias.
- When resources are explicitly available, run the release-backed GROMACS
  smoke/physics scripts listed in `.github/workflows/real-trajectory-smoke.yml`,
  then the separate six-production-replica sensitivity/convergence review.
  Retain their numerical/resource evidence, not just a successful exit status.
- Use [CITATION.cff](../../CITATION.cff) with the actual version and commit.
  This repository currently declares no redistribution license. Record the
  owner's license decision before a public distribution; citation metadata
  neither grants a license nor chooses one. No release is authorized by this
  checklist.

For the small reference project, use a review budget of one minute and 512 MiB
peak RSS on an ordinary development CPU, with `MDSCOPE_STATS_WORKERS=1`.
These are investigation thresholds, not flaky timing assertions or a production
throughput promise. Optional kernel benchmarks remain opt-in (`pytest -m
benchmark -s tests/test_performance_benchmarks.py`). Real-data downloads and
full-production analysis have a separate resource budget; do not trigger them
as part of routine fixture generation or offline verification.

For the current quota-limited work, GitHub Actions, the Python 3.11/3.13 matrix,
release-backed jobs and six-production-replica checks are intentionally deferred
at the owner's instruction. Document numerical implementation completion and
those unverified release gates separately.
