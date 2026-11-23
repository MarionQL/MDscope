# Scientific audit

Numerical oracles and short release tests exist; **six-production-replica convergence and biological conclusions remain unverified**. A passing numerical test is not evidence of physical mechanism.

## Priority findings

| Area | Established | Still unresolved |
|---|---|---|
| Lipid specificity | Matched local contact numerator/denominator corrected (PR #117); old unmatched metric retained separately | Six-replica effect and lipid-size/area/cutoff sensitivity |
| Residence | Observed span, interval-imputed duration, frame occupancy, censoring separated | Microscopic lifetime between saved frames; irregular-time occupancy |
| Solvation RDF | Correctly labeled **full-cell** reference; optional aqueous-slab reference | Solvent-excluded shell volume, real bulk plateau |
| SASA | Isolated-atom FreeSASA crash corrected analytically; sphere/cap oracles | Radii, topology element assignment and whole-target PBC |
| Shared states | Replica-balanced k-means; missing-feature/gap diagnostics | k/seed/feature stability, metastability, causality, nonrandom missingness |
| Time statistics | Invalid/gapped/irregular traces retain descriptive means but not unsupported temporal CIs | Missingness bias, low independent-replica count |
| Solvation memory | Event extraction streamed; survival risk-count computation batched | Other bridge/survival/SASA derivatives still materialize tables |
| Leaflets and MSD | Fixed cohorts distinguished from dynamic leaflets; direct-sum/PBC/drift tests | Flip-flop, stride aliasing, long-lag fit and finite-size effects |
| Membrane geometry | Planar/PBC reference tests and unsupported-box errors | Strongly curved/tilted cells, headgroup-height mismatch |
| Local bilayer normals | Whole-box-shift PBC invariance on every real-data frame; hard-cutoff ties (a neighbour within 2.5e-4 Å, ~10× measured float32 shift error, of the 15 Å cutoff) reported separately (#250) | Per-lipid tilt is discontinuous at the cutoff: one neighbour crossing it moved a ~16-neighbour tilt by ~1° (carb test data, 250 ns). Ensemble means are robust; per-lipid tilts and their tails are cutoff-sensitive |
| Lipid tail naming | Tail mappings checked against topology elements/bonds; per-species/tail coverage in `mdscope validate`; corrected default and optional `charmm36` preset use sn-1=C3x, sn-2=C2x (#272) | **Pre-#272 saved outputs retain reversed identities** (legacy `sn1`=C2x/sn-2, `sn2`=C3x/sn-1); rerun in fresh output directory. PSM/ceramide sphingosine and n_acyl tails are mapped separately by the opt-in charmm36 preset (#276), not pooled into sn chains |
| H bonds and bridges | Bond-verified molecule IDs, geometry thresholds | Protonation/chemical assignments; proximity is not cooperativity |

## Other model assumptions

- **RDF:** per-bin replica mean/SD descriptive; selected solvent divided by full periodic-cell volume, **not** accessible water volume.
- **SASA:** geometric element radii, not CHARMM LJ; burial is isolated minus environment, nm².
- **Footprint:** projected free grid area, not molecular free area or affinity; check resolution, exclusion radius, slab width.
- **Solvation:** unique solvent entities per target, whereas RDFs pool atom-site pairs; site sizes affect counts.
- **Density:** variable NPT cells, split centers and oblique box axes can violate a nominal slab reference.
- **Structure:** RMSD/RMSF depend on fitting, Rg on masses, DSSP on geometry; proline regions are angle bins.
- **SDC4:** NxIP/GGIVG motifs, six HS attachments, GM1a/BGAL `chainID`/CER160 identities are topology-specific.
- **Glycans:** discovery trusts force-field residue names and bonds; only linear chains are numbered. Non-SDC4 glycan chemistry is validated on synthetic topologies only, and `O*` names are assumed to be the oxygen hydration sites.
- **Legacy outputs:** readable tables can predate numerical corrections. Missing corrected fields remain NaN.

Definitions and numerical contracts: [Scientific contracts](scientific-contracts.md) · [Analyses](../analyses.md).

## Verification coverage

CI uploads an `mdscope-coverage` report on the Python 3.11 job; line coverage is not a scientific-validity score.

### Hand-tested mutations

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

### Highest-value remaining checks

| Gap | Next check |
|---|---|
| Subprocess CLI code not covered in parent process | Trace spawned executions before inferring missing verification |
| SDC4/HS chemistry and selection preflight | Validate real TPRs; toy topologies are not exhaustive |
| Non-SDC4 glycan topologies | Audit a real N-/O-glycoprotein TPR/PSF (naming, linkage atoms, branches) |
| Long-event and large interaction memory | Measure peak RSS, output size, wall time on real six-replica work |
| Convergence/sensitivity | Inspect independent replicas, time windows, cutoff, RDF and diffusion behavior |
| Completed-stage artifact failure | Fault-inject missing/corrupt companion; require safe failure or rebuild |
| Irregular/missing time grids | Preserve partial valid observations, distinguish NaN from measured zero |
| Scientific interpretation | Challenge physical assumptions with independent reference models |

## Six-replica acceptance on ARC

Use a **fresh output directory** for changed estimators. Do not replace historical results.

```bash
git pull origin main
python -m pip install -e '.[test,sasa]'
python - <<'PY'
from pathlib import Path
import yaml
cfg = yaml.safe_load(Path("synd4_full_mechanisms.yaml").read_text())
cfg["project"]["output"] = "results_scientific_validation_<date>"
Path("synd4_scientific_validation.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
PY
mdscope validate synd4_scientific_validation.yaml
mkdir -p slurm-out/mdscope
sbatch slurm/run_mdscope.sh synd4_scientific_validation.yaml
# Submit slurm/finalize_mdscope.sh with an afterok dependency on the returned job ID.
```

Review across **all six** trajectories:
1. Matched/legacy enrichment, anchor populations, cutoff/footprint sensitivity and per-replica agreement.
2. Observed/imputed duration by censor state and physical stride; compare time-based gap limits.
3. RDF aqueous references, hydration/SASA, K steric/null controls and shell stability.
4. Shared-state missing fractions, seed/k/feature sensitivity and transitions.
5. Diffusion lag windows, initial/dynamic leaflets, trajectory provenance.
6. Full event-table/finalizer peak memory, runtime and interaction scaling.

Do not pool frames as independent replicas. Historical issue closure does not replace unverified six-replica evidence.
