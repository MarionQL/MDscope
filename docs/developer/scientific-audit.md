# Scientific audit

**2026-11-02 checkpoint:** numerical oracles and short release tests exist; **six-production-replica convergence and biological conclusions remain unverified**. A passing numerical test is not evidence of physical mechanism.

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
| H bonds and bridges | Bond-verified molecule IDs, geometry thresholds | Protonation/chemical assignments; proximity is not cooperativity |

## Other model assumptions

- **RDF:** per-bin replica mean/SD descriptive; selected solvent divided by full periodic-cell volume, **not** accessible water volume.
- **SASA:** geometric element radii, not CHARMM LJ; burial is isolated minus environment, nm².
- **Footprint:** projected free grid area, not molecular free area or affinity; check resolution, exclusion radius, slab width.
- **Solvation:** unique solvent entities per target, whereas RDFs pool atom-site pairs; site sizes affect counts.
- **Density:** variable NPT cells, split centers and oblique box axes can violate a nominal slab reference.
- **Structure:** RMSD/RMSF depend on fitting, Rg on masses, DSSP on geometry; proline regions are angle bins.
- **SDC4:** NxIP/GGIVG motifs, six HS attachments, GM1a/BGAL `chainID`/CER160 identities are topology-specific.
- **Legacy outputs:** readable tables can predate numerical corrections. Missing corrected fields remain NaN.

Definitions and numerical contracts: [Scientific contracts](scientific-contracts.md) · [Analyses](../analyses.md).

## Historical performance evidence

| Evidence | Observation | Cannot establish |
|---|---|---|
| PR #117, `1110e89` | Real release smoke, 4 mechanism frames/topology | Production condition effect |
| `936dfe5` | 257 local tests; 120-frame solvation, ~1,074/1,056 MiB RSS and ~4.85 fps | Whole MDscope production speed |
| 100,000 episodes | Postprocessing 4.714 → 0.815 s (5.79×); 47.30 → 48.58 MiB Python allocations | Whole-process RSS, physical lifetimes |
| Sparse motif search | 160 → 80 queries; 2.384 → 1.105 s (2.16×) | Full-run speedup |
| 1,000-carbon SASA | 0.0245 → 0.0325 s with crash prevention | Whole-run cost |

References include PBC, contact populations, SASA sphere/cap, sampled residence, RDF excluded-volume counterexamples, MSD and exact permutation/bootstrap arithmetic. Synthetic results do not certify convergence.

## Six-replica acceptance on ARC

Use a **fresh output directory** for changed estimators. Do not replace historical results.

```bash
git pull origin main
python -m pip install -e '.[test,sasa]'
python - <<'PY'
from pathlib import Path
import yaml
cfg = yaml.safe_load(Path("synd4_full_mechanisms.yaml").read_text())
cfg["project"]["output"] = "results_memb_scientific_validation_20261102"
Path("synd4_scientific_validation.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
PY
mdscope validate synd4_scientific_validation.yaml
mkdir -p slurm-out/mdscope
sbatch run_mdscope.sh synd4_scientific_validation.yaml
# Submit finalize_mdscope.sh with an afterok dependency on the returned job ID.
```

Review across **all six** trajectories:
1. Matched/legacy enrichment, anchor populations, cutoff/footprint sensitivity and per-replica agreement.
2. Observed/imputed duration by censor state and physical stride; compare time-based gap limits.
3. RDF aqueous references, hydration/SASA, K steric/null controls and shell stability.
4. Shared-state missing fractions, seed/k/feature sensitivity and transitions.
5. Diffusion lag windows, initial/dynamic leaflets, trajectory provenance.
6. Full event-table/finalizer peak memory, runtime and interaction scaling.

Do not pool frames as independent replicas. Historical issue closure does not replace unverified six-replica evidence.
