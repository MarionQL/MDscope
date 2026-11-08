# Scientific methodology audit — 2026-11-02

This is a code-level methodology audit with independent numerical references,
release-backed CI evidence and targeted corrections. It is **not** validation
of equilibrium sampling or biological conclusions across the six production
replicas. Source paths below are relative to `src/mdscope/`; function names
identify the evidence. An unresolved limitation is not an acceptance pass.

## Prioritized findings and disposition

| Severity | Observable / source evidence | Definition and assumption | Disposition / scientific impact |
|---|---|---|---|
| High | Lipid specificity, `mechanisms/sdc4_mechanisms.py:update/finalize` | Unique contacting molecules versus anchor-in-disk availability; ratios are dimensionless. Numerator and denominator originally had different membership. | Corrected in PR #117; matched numerator cannot exceed availability. Whole-leaflet and legacy local outputs remain distinct. Four-frame real examples change other species' contact shares even when only POPC/cholesterol membership changes. Full production effect remains unknown. |
| High | Lipid/water/ion residence, `core/residence.py`, `mechanisms/motif_residence.py`, `observables/solvation.py:contact_episodes` | Contact spans in ps, interval-imputed estimates, frame occupancy and boundary censoring are different quantities. Irregular sampling does not make frame occupancy a time fraction. | Shared duration convention and explicit columns; backward-compatible historical durations. Batched time-axis processing removes repeated scans. Continuous/intermittent oracles pass. Interval-censored microscopic on/off times remain unresolved; no physical binding lifetime claim. |
| High | Solvation RDF and 3D density, `observables/solvation.py:update/finalize` | Pair counts divided by spherical shell volume and sum of target-atom count × selected-solvent count / instantaneous cell volume. It is a *full-cell average* reference, not measured aqueous bulk. | Corrected metadata, documentation and figure labels; historical `g_bulk` values/names preserved. Independent hemisphere example shows excluded volume and reference-density bias can cancel. A true bulk-water-region reference and excluded-volume estimator remain unresolved. Do not rank intrinsic hydration preference from g(r) alone. |
| High | Native SASA, `observables/sasa.py:sasa_residue_rows` | Probe-inflated sphere union; areas in nm². Isolated spheres have exact area 4π(r+probe)². | A real FreeSASA 2.2.1 Shrake–Rupley crash was reproduced for isolated atoms, including isolated+connected mixtures. Isolated atoms are now evaluated analytically and removed from native input; connected atoms still use the configured algorithm. Real sphere/cap references run in Python 3.11 CI. |
| High | Missing state-feature frames, `mechanisms/shared_states.py:fit_shared_states` | State populations condition on finite feature vectors; removed rows are unknown, not adjacent observations. | Report input/complete counts and missing fractions; exclude transitions across missing gaps. Suppress temporal block intervals and circular-shift nulls for incomplete or irregular traces. Population bias from nonrandom missingness remains unresolved; means and replica intervals are explicitly conditional. |
| High | Shared conformational model, `mechanisms/shared_states.py:_lloyd/fit_shared_states` | Global replica-balanced k-means of standardized geometric features; geometric clusters, not proven metastable states. k, scaling and feature choice define the partition. | Existing deterministic fit and replica weighting retained. Transition/half-window diagnostics and conditional samples are available. Multi-seed, feature-subset, k-sensitivity and lagged persistence analysis on production data remain required. No mechanistic or causal conclusion from state-conditioned contacts. |
| High | General temporal statistics, `statistics/time_series.py:analyze_time_series`, `statistics/common.py:finite/autocorrelation/_temporal` | Nonfinite observations are removed before ACF/block resampling. Retained rows can bridge missing times. Automatic positive-sequence ACF uses a finite lag window. | **Unresolved** for gapped/irregular scalar series: effective sample size and temporal intervals may be optimistic. Existing direct-sum ACF, exact permutation and hierarchical-bootstrap references establish numerical correctness only for their stated regular-data contracts. Validate regular sampling and missingness before interpreting these intervals; do not apply the shared-state temporal guard as a claim that all statistics are fixed. |
| High | Long-run solvation finalization, `observables/solvation.py:finalize`, `observables/solvation_postprocess.py:kaplan_meier` | Frame updates/chunk writes stream; finalization still concatenates persisted events/episodes. Survival risk counts scan groups per distinct duration. | **Unresolved scaling limit**: input streaming does not imply bounded finalization RAM. Release benchmarks cover short assets, not millions of identity events. Use class-restricted tracking and benchmark long kinetics runs before allocating production jobs. Pair-wise streaming episode reduction is a separate correction. |
| Medium | Replica RDF coverage, `statistics/solvation.py:analyze_solvation_rdfs` | Equal-weight replica mean and SD at each radial bin, not independent-bin inference. | Corrected coverage denominator to include every requested replica, including an entirely missing file/class. New independent missing-replica regression. Class absent for physical reasons remains distinguishable from complete data coverage only through selections/metadata. |
| Medium | Redundant contact queries, `mechanisms/sdc4_mechanisms.py:update` | Motif occupancy, occlusion and specificity use identical motif/lipid atoms, box and cutoff in the same frame. | Reuse the sparse motif pairs/distances for specificity, chemistry and residence. TM, HS, hydration and protein-interface searches retain their different selections/definitions. Seven synthetic observer Parquet outputs are numerically unchanged versus the pre-reuse code. |
| Medium | Leaflets, `observables/membrane.py:initialize_membrane`, `core/lipid_geometry.py:classify_leaflets/aligned_leaflets` | Legacy APL/MSD uses fixed initial median cohorts; enrichment/SDC4 uses circular instantaneous labels aligned to initial labels. | These are intentionally different populations. Released control data include a cholesterol leaflet change; never equate initial-cohort diffusion with instantaneous-leaflet diffusion. Median initialization assumes a centered, resolved bilayer and can fail for asymmetric/wrapped systems. General cohort and surface configuration remains unresolved. |
| Medium | Diffusion, `observables/msd.py`, `observables/multi_origin_msd.py`, `observables/multi_origin_pipeline.py` | Lateral MSD in Å²; D = slope × 2.5 in nm²/ns. FFT sums average overlapping origins; source coordinates are not protein-fitted. | Existing direct-sum, drift, PBC and source-reader references retained; new exact linear-curve unit oracle. Minimum-image intersample reconstruction assumes displacements are unaliased at the saved stride; intra-residue COM also assumes whole selected residues. Long-lag fit stability, collective drift choices and hydrodynamic finite-size bias remain production interpretation limits. Origins are not independent replicas. |
| Medium | SASA chemistry and environment, `observables/sasa.py:sasa_instances/sasa_residue_rows` | Explicit elemental radii/probe/occluders; periodic mode requires a compact whole target. Heterogeneous elemental radii are a geometric model. | Target/reference compatibility and exact topology/time joins already enforced. Real isolated-sphere and two-sphere cap tests added. CHARMM-name element inference is not a universal chemistry classifier; validate topology elements/radii before cross-force-field use. Reuse only when selections and radii/probe/sampling match. |
| Medium | Projected accessible area, `mechanisms/sdc4_mechanisms.py:_accessible_area_grid` | Count free xy grid cells after protein-disc exclusion inside a leaflet slab; nm² after conversion. | Analytic grid/PBC oracles retained. This is not exact lipid free area or solvent-accessible area. Grid spacing, exclusion radius and slab width require sensitivity checks. Matched composition uses molecular counts, not this area as an affinity correction. |
| Medium | Hydrogen bonds, `interactions/hydrogen_bonds.py`, `interactions/lipid_hydrogen_bonds.py`, `observables/solvation.py:water_donor_hbonds` | Donor–acceptor distance and donor–hydrogen–acceptor angle thresholds. Lipid identity follows bonded molecules; same-molecule events are excluded for intermolecular statistics. | Real topology parity and angular/PBC tests retained. Protonation, charge-based guessed chemistry and explicit donor/acceptor selections remain force-field inputs. Water bridge mode tests water-donor geometry; it is not a general all-direction H-bond bridge detector. Coordination bridges are proximity only. |
| Medium | Solvation counts/shells, `observables/solvation.py:unique_target_solvent_pairs/_targets`, `observables/solvation_postprocess.py:suggest_shell_boundaries` | Unique solvent entities per target within nearest selected-site cutoff. RDF pools atom-site pairs and is a different weighting. | Existing unique-count, overlap, ion-species and RDF-minimum tests retained. Cutoffs are explicit, minimum suggestions are not automatically imposed, and atom-count/target-size differences affect cross-species comparisons. Compare counts with target size, SASA and shell sensitivity. |
| Medium | Density profiles, `observables/density_1d.py:DensityProfileAccumulator` | Number or mass per physical slab volume; bins are initialized from the first frame and centering can use arithmetic COM/COG. | **Unresolved** for strongly varying NPT axis lengths, triclinic axis projections or a center group split across PBC. Check integrated density/count conservation and wrapped-center invariance on the intended input. Existing nominal-box tests do not establish these general cases. |
| Medium | APL/thickness/local normals/order, `observables/membrane.py`, `observables/bilayer_normal.py`, `observables/membrane_order.py` | Projected box area/leaflet count; anchor height differences; local height-field slopes; selected C–H orientational order. | These operational membrane metrics assume a resolved, approximately planar membrane with comparable anchors. SDC4 explicitly rejects nonorthorhombic boxes and non-z normals. Strong curvature, vesicles and mixed anchor heights need a different surface model; no general membrane-shape claim. |
| Medium | Structural profiles, `core/alignment.py`, `observables/rmsd.py`, `observables/rmsf.py`, `observables/radius_of_gyration.py`, `observables/dssp.py` | RMSD/RMSF depend on fit selection/reference; Rg is mass-weighted; DSSP is a secondary-structure assignment. | Existing alignment/translation/structure regressions retained. Whole-molecule imaging, correct masses and missing-residue mappings remain input contracts. Proline Ramachandran region labels are local angle bins, not proof of a helix/turn transition. |
| Medium | Biological mapping, `mechanisms/sdc4_mechanisms.py`, `mechanisms/shared_states.py` | GGIVG/NxIP/VVHPLVPL and A/B-chain interpretation are SDC4-specific; glycolipid chain IDs and CER160 chemistry are topology-specific. | Keep these in mechanism/configuration layers. Numerical geometry, residence, SASA and MSD kernels are reusable. Do not present SDC4 conditionals as generic protein-motif code; some chain/first-match assumptions still limit transfer. |
| Low | Output/provenance, `pipeline/runner.py:initialize`, `reporting/report.py`, `reporting/figure_inventory.py` | Checkpoint structure and replica identity protect restart; PDF is a bounded presentation with full numerical bundle. | Corrected metrics now reach figures/tables and full bundle; legacy names remain explicit. Missing old fields stay NaN. Figure presence is coverage, not scientific validation. Removed one unused duplicate COM helper from `observables/membrane.py`. |

## Numerical and performance evidence

- PR #117 historical checkpoint `1110e89`: Python 3.11/3.13, real GROMACS
  smoke and release solvation benchmark passed. Four mechanism frames per
  topology support the enrichment impact numbers in [scientific contracts](scientific-contracts.md).
- PR #117 follow-up `936dfe5`: 257 local tests passed. Release solvation
  benchmark passed: 120-frame standard runs used 1,074 MiB and 1,056 MiB peak
  RSS for glycolipid/control topologies (~1.02 million atoms); update rates
  were 4.857 and 4.854 fps. These measure the solvation subset, not the full
  production analysis. Native allocation and topology-reader memory dominate.
- Duration postprocessing: 100,000 episodes over 2,000 observations, 4.714 s
  before / 0.815 s after, 47.30 / 48.58 MiB Python allocation peak. Numerical
  output equality checked. This is 5.79x for the timed postprocessing workload.
- Motif search reuse: deterministic 80,000-atom lipid cloud, eight 40-atom
  motifs, ten frame evaluations, 5 Å cutoff. 160 searches become 80;
  2.384 s becomes 1.105 s (2.16x). Python allocation peak grows 0.617 → 0.892
  MiB because sparse pairs are kept until aggregation. Exact contact/distance
  checksum retained. Reproduce with `python scripts/benchmark_contact_reuse.py`.
  This is a query microbenchmark, not a full-trajectory speedup claim.
- SASA crash protection adds a neighbor query. A connected 1,000-carbon grid,
  1,000 surface points, one thread: minimum of three evaluations increased
  0.0245 → 0.0325 s; areas agree to 1e-12 nm². The extra ~8 ms is accepted for
  exact isolated-atom handling and native-crash prevention. No whole-run RSS
  or speed claim is made from this small benchmark.

Independent references include minimum-image contacts, chemical moieties,
matched species shares, irregular residence intervals, right/left censoring,
hemisphere RDF/reference-density cancellation, real SASA sphere/cap geometry,
MSD direct sums/unit conversion, analytic order/normal geometry, exact
permutation enumeration, BH correction and deterministic hierarchical draws.
Passing these establishes the tested estimators, not adequate trajectory sampling.

## Remaining production validation and exact launch commands

The six production trajectories at `/projects/bevanlab/...` are unavailable in
this execution environment. The downloaded CI artifacts contain short results,
not those trajectories. No production convergence or new biological conclusion
has been claimed. Keep issues #100–#102 open for outstanding scientific evidence.

On ARC, after the validated changes are merged and installed, create a fresh
output configuration so historical outputs remain intact:

```bash
git pull origin main
python -m pip install -e '.[test,sasa,report]'
python - <<'PY'
from pathlib import Path
import yaml
source = Path('synd4_full_mechanisms.yaml')
cfg = yaml.safe_load(source.read_text())
cfg['project']['output'] = 'results_memb_scientific_validation_20261102'
Path('synd4_scientific_validation.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
PY
python -m mdscope.cli validate synd4_scientific_validation.yaml
mkdir -p slurm-out/mdscope
validation_job=$(sbatch --parsable run_mdscope.sh synd4_scientific_validation.yaml)
sbatch --dependency="afterok:${validation_job}" finalize_mdscope.sh synd4_scientific_validation.yaml
```

Finalization consumes completed replica outputs for summaries/figures; do not
substitute the finalize-only YAML when requesting a new estimator recomputation.
The resolved YAML and available run metadata are part of the evidence bundle.

After completion, compare legacy versus matched enrichment within each replica,
region, leaflet and species, including valid-frame counts and contact-population
mismatch frequency. Compare observed-span and imputed durations separately,
stratified by censor flags. Do not pool frames across replicas. Repeat at
prespecified reasonable contact cutoffs, local radii and trajectory strides in
separate output directories. Frame-gap tolerance changes physical gap time
when stride changes; align tolerances in time before comparing kinetics.

Review first/second-half estimates, missing-feature rates, transition counts,
condition-specific state support and replica consistency. Repeat k/seed/feature
choices before mechanistic interpretation. Audit general scalar-series temporal
intervals on any gapped output, long-event finalization memory and physical
bulk-water references before treating the remaining high-severity findings as
resolved. They are explicit research/implementation work, not silently waived.
