# ADR-0002: Evidence-first architecture and cohesive modules

**Status:** accepted design criteria and preferred direction (2026-11-22); implementation decisions require measured evidence. **Supersedes ADR-0001's mandatory directory migration**, not its scientific or dependency-contract goals. **Tracker:** [#269](https://github.com/MarionQL/MDscope/issues/269); execution work: [#279](https://github.com/MarionQL/MDscope/issues/279).

**Implementation note (issue #269, P0/P1):** The onboarding baseline and exemplar traces are in [analysis-onboarding-baseline.md](../analysis-onboarding-baseline.md). The proposed boundary fix centralizes enabled selectors in `config/runs.py`, moves topology-specific preview out of core alignment, and makes persisted scalar sources contract-owned; historical module inventories below intentionally retain their original snapshot identifiers. Do not confuse these import changes with a measured performance improvement.

## Why the earlier refactors felt incomplete

Moving a module is not inherently an architectural improvement. MDscope's highest costs arise when a scientific analysis needs many separate changes to become runnable, resumable, discoverable in Studio, statistically selectable, and correctly represented in figures. The previous ADR prescribed an entire `analyses/`, `engine/`, `interfaces/` layout before proving it would reduce that cost. Conversely, micro-modules and compatibility-only forwarding increased navigation without creating independently owned behavior.

**Decision:** choose the smallest cohesive package structure that measurably improves analysis-change locality, execution/recovery ownership and dependency direction. Preserve distinct scientific estimands and coordinate semantics. Folder names, lines per file, Ruff C901, and number of classes are *not* success criteria.

## Evidence and scope

- **Structural inventory:** stacked PR #283 snapshot `5706909b361ecafb37d426e25d8ac15064dc1bca`: **136 Python files in `src/mdscope/`**, approximately **1.30 MB**. All paths are enumerated in the appendix. Size is a *triage signal*, not a cohesion measure.
- **Actual static graph:** `mdscope-architecture-map` from #282's Python 3.11 GitHub Actions run #45060063221, recorded source SHA `9119fee4d63d46f4c2ae129a7b3db1d7e81ad1c0`: 134 modules, 425 potential directed module-import edges, 665 indexed functions/methods and 1 multi-module strongly connected component: `reporting.figure_inventory` ↔ `reporting.table_reports`, both edges function-local/deferred. This is **not evidence of an import-time cycle**. This graph predates #283's final SDC4 relocation: regenerate before quoting *current-head* edge counts.
- **Inspected code paths:** `contracts/outputs.py`, `registry.py`, `config/{__init__,runs,validators}`, `config_schema.py`, `planning.py`, `pipeline/{analysis_selection,runner,streaming,fused_all,persisted_inputs,replica_tasks}`, `studio/{app,model}`, solvation/membrane/SDC4 observers, statistics selection/state models, and reporting. The full inventory below is structural; this is **not a line-by-line scientific audit of every file**.
- **Baseline limitation:** no measured refactor speedup, import timing, new-analysis edit-count experiment, full current-head graph, or scientific real-data parity is asserted by this review. Use the existing `architecture_map.py`, synthetic fixtures, and local real-data checks when implementing behavior changes.

## Module and package assessment

| Area | Modules | Current role and quality | Proposed action |
|---|---:|---|---|
| Root: `registry.py`, `planning.py`, `config_schema.py`, `cli.py`, others | 8 | `registry.py` (~412 lines) correctly names outputs and companion paths but gets pure selection predicates from `pipeline/`. `planning.py` (~989 lines) combines output plan, saved-data inspection, metadata, figures and authoring. Two valid but broad read-only responsibilities, not an excuse for wholesale splitting. | **High priority:** remove inverted selection dependency and metadata duplication; separate read-only planning concerns only at demonstrated change/lifecycle boundaries. |
| `contracts/` | 5 | `outputs.py` owns primary/companion/schema facts. Its descriptor does *not* fully own figure/stats eligibility; `artifact_schema()` still has analysis branches. `state_feature_spec.py` imports `output_paths.py`, which itself reexports contracts—an unnecessary upward contract dependency. | Keep one data-only contract owner; move saved-source static identity down, and make missing capability facts explicit. Do **not** create a plugin class per analysis. |
| `config/` | 15 | YAML resolution and per-analysis validation now separate. Some tiny validators and contract wrappers are costly relative to their logic, but domain grouping is substantially clearer than the prior monolith. | Consolidate micro-validators when they change together; retain established validation order/messages and one public config API. Give analysis enablement one config-owned home. |
| `core/` | 17 | Geometry, PBC, glycan/sequence identity, persistence and checkpointing are meaningful independent contracts. `core/alignment.py` imports `config.analysis_groups`, reversing the intended foundation direction. | Keep scientific kernels separate; pass resolved groups or move config-specific selection adaptation above core. Do not merge unrelated kernels because they are short. |
| `observables/` | 21 | Generic observables coexist with fused observers and legitimate saved-data postprocessors. `solvation.py` (~837 lines) and `solvation_postprocess.py` (~537) have different lifecycle roles. `multi_origin_pipeline.py` is a standalone trajectory-execution adapter living among observables. | Keep observer and postprocess boundaries when responsibilities differ. Consider moving *only the standalone runner* into pipeline during #279, without moving the MSD numerical kernels. |
| `interactions/` | 6 | Cohesive detection and persistence identity functions; clear separate scientific estimands (H bonds versus salt bridges versus lipid/ligand contacts). | **Keep.** Only merge truly duplicate infrastructure, not chemically different kernels. |
| `mechanisms/` | 6 | Generic local environment, motif residence and occlusion are appropriately separated. `state_outcomes.py` is an SDC4-specific downstream consumer in a generic-named package. | Inspect call sites before relocating SDC4-specific outcome presentation; no generic base class or forced move. |
| `systems/` | 4 | #283 gives SDC4-specific geometry/potassium audit a clear home. `systems/sdc4/mechanisms.py` is ~1,366 lines but largely one large stateful observer; `enabled` currently comes from `pipeline.analysis_selection`. | Keep SDC4 layer; remove its dependency on pipeline. Refactor the observer only when independent preparation, frame update, flush, or publication seams can be proven and tested. |
| `pipeline/` | 11 | `runner.run_replicate` (~210 lines) owns many pending-output/fallback decisions; `fused_all.run_fused_only` already replays compatible observers in one source-frame pass; `streaming` owns delicate source-before-fit and checkpoint state. `persisted_inputs` has distinct validity semantics. | **Highest-value implementation:** issue #279's one derived execution decision, shared compatible replay and explicit reusable/incomplete/invalid states. Avoid a second analysis registry. |
| `statistics/` | 22 | Domain estimators are already separately testable; `common.py` centralizes inference primitives and `dispatch.py` handles finalization. Module granularity mostly follows estimands. | **Keep.** Preserve simulation-replica sampling units; don't combine estimators or introduce a generic estimator framework to reduce file count. |
| `reporting/` | 15 | `plotting.py` (~1,006 lines) covers time series, profiles and interactions while specialist figure modules already exist. PDF/inventory modules show deferred bidirectional imports. | Examine data-loading/rendering/publication seams and mutual function-local imports; split by recurring figure changes if it reduces edit surface, not merely because of LOC. |
| `studio/` | 6 | `app.py` (~1,007 lines) renders UI; `model.py` (~677) owns draft lifecycle and edit rules. The separation is useful. Preflight is reused via pipeline rather than reinvented. | **Keep** the app/model distinction. Improve read-only service calls and caching only where tested; do not build a second GUI schema/catalog. |

### Highest-value evidence-backed seams

1. **Activation/configuration:** `registry.fused_output_plan()` imports `pipeline.analysis_selection`; `systems/sdc4/mechanisms.py` imports the same module. Move *pure* enabled-family/run selection under existing `config/runs.py` or equivalent lower config owner; migrate consumers and remove the old internal path, preserving selection order.
2. **One authoritative output identity:** `contracts/outputs.py` owns static outputs, yet `registry.required_artifacts` has conditional branches and `planning.ANALYSIS_METADATA` / `analysis_figure_candidates` own additional capabilities. Do **not** simply stuff all statistics and figure execution into a descriptor. Instead define which facts are intrinsic/pure and which are computed by distinct production or rendering owners; tests must verify that their views agree.
3. **Pending/replay decisions:** `runner.run_replicate` repeatedly checks files, chooses fallback calls, writes markers, and triggers downstream summary invalidation. `persisted_inputs`, registry, and `RunState` must agree on the distinction between absent, verified reusable, recoverable incomplete, and invalid existing products. #279 is the implementation owner.
4. **Foundation direction:** `core.alignment` imports config for `analysis_groups`; make it a selection/config adapter above core or pass resolved group specs. `contracts.state_feature_spec` must not depend on root `output_paths`.
5. **Specialized execution and reporting:** `observables.multi_origin_pipeline` opens trajectories and publishes results; distinguish this orchestration from the FFT/MSD kernels. Reporting's `figure_inventory`/`table_reports` deferred two-way imports deserve a focused ownership review, not an urgent claim of a runtime cycle.
6. **Cohesion candidates, not automatic edits:** tiny `observables/rmsd.py` and `radius_of_gyration.py`, `statistics/families.py`, and several `config/validators` files should be retained if they establish real reuse, lightweight import, or scientific semantic boundaries. Merge only after checking callers and the likely future reason to change.

## Python package design criteria

**Hard gates: any design failing one is rejected regardless of weighted score.**

- Preserve CLI/YAML keys, defaults, resolved config hashes, `mdscope plan --format json`, persistence schemas/paths and documented public APIs. Document intentional breaking changes separately; unsupported internal paths need no permanent adapters.
- Preserve PBC and coordinate frames (especially source membrane frames before fitted protein frames), atom/residue/glycan identity, units, floating precision, temporal order, replica independence, bootstrap semantics and missing-versus-zero meaning.
- Preserve exact restart/checkpoint safety, atomic publication, failure-before-mutation, no-op zero reread, and existing fused trajectory-pass opportunities. No additional full XTC traversal without evidence and explicit acceptance.
- Maintain import-light planning and CLI entry points; heavy MDAnalysis/plotting/estimators load at use sites as justified.
- Validate a built wheel/sdist outside the checkout and Python 3.11–3.13. Keep conventional `src/` packaging; do not use a new plugin system absent demonstrated external extension needs.

**Weighted architecture score (1 poor to 5 excellent; revised 2026-11-22):** Migration complexity is no longer scored as an architecture-quality dimension. Its previous 10 percentage points are split equally between cohesion and runtime/resource efficiency. **Implementation safety and reviewability still govern PR merging**; excluding migration effort from a design score does not waive those gates.

| Dimension | Weight | Operational measure |
|---|---:|---|
| Analysis-change locality | 25% | Named-run and new-observable change-map: independent owners touched, duplicated switches/catalog entries, tests required |
| Responsibility ownership and dependency direction | 20% | Unique decision owner, AST module graph, inward imports, SCC reasons, no duplicated persisted/selection facts |
| Scientific, restart and compatibility safety | 20% | Independent oracle and lifecycle tests, snapshot parity, negative-path and interrupted-publication checks |
| Cohesion and navigation | **20%** | Functions grouped by reason to change, number of forwarding hops/indirect imports, discoverable ownership, meaningful seams |
| Runtime and resource behavior | **15%** | Real XTC opens/passes, no-op reads, wall/CPU time, peak RSS and disk write amplification; no speed claims without measurement |

These are **reasoned design estimates, not measured benchmarks**; dimension ratings are held fixed to isolate the effect of reweighting. In particular, runtime scores are hypotheses until measured on representative workloads.

| Candidate | Locality | Ownership | Safety | Cohesion | Runtime | Previous / 5 | Revised / 5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| **A. Keep existing type-based packages; clarify ownership** | 2.5 | 3.5 | 4.5 | 3.5 | 4.0 | 3.65 | **3.53** |
| **B. Enforce entire ADR-0001 `analyses/engine/interfaces` layout** | 3.5 | 3.5 | 4.0 | 3.5 | 4.0 | 3.45 | **3.68** |
| **C. Hybrid: fix ownership/change locality; regroup domains only where proven** | 4.5 | 4.5 | 4.5 | 4.0 | 4.5 | 4.33 | **4.40** |
| **D. Dynamic plugin/registry framework with per-analysis classes** | 4.5 | 3.0 | 3.0 | 2.0 | 3.0 | 3.08 | **3.18** |

**Decision unchanged: C (hybrid).** A is the low-risk fallback. B rises above A when its migration burden is excluded, but directory relocation alone has no demonstrated runtime advantage. B remains available for a *specific* domain when measured feature-change navigation justifies it. D is premature: MDscope currently ships analyses within one distribution, not a proven third-party plugin ecosystem.

**Disposition of in-flight PRs:** retain [#282](https://github.com/MarionQL/MDscope/pull/282) (real config/persistence ownership improvement, review micro-module cohesion before merging) and [#283](https://github.com/MarionQL/MDscope/pull/283) (legitimate SDC4-versus-generic-hydration boundary). Neither requires rollback merely because mandatory package migration was rejected. Merge #282 before #283 only when final-head checks and applicable contracts pass; do not conflate structural suitability with verified runtime or scientific parity.

## Preferred organization, without gratuitous moves

```text
mdscope/
  contracts/       # static data-only persisted identity/schema/option contracts
  config/          # parse, defaults, selectors, per-analysis validity
  core/            # geometry, source-frame contracts, identity, persistence/checkpoints
  observables/     # generic scientific kernels and cohesive streaming observers
  interactions/    # chemically distinct contact detectors
  mechanisms/      # generic mechanistic derived measurements
  systems/sdc4/    # actual SDC4-specific models and audits (not engine imports)
  pipeline/        # one work-plan/dispatcher/replay, validated input/recovery, streaming
  statistics/      # saved-data estimators and replica-aware inference
  reporting/       # rendering, table publication and findings
  studio/          # UI and detached draft model; reuse canonical plans
  cli.py           # thin user entry point
  registry.py      # read-only resolved outputs and producer contracts
  planning.py      # read-only inspection/catalogue; split only with real cohesion evidence
```

This is a **dependency DAG, not a strict horizontal ordering**: `contracts/` and numerical core are foundations; configuration-specific adaptation can depend on core but core must not depend on configuration; observers depend on foundations/config; system-specific science depends on generic analyses, not pipeline; pipeline depends on analyses; statistics consume saved contracts, not pipeline; reporting consumes validated saved results; CLI/Studio can call *lazy* service entry points. Document and test exceptions rather than relying on a misleading single linear order.

### What adding analyses should look like

For a **new source-frame scalar observable**: (1) define YAML/defaults and optional validation (only when it truly has parameters); (2) add the scientific kernel/observer near related science; (3) declare persisted identities and required companions in the existing output contracts; (4) select it through one runtime planning/dispatch owner; (5) let planning and Studio derive eligibility from the real registration. **Statistics/figures are explicit optional extensions**, not mandatory fake implementations. The number of legitimate responsibility owners touched may be greater than one; **no metadata fact or scheduling decision should require an independent duplicate edit**.

For a **named-run profile with restart**: repeat that exercise with two differently configured instances and a disabled instance. Assert exact run IDs, only-the-missing-run replay, output-path/column parity, failure on corrupt committed companions, and conditional figure/stat eligibility.

Before changing layout, run *both* onboarding exercises against the baseline and prospective branch; record changed modules, duplicate sources of truth, the number of manual registrations, owner hops and test failures caught. A nominal reduction in file count does not win if the edit surface and failure risk increase.

## Specific refactor candidates, graded

Benefit = potential improvement in change locality/ownership (5 highest); cost = migration effort (5 most expensive); risk = likelihood/severity of unintentional behavior changes (5 highest); confidence = strength of source evidence that this is worth changing (5 strongest). These are **review estimates**, not runtime measurements.

| Candidate | Benefit | Cost | Risk | Confidence | Recommendation |
|---|---:|---:|---:|---:|---|
| Move pure analysis-activation predicates from pipeline to existing config ownership | 5 | 2 | 2 | 5 | **Do first** after baseline; simultaneously removes registry and SDC4 reverse dependencies |
| Remove `core.alignment → config` and `contracts.state_feature_spec → output_paths` | 4 | 2 | 2 | 5 | **Do early**, preserve exact group validation and saved-source identity |
| Derive one execution/replay plan and remove duplicated runner fallback decisions (#279) | 5 | 5 | 4 | 5 | **High payoff, high risk**; staged ownership and restart tests are mandatory |
| Reduce duplicate metadata across descriptors, registry, planner and Studio | 5 | 4 | 3 | 4 | **High payoff** after execution decision baseline; prove onboarding gains |
| Move standalone multi-origin runner out of scientific observables | 3 | 2 | 3 | 4 | **Conditional** on #279 work; keep FFT/MSD kernels unchanged |
| Break SDC4's large streaming class into many components | 2 | 4 | 5 | 2 | **Do not do by file size**; extract only a demonstrated independent lifecycle seam |
| Split all files into `analyses/`, `engine/` and `interfaces/` | 2 | 5 | 4 | 2 | **Defer / not an objective**; may be justified for one domain later |
| Introduce plugin APIs / base classes for every analysis | 1 | 5 | 5 | 1 | **Reject for now**, no demonstrated third-party plugin requirement |

These candidate ratings complement the architecture-alternative weighted score above; do not add them as if they were measured improvements.

## Staged work and measurable gates

**P0 — Baseline and tests before implementation** (#269, coordinate with #279). Regenerate current HEAD graph with deferred/module-load flags and source hash. Record module inventory, current output/analysis selection registrations, sample extension change-maps, import/CLI startup cost, XTC pass/read counts, and before/after YAML/plan/config hashes. Establish feature-onboarding exercises and an actual *single source of truth* audit. Do not code a new framework to support the audit.

**P1 — Correct actual dependency inversions first** (#269, after #282/#283 settle). Move enabled selection lower without duplicating the selection list; remove `core.alignment → config` and `contracts.state_feature_spec → output_paths`; check system-science `→ pipeline`; remove superseded paths in the same PR. Preserve import-light planning.

**P2 — One execution/replay decision** (#279). Implement one pure derived per-replica work plan using existing output identities, then eliminate duplicate fallback decisions; consolidate compatible observer-only replay with existing `fused_all` rather than a second registry. Retain scientific stream and checkpoint-specific behavior. Validate fresh/no-op/mixed named-run/restart/corruption scenarios.

**P3 — Reduce analysis onboarding friction** (#269 after P2). Reconcile static output descriptors, conditional companion ownership, analysis activation, figure eligibility and statistics selectors. Explicitly document where consumer-specific logic legitimately lives. Add tests that a synthetic new analysis cannot be advertised as runnable without an execution registration and that an enabled persisted product is not omitted from planning. No generic plugin API required.

**P4 — Targeted domain consolidation** (only after P0–P3 evidence). If solvation, membrane, or structural features still require many unrelated edits, move one cohesive scientific domain at a time. Do not relocate all existing modules to `analyses/` as a standalone milestone. Assess the giant SDC4 observer and reporting/Studio files by reasons to change and test seams, not line counts.

**P5 — Final objective review.** Compare onboarding change-map, no-op runtime, pass counts, dependency violations, SCC explanations, test/production parity, packaging and clean install. Close an issue with an explicit "no structural move needed" outcome when the design already meets the criteria.

### Anti-bloat / API decision test

Before adding a new module, ask: *who owns this? which independent change would cause only this module to change? does it remove real duplication or a forbidden dependency? would a teammate find the feature faster? will it force adapters?* If answers are not specific, keep the function in its cohesive owner. A short module is fine when it protects a true public import, performance/optional dependency, or scientific contract.

Use the [PyPA src-layout guidance](https://packaging.python.org/en/latest/discussions/src-layout-vs-flat-layout/) and [pytest good integration practices](https://docs.pytest.org/en/stable/explanation/goodpractices.html) for packaging/test import isolation. Plugin discovery/entry points exist but are **not** an MDscope requirement absent actual independently distributed plugins.

## Appendix — full current Python module inventory

The inventory below covers every Python file in the stacked #283 snapshot; it is **inventory, not a claim that each individual module received a line-by-line review**.

- **(root) (8)**: `__init__.py`, `biophysical_common.py`, `cli.py`, `config_schema.py`, `errors.py`, `output_paths.py`, `planning.py`, `registry.py`
- **config (15)**: `config/__init__.py`, `config/defaults.py`, `config/runs.py`, `config/validation.py`, `config/validators/__init__.py`, `config/validators/common.py`, `config/validators/hydration.py`, `config/validators/inputs.py`, `config/validators/membrane.py`, `config/validators/plots.py`, `config/validators/residue_response.py`, `config/validators/sdc4.py`, `config/validators/solvation.py`, `config/validators/statistics.py`, `config/validators/structure.py`
- **contracts (5)**: `contracts/__init__.py`, `contracts/metric_selection.py`, `contracts/outputs.py`, `contracts/plot_validation.py`, `contracts/state_feature_spec.py`
- **core (17)**: `core/__init__.py`, `core/alignment.py`, `core/atom_selection.py`, `core/checkpoints.py`, `core/coordinate_contracts.py`, `core/geometry.py`, `core/glycans.py`, `core/leaflet_partition.py`, `core/lipid_geometry.py`, `core/motifs.py`, `core/persistence.py`, `core/profiling.py`, `core/progress.py`, `core/provenance.py`, `core/residence.py`, `core/sequence.py`, `core/trajectory_sources.py`
- **interactions (6)**: `interactions/__init__.py`, `interactions/hydrogen_bonds.py`, `interactions/interaction_tables.py`, `interactions/ligand_contacts.py`, `interactions/lipid_hydrogen_bonds.py`, `interactions/salt_bridges.py`
- **mechanisms (6)**: `mechanisms/__init__.py`, `mechanisms/local_environment.py`, `mechanisms/mechanistic_common.py`, `mechanisms/motif_residence.py`, `mechanisms/residue_occlusion.py`, `mechanisms/state_outcomes.py`
- **observables (21)**: `observables/__init__.py`, `observables/bilayer_normal.py`, `observables/density_1d.py`, `observables/dssp.py`, `observables/hydration.py`, `observables/lipid_contacts.py`, `observables/lipid_enrichment.py`, `observables/membrane.py`, `observables/membrane_order.py`, `observables/msd.py`, `observables/multi_origin_msd.py`, `observables/multi_origin_pipeline.py`, `observables/radius_of_gyration.py`, `observables/residue_membrane_distance.py`, `observables/rmsd.py`, `observables/rmsf.py`, `observables/sasa.py`, `observables/solvation.py`, `observables/solvation_common.py`, `observables/solvation_postprocess.py`, `observables/structure_profiles.py`
- **pipeline (11)**: `pipeline/__init__.py`, `pipeline/analysis_selection.py`, `pipeline/fused_all.py`, `pipeline/fused_distance.py`, `pipeline/fused_extensions.py`, `pipeline/fused_membrane.py`, `pipeline/persisted_inputs.py`, `pipeline/preflight.py`, `pipeline/replica_tasks.py`, `pipeline/runner.py`, `pipeline/streaming.py`
- **reporting (15)**: `reporting/__init__.py`, `reporting/figure_inventory.py`, `reporting/labels.py`, `reporting/lipid_contact_bootstrap.py`, `reporting/membrane_figures.py`, `reporting/plot_ids.py`, `reporting/plotting.py`, `reporting/report.py`, `reporting/residue_response.py`, `reporting/sdc4_condition_figures.py`, `reporting/sdc4_exports.py`, `reporting/sdc4_figures.py`, `reporting/shared_scales.py`, `reporting/solvation_figures.py`, `reporting/table_reports.py`
- **statistics (22)**: `statistics/__init__.py`, `statistics/common.py`, `statistics/contracts.py`, `statistics/coupling.py`, `statistics/coverage.py`, `statistics/diffusion.py`, `statistics/dispatch.py`, `statistics/families.py`, `statistics/interactions.py`, `statistics/lipid_contacts.py`, `statistics/manifest.py`, `statistics/permutation_resolution.py`, `statistics/profiles.py`, `statistics/replicate_contrasts.py`, `statistics/residue_regions.py`, `statistics/residue_response.py`, `statistics/selection.py`, `statistics/shared_states.py`, `statistics/solvation.py`, `statistics/state_features.py`, `statistics/time_series.py`, `statistics/windows.py`
- **studio (6)**: `studio/__init__.py`, `studio/app.py`, `studio/launch.py`, `studio/model.py`, `studio/project_editor.py`, `studio/topology_preview.py`
- **systems (4)**: `systems/__init__.py`, `systems/sdc4/__init__.py`, `systems/sdc4/mechanisms.py`, `systems/sdc4/potassium_bridge_audit.py`

**Future recheck:** `python scripts/architecture_map.py --output /tmp/mdscope-architecture`, inspect `architecture.json` for actual source SHA, import sites and SCCs; pair the graph with a developer change-map and the existing execution/restart parity tests. No arbitrary edge/LOC threshold in CI.
