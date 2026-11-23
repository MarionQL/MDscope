# Architecture

How `src/mdscope/` is organized, why, and where each decision lives. Design criteria and the rejected alternatives are recorded in [ADR-0002](adr/0002-evidence-first-architecture.md). **Internal module aliases are not retained**: CLI commands, YAML keys, `mdscope plan --format json` and persisted outputs are the public contracts.

## Data flow

```mermaid
flowchart TD
    A["YAML + TPR/XTC"] --> B["config: resolve + validate"]
    B --> C["pipeline.preflight: topology and selections"]
    C --> D["pipeline.runner: per replica"]
    D --> E["One shared source-frame pass (pipeline.streaming)"]
    E --> F["Fused observers: membrane, extensions, distances, SDC4"]
    D --> G["One selected replay for missing fused outputs (pipeline.fused_all)"]
    D --> H["Specialized tasks: DSSP/SASA profiles, potassium audit, multi-origin MSD"]
    F --> I["Per-replica raw / cache / summary / profiles / extensions"]
    G --> I
    H --> I
    I --> J["Offline statistics (saved products only)"]
    J --> K["Findings and optional figures"]
```

Trajectories are read only inside `pipeline/` and by the specialized tasks it schedules. Statistics, reporting, planning and Studio work from YAML and saved products and never open an XTC. Membrane observers consume **source coordinates before protein fitting**; never rotate or recenter the shared frames ([Scientific contracts](scientific-contracts.md)).

## Packages and dependency direction

Packages are grouped by role. Lower layers never import higher ones; `tests/test_import_layering.py` enforces this, including function-level imports.

```text
contracts/      dependency-free: output descriptors, statistics/figure family lists, option schemas
config/         YAML loading, defaults, named runs, enabled-analysis predicates, ordered validators
core/           geometry/PBC, selections, sequence and glycan identity, persistence, checkpoints
observables/    physical measurements (structure, membrane, lipid, SASA, DSSP, MSD, density, hydration)
interactions/   chemically distinct contact detectors (H bonds, salt bridges, ligand/lipid contacts)
mechanisms/     generic derived mechanisms (local environment, motif residence, occlusion, states)
systems/sdc4/   SDC4-specific geometry, HS mechanisms and potassium-bridge audit
statistics/     replica-aware inference from saved products
reporting/      findings, tables and figures from saved results
pipeline/       when and how analyses run: streaming, replay, tasks, restart, finalization
cli.py, studio/ user entry points; planning.py, registry.py, config_schema.py are read-only views
```

Allowed direction: `contracts` ← `config`/`core` ← `observables`/`interactions` ← `mechanisms` ← `systems` ← `statistics`/`reporting` ← `pipeline` ← `cli`/`studio`. `contracts/` imports nothing from `mdscope`, which keeps planning, validation and Studio import-light; `config/` may import only `core` and `contracts`. Generic sampled contact episodes live in `core/residence.py`, below both solvation and motif residence.

The diagram below is generated from the source, so it cannot drift: an arrow means the upper unit imports the lower one (dashed: only inside functions), and an import already implied by a longer path is not drawn. `tests/test_architecture_map.py` fails when it is stale; regenerate it with `python scripts/architecture_map.py --write-docs`.

<!-- BEGIN GENERATED package dependencies; regenerate with: python scripts/architecture_map.py --write-docs -->
```mermaid
flowchart TD
    cli["cli.py"]
    config_config_schema["config/ + config_schema.py"]
    contracts["contracts/"]
    core["core/"]
    errors["errors.py"]
    interactions["interactions/"]
    mechanisms["mechanisms/"]
    observables["observables/"]
    pipeline["pipeline/"]
    planning["planning.py"]
    registry["registry.py"]
    reporting["reporting/"]
    statistics["statistics/"]
    studio["studio/"]
    systems["systems/"]
    cli -.-> pipeline
    cli -.-> planning
    config_config_schema --> contracts
    config_config_schema --> core
    config_config_schema --> errors
    interactions --> core
    mechanisms --> config_config_schema
    mechanisms --> interactions
    observables --> config_config_schema
    pipeline --> registry
    pipeline --> systems
    planning --> registry
    registry -.-> reporting
    reporting --> statistics
    statistics --> mechanisms
    statistics --> observables
    studio -.-> pipeline
    studio --> planning
    systems --> mechanisms
    systems --> observables
```

- `config/ + config_schema.py` import each other and are drawn as one node.
<!-- END GENERATED package dependencies -->

Design rationale: [ADR-0002](adr/0002-evidence-first-architecture.md). Keep role-based packages, one owner per fact, system-specific science above generic kernels, and scheduling separate from measurement.

## One owner per fact

| Fact | Owner | Consumers |
|---|---|---|
| YAML defaults, named runs, group-filter and multi-run capabilities (`GROUP_ANALYSES`, `MULTI_RUN_ANALYSES`) | `config/defaults.py` | validators, option schema, streaming, planning, Studio |
| Whether an analysis or family is enabled | `config/runs.py` | registry, pipeline, `systems/sdc4` |
| Validation order and messages | `config/validation.py` dispatching to `config/validators/` | `mdscope validate`, Studio |
| Editor option schema | `config_schema.py` | Studio, unknown-key validation |
| Persisted paths, companions, column/array schemas, producer constants | `contracts/outputs.py` | producers, registry, persisted-input checks |
| Configuration-dependent outputs and required companions | `registry.py` (`analysis_output_plan`, `fused_output_plan`, `saved_output_plan`, `required_artifacts`) | runner, planning, Studio |
| Statistics families → output folders; figure families | `contracts/statistics.STATISTICS_FAMILIES`, `contracts/plot_validation.FIGURE_FAMILIES` | YAML validation, `statistics/dispatch.py`, planning |
| Saved statistics catalog and analysis display metadata | `contracts/statistics.py` (statistics metadata), `planning.py` (`ANALYSIS_METADATA`, `authoring_catalog`) | `mdscope plan`, Studio |
| Figure IDs and renderer destinations | `reporting/plot_ids.py` | plotting, inventory, planning |
| Validity of existing outputs (schema, time grid, provenance) | `pipeline/persisted_inputs.py` | runner, Studio, finalization |
| Atomic writes, archiving, checkpoint inventory | `core/persistence.py`, `core/checkpoints.py` | all producers (kernels never write files directly) |
| Covalent glycan chains; sequence motifs | `core/glycans.py`; `core/motifs.py` + `core/sequence.py` | HS/SDC4, hydration, solvation, K-bridge audit, Studio |

Tests pin the single-owner rules that have regressed before: every statistics and figure family is configurable and catalogued, the group-filter and named-run capability sets match the schema, and analysis metadata names only contract-owned families.

Static statistics schemas, identities, estimator capabilities and candidate filenames live together in `contracts/statistics.py`; planning and runtime selection consume that owner directly.

## Replica execution and recovery

`pipeline/runner.run_replicate` makes one registry-derived decision per replica:

1. **Preflight** validates every existing output and companion (`persisted_inputs`); an invalid committed product fails before anything is rewritten.
2. **Primary stream** (`streaming.run_streaming_pass`) runs when core per-frame products are missing. It drives core measurements and every enabled fused observer from the same frames, with chunked, checkpointed commits.
3. **Fused replay** (`fused_all.run_fused_only`) recomputes exactly the pending `(analysis, run ID)` outputs in **one** source-frame pass, whatever families they belong to. Completed siblings are never recomputed or rewritten. This is the only fused replay loop.
4. **Specialized tasks** (`replica_tasks.specialized_tasks`) run work that cannot share the stream (DSSP/SASA profiles, potassium-bridge audit). Each `AnalysisTask` declares the registered outputs it `produces` and runs only when one of them is absent, the same readiness rule as fused outputs; a missing marker alone does not force recomputation.
5. **Summary** consolidates saved chunks; `finalize_project` runs statistics from saved products.

The primary stream and fused replay close their readers; MDAnalysis keeps Universes reachable, and an open reader locks the file on Windows. Measured on real data, deleting three products from a completed project and resuming decodes 44 frames in one pass instead of 126 frames in four passes, with identical outputs ([#279](https://github.com/MarionQL/MDscope/issues/279)).

## Failure boundaries

- Invalid existing primary/companion artifacts fail before replacement; absence can schedule missing work.
- Completed checkpoints must pass inventory validation. Incomplete lipid footprints without a commit reference can be rebuilt; SDC4 residue occlusion supports companion-only completion.
- Replay checks the saved frame grid when available, stages scratch chunks separately and removes them after successful publication. Completed sibling runs remain untouched.
- Topology preflight reads no XTC; statistics and plotting read saved products.

| Boundary | Condition | What to inspect / safe next action |
|---|---|---|
| `validate`, `plan --check-topology`, `replicate` | Unreadable topology or invalid/empty MDAnalysis selection | Check the reported topology path, system and selection; correct the YAML/topology before running |
| `replicate` resume | Corrupt checkpoint or committed chunk inventory | Inspect the reported replica checkpoint and committed artifacts; archive the affected replica before any clean restart rather than deleting chunks to bypass validation |
| `statistics` / `finalize` | Missing, invalid or nonmonotonic saved replica summary/required companion | Finish or inspect the named replica and file; never interpret missing data as zero or silently overwrite earlier results |
| `plotting` | Selected statistical figure lacks completed statistics | Complete statistics first, or select a descriptive-only figure supported by the saved replicas; plotting does not infer missing CIs |
| `report STATISTICS_DIR` | Missing/invalid manifest or contrast-coverage file | Inspect the named statistics directory; complete the statistical export before regenerating findings; preserve existing findings |

## Adding an analysis

A **source-frame scalar** (shaped like `radius_of_gyration`):

1. Default and enabled flag in `config/defaults.py`; a validator in `config/validators/` only if it has parameters.
2. The kernel in the package that matches its science (`observables/`, `interactions/`, …), free of file I/O.
3. Bind it to the shared pass in `pipeline/streaming.py`.
4. Declare its saved identity, units and summary storage in `registry.py` / `contracts/outputs.py`.
5. Add display metadata to `planning.ANALYSIS_METADATA`. Generic time-series statistics and plots then pick up the new column; custom statistics or figures are optional extensions.

A **named-run profile** (shaped like `residue_membrane_distance`) additionally needs `MULTI_RUN_ANALYSES`, option schema and validator entries, a fused observer in `pipeline/fused_*.py`, a descriptor with companions in `contracts/outputs.py`, and run-scoped figure patterns. Test with two configured runs and one disabled run: exact run IDs and paths, no output for the disabled run, and replay of only the missing run after deleting one.

No step should require copying a list that another module already owns; if it does, move the fact to its owner instead.

## Files and lifecycle

| `project.output` path | Purpose | Keep? |
|---|---|---|
| `resolved_config.yaml`, `run_metadata.json` | Reproducibility and provenance | Yes |
| `.state/rep-*.json` | Completion markers and required-artifact inventory | While restarting |
| `<system>/<replicate>/raw/` | Chunked primary frame series and events | Until finalized |
| `<system>/<replicate>/cache/` | Stream checkpoint, fused chunks, optional aligned arrays | **Keep for restart** |
| `<system>/<replicate>/summary/` | Consolidated scalar series, diffusion, per-replica statistics | Scientific output |
| `<system>/<replicate>/profiles/` | SASA, DSSP, density, RMSF, membrane order, residue distances | Scientific output |
| `<system>/<replicate>/interactions/` | Events, occupancy, lifetimes, presence | Scientific output |
| `<system>/<replicate>/extensions/` | Lipid, HS, solvation, SDC4 and K⁺ measurements | Scientific output |
| `project_summary/statistics/` | Comparisons, uncertainty, coverage, manifest | Regenerable from valid replicas |
| `project_summary/shared_states/` | State model, assignments, populations | Derived output |
| `project_summary/figures/` | Presentation | Regenerable |
| `.archives/` | Previous numerical products moved aside on rerun | Audit for storage |

A file's existence or a completion marker is not proof of validity: `registry.required_artifacts` names the companions and `persisted_inputs` checks identity, schema and sampling. Relative paths written into saved tables and manifests always use `/`, on every OS.

## Statistics

- Finalization builds inference, shared states, diffusion summaries and findings in a same-filesystem staging directory. Successful publication archives the previous `project_summary`; computation failures leave it intact.
- `.state/finalization.json` tracks the latest attempt. An unfinished attempt keeps prior files recoverable but prevents readers from advertising them as current. Rerun finalization to complete publication.
- Serialize finalization/plotting writers. Directory publication uses two renames with rollback through the completion-record commit, not a cross-process atomic swap.

All inference reads saved replica products. Replicas, not frames, are the independent units. `statistics/__init__.py` imports nothing heavy; `statistics/dispatch.py` owns project inference and imports estimators only when a family runs. `statistics/manifest.py` is the read-only manifest reader: planning, plotting and findings require an explicitly `complete` manifest.

| Family | Input | Output folder |
|---|---|---|
| `time_series` | Scalar and group series (plus the optional `matched_time_series` window) | `statistics/time_series/` |
| `profiles` | RMSF, density, DSSP/SASA, residue distances, lipid order, occlusion | `statistics/profiles/` |
| `interactions` | Saved presence, events, lifetimes | `statistics/interactions/` |
| `solvation_rdfs` | RDF curves and shell integrals | `statistics/solvation/` |
| `diffusion` | Multi-origin curves refitted on common lags | `statistics/diffusion/` |
| `coupling` | Matched-grid replicate correlations | `statistics/coupling/` |
| `lipid_contacts` | Lipid contact matrices | `statistics/lipid_contacts/` |

Contrast coverage is written to `statistics/coverage/contrast_coverage.parquet`. `statistics/common.py` holds the shared inference primitives; `statistics/shared_states.py` fits one state model across replicas from saved features. Estimator assumptions: [Scientific contracts](scientific-contracts.md).

## Reporting and Studio

- `reporting.plot_ids.describe_plot` owns figure prerequisites; descriptive lipid species matrices can use raw saved products, while class, pair-count and comparison figures require completed statistics.
- `reporting/` renders saved results only. `table_reports` owns `statistics_tables.pdf` and its manifest; `record_generated_figures` alone writes `figure_provenance.json`; `figure_inventory` lists published figures and marks legacy PDFs without publication identity as unverified. Panels of the same observable share axis limits via `reporting.shared_scales`.
- Studio authors YAML and never runs trajectories, inference or production rendering. `studio.model` owns the detached draft, group/run edits and export reconciliation; `studio.app` renders widgets; figure candidates come from `planning.authoring_catalog`, topology checks from `pipeline.preflight`, and saved-output checks from `pipeline.persisted_inputs`. User guide: [Studio](../studio.md).

## Change map

| Editing | Check together | Target tests |
|---|---|---|
| Scalar/profile observable | Kernel, streaming binding, registry, summary, statistics family | `test_streaming_smoke.py`, `test_statistics_modes.py`, the numerical oracle |
| Contacts/chemistry | Detector, interaction tables, identity fields, statistics | `test_interaction_statistics_reference.py` |
| HS, lipid, solvation or SDC4 | Mechanism/observable, selections, YAML, companion plan | `test_sdc4_geometry_references.py`, `test_solvation.py`, `test_scientific_reference_oracles.py` |
| Statistical estimator | Family, common kernels, coverage | `test_publication_references.py`, `test_statistics_controls.py` |
| Restart and completeness | Checkpoints, runner, registry, persistence | `test_checkpoint_integrity.py`, `test_mixed_observer_restart.py`, `test_issue128_restart_artifact_contract.py` |
| Figures, CLI, packaging | Reporting or CLI, scripts | `test_report.py`, `test_cli_end_to_end.py`, `test_package_contract.py` |

Which checks a change needs before merge: [Check policy](check-policy.md).

## Inspecting the code graph

```bash
python scripts/architecture_map.py --output /tmp/mdscope-architecture   # modules, imports, calls, cycles
python scripts/architecture_map.py --write-docs                         # regenerate the package diagram above
```

`architecture.json` records the source revision, module imports (including deferred ones), statically resolved calls and multi-module import cycles; CI publishes it as the `mdscope-architecture-map` artifact. The current map has no import cycles. Edges are *possible* dependencies, not observed execution; do not gate CI on edge or module counts.

## Memory and integrity reuse

- Persistence streams Parquet record batches and publishes NPZ archives atomically; H-bond tabulation processes one committed chunk at a time.
- Large individual order profiles use disk-backed staging; normal profiles retain the batched fast path. Other scientific derivatives may still materialize tables.
- Bootstrap weights preserve eight-key eviction and RNG scheduling, with a 128 MiB memory budget and disk-backed overflow. Batched reductions avoid whole-matrix float temporaries.

- Preflight/final audit reuse successful integrity scans only within one invocation while device/inode, size, mtime and ctime match. Changed files are rescanned; schema/time/missing-product checks remain mandatory. No cache persists across runs.
