# Architecture

MDscope uses a `src/mdscope/` package. [Package map](../../src/mdscope/README.md) lists ownership; [Contributing](contributing.md) defines refactoring rules. **Internal module aliases are not retained**; CLI/YAML/output compatibility is a separate contract.

## Execution flow

```mermaid
flowchart TD
    A["YAML + TPR/XTC"] --> B["CLI / config / topology preflight"]
    B --> C["pipeline.runner"]
    C --> D["Shared source-frame pass"]
    D --> E["Fused membrane, extensions, distances, SDC4 + core observables"]
    C --> F["Optional additional passes: multi-origin MSD, K audit, DSSP variants"]
    E --> G["Per-replica raw/cache/profiles/extensions"]
    F --> G
    G --> H["Replicate summary + registered artifact validation"]
    H --> I["Offline statistics + shared states + coverage"]
    I --> J["Findings and optional figures"]
```

**Not every configuration is single-pass.** Membrane observables consume original source coordinates before protein fitting; optional `msd_trajectory` is synchronized separately.

| Command | Effect |
|---|---|
| `mdscope init --output X.yaml` | Create starter configuration |
| `mdscope validate X.yaml` | Check configuration and topology selections; no numerical results |
| `mdscope replicate X.yaml --index N` | One indexed replica; no project statistics |
| `mdscope run X.yaml` | All replicas and numerical finalization |
| `mdscope statistics X.yaml` / `mdscope finalize X.yaml` | Recompute project statistics from saved replica products |
| `mdscope plotting X.yaml` | Optional saved-data figures |
| `mdscope report STATISTICS_DIR` | Regenerate findings from a completed statistics export |

Key interfaces: `cli.py`, `config.py` (`analysis_runs`, `iter_replicates`), `registry.py` (`analysis_output_plan`, `required_artifacts`), `pipeline/runner.py`, `pipeline/streaming.py`.

## Inspect the current module and call graph

Run from a checkout with Python 3.11+ (and optionally an installed `mdscope`
distribution):

```bash
python scripts/architecture_map.py --output /tmp/mdscope-architecture
# Optional: Graphviz turns DOT graphs into browsable SVGs:
dot -Tsvg /tmp/mdscope-architecture/modules.dot -o /tmp/mdscope-architecture/modules.svg
dot -Tsvg /tmp/mdscope-architecture/calls.dot -o /tmp/mdscope-architecture/calls.svg
```

`architecture.json` contains sorted module dependencies, statically resolved
function calls, multi-module strongly connected components (SCCs), source
revision, Python version, platform and installed distribution version.
The corresponding `modules.dot` and `calls.dot` are generated from
**the current `src/mdscope/` tree**. The CI checks also publish these
files as an `mdscope-architecture-map` artifact. Do not review the old
`analysis/` pyan3 graph as though it represented the current package.

**Static limitations:** imports within conditionals/functions count as
*possible* dependencies, not observed execution; calls are only recognized
when they resolve to a known function or direct method. Decorators, callbacks,
dynamic imports, aliased attributes and reexports can be missed or misleading.
SCCs flag places to inspect; they do **not** prove a runtime circular-import
failure. No trajectory, checkpoint or scientific calculation is executed.

### Static import-cycle review

The first current-source CI inventory for #180 identified two two-module SCCs.
Both involve **function-local/deferred** imports, not a demonstrated import-time
cycle. The generated `architecture.json.import_sites` includes the file line
and `module_load` versus `deferred` classification so reviewers can confirm
the context before changing code.

| Static cycle | Why it appears | Decision |
|---|---|---|
| `observables.lipid_enrichment` ↔ `pipeline.fused_membrane` | The optional `run_lipid_enrichment` wrapper imports the fused runner *inside the function*; the runner imports `LipidFootprint` at module load. | No eager circular import; investigate wrapper ownership only if an actual consumer/API refactor needs it. |
| `observables.solvation` ↔ `observables.solvation_postprocess` | Solvation finalization lazily imports postprocessors; postprocessors lazily import `rdf_from_counts` and `EPISODE_COLUMNS` from the solvation implementation. | No eager circular import; preserve numerical and episode contracts until a concrete reason to extract shared definitions. |

In this inventory, `pipeline.runner` has many outgoing imports because it
coordinates stages; `config.py` is imported widely because it owns resolved
settings. Neither degree count alone justifies splitting a module. Do not add
CI thresholds for edge counts or SCC counts.

| Responsibility | Canonical owner | Avoid |
|---|---|---|
| YAML defaults, schema, user selection configuration | `config.py`, `config_schema.py` | Another Studio-only configuration engine |
| Topology/selection checks | `core/validation.py`, `core/atom_selection.py` | GUI parsing MDAnalysis expressions independently |
| Atomic writes, archiving, checkpoints | `core/persistence.py`, `core/checkpoints.py` | New output writers inside scientific kernels |
| Numeric observables and interaction identity | `observables/`, `interactions/interaction_tables.py` | Re-implementing measurements in the planner |
| Streaming, scheduling, restart and project finalization | `pipeline/streaming.py`, `pipeline/runner.py` | An extra trajectory pass for each observer |
| Replica-aware inference, resampling, comparisons | `statistics/` | Treating frames as independent replicates |
| Figure generation and provenance | `reporting/` | Synthetic-preview-only rendering code |
| Run/output introspection | `planning.py` over `registry.py` | Duplicate expected-output filename lists |

### Follow an RMSD group from YAML to figure

1. In YAML, define a named `selections.groups.acidic_patch` containing an
   **MDAnalysis** selection such as `segid PROA and resid 35:62`, with optional
   independent `alignment_selection`. Enable `analyses.rmsd` and choose the
   group using its ID. This example requires a matching user topology.
2. `mdscope plan analysis.yaml --format json --check-topology` calls the
   existing configuration/planning and topology preflight; no XTC is read.
3. `mdscope replicate analysis.yaml --index 0` validates the topology before
   `pipeline/streaming.py::run_streaming_pass`. Membrane observers see
   **source coordinates before fitting**. RMSD uses fitted features when
   alignment is enabled (otherwise source positions); individual groups may
   have separate fitting atoms. Sharing this loop avoids redundant XTC reads.
4. The streaming pass commits chunked raw data. `pipeline/runner.py` merges
   these into the replica's `summary/group_timeseries.parquet`, including
   `group`, `time_ps` and `rmsd_A` among supported columns. Atomic writes
   and checkpoint inventory support safe resume.
5. Offline `statistics/time_series.py` can estimate replica-aware summaries
   and condition comparisons from saved products. It does **not** revisit XTC.
   `reporting/plotting.py` can also render a selected descriptive RMSD
   time-series PNG directly from saved replica data; inferential plots require
   their corresponding completed statistical products.

**Import cleanup rule:** `core.persistence` owns filesystem utilities, but
`config.py` currently reexports some of them. Keep compatibility until all
consumers have been audited; prefer direct canonical imports in *new* code.
`interactions/__init__.py` explicitly reexports the existing interaction
function API rather than using `import *`; new consumers should import from
the owning submodule directly. Neither boundary justifies a second implementation.

## Files and lifecycle

| `project.output` path | Purpose | Keep? |
|---|---|---|
| `resolved_config.yaml`, `run_metadata.json` | Reproducibility/provenance | Yes |
| `.state/rep-*.json` | Completion and required artifact inventory | Yes while restarting |
| `<system>/<replicate>/raw/` | Chunked primary frame series/events | Until safely finalized/reusable |
| `<system>/<replicate>/cache/` | Stream checkpoints, fused chunks, optional aligned arrays | **Keep for restart** |
| `<system>/<replicate>/summary/` | Consolidated scalar series, diffusion, per-replica statistics | Scientific output |
| `<system>/<replicate>/profiles/` | SASA, DSSP, density, RMSF, membrane order, residue distances | Scientific output |
| `<system>/<replicate>/interactions/` | Events, occupancy, lifetime, presence | Scientific output |
| `<system>/<replicate>/extensions/` | Lipid/HS/solvation/SDC4/K measurements | Scientific output |
| `project_summary/statistics/` | Comparisons, uncertainty, coverage, complete manifest | Regenerable from valid replica products |
| `project_summary/shared_states/` | Model, assignments and populations | Derived scientific output |
| `project_summary/figures/` | Diagnostic presentation | Regenerable |
| `.archives/` | Previous numerical products moved on rerun | Audit periodically for storage cost |

**Completeness:** file existence or a completion marker alone is insufficient. `registry.required_artifacts` defines companions; runner/checkpoint validation checks identity, structure, sampling and safe replay. Avoid hash-based gating of compatible historical outputs. Do not delete restart chunks or replace observable formats as part of unrelated cleanup.

## Statistics inputs and outputs

All inference reads **saved** replica products. The manifest records methods, enabled families and completion; source-file snapshots detect input changes during finalization.

| Family | Input | Main output folder |
|---|---|---|
| `time_series` | Scalar/group series; selected solvation/lipid/bridge metrics | `statistics/time_series/` |
| `profiles` | RMSF, density, DSSP/SASA, residue distances, lipid order, occlusion | `statistics/profiles/` |
| `interactions` | Saved presence/events/lifetimes | `statistics/interactions/` |
| `solvation_rdfs` | RDF curves and shell integrals | `statistics/solvation/` |
| `diffusion` | Multi-origin curves refitted on common lags | `statistics/diffusion/` |
| `coupling` | Matched-grid replicate correlations | `statistics/coupling/` |
| `matched_time_series` | Optional second scalar window | `statistics/matched_time_series/` |
| `coverage` | Family-level comparisons | `statistics/coverage/contrast_coverage.parquet` |

- `statistics/common.py`: autocorrelation, block bootstrap, replicate hierarchy, paired/unpaired contrasts, permutation and FDR.
- `statistics/replicate_contrasts.py`: replica-level RDF-shell/diffusion/coupling contrasts; do **not** replace with frame-level resampling.
- `mechanisms/shared_states.py`: shared state fit/populations; not every conditioned table enters the central seven-family coverage audit.
- A 3+3 exact unpaired two-sided test has minimum p = 0.10; pointwise RDF mean/SD is **descriptive**, not a CI. See [Scientific contracts](scientific-contracts.md).

## Change map

| Editing | Check these together | Target tests |
|---|---|---|
| Scalar/profile observable | Observer, registry, summary collector, statistics family | `test_streaming_smoke.py`, `test_statistics_modes.py`, relevant numerical oracle |
| Contacts/chemical entities | Detector, interaction tables, identity fields, statistics | `test_interaction_statistics_reference.py`, contact/chemical tests |
| HS, lipid, solvation or SDC4 | Mechanism/observable, selections, YAML, companion output plan | `test_sdc4_geometry_references.py`, `test_solvation.py`, `test_scientific_reference_oracles.py` |
| Statistical estimator | Family and common kernels, coverage | `test_publication_references.py`, `test_statistics_controls.py` |
| Resume/output completeness | Checkpoints, runner, registry, persistence | `test_checkpoint_integrity.py`, `test_mixed_observer_restart.py`, `test_issue128_restart_artifact_contract.py` |
| Figures/CLI/packaging | Reporting or CLI and scripts, no unneeded science changes | `test_report.py`, `test_cli_end_to_end.py`, `test_package_contract.py` |

CI: Python 3.11/3.13, wheel/sdist from outside checkout, applicable real-trajectory smoke/physics checks. Document tests **actually run**. A synthetic or release smoke pass does not establish six-replica convergence.

## Planned: discover differences without reprocessing XTC

**Proposal, not implemented:** a read-only catalog layered over existing Parquet/NPZ products.

| View | One row per |
|---|---|
| `observable_catalog` | Identity/selection/units, source, availability, supported estimator |
| `estimates` | Replica or system estimate with status, CI and provenance |
| `comparisons` | Condition pair, system-2-minus-system-1 effect, CI, p/q and coverage |
| `discovery` | Physically meaningful region/family effect, support, exploratory rank |
| `coverage` | Every enabled output: tested, descriptive, excluded or unavailable **with reason** |

Approach:
1. Inventory configured and saved observables **before** expensive CI work; preserve species, chain, residue *topology index*, units, leaflet and normalization.
2. Screen candidate differences and replica consistency cheaply; group ordered profiles into meaningful regions.
3. Run expensive per-coordinate inference on supported/preregistered priorities; mark others `not_ci_evaluated`, **not** nonsignificant.
4. Treat selected-on-same-data intervals as exploratory unless independently confirmed or appropriately calibrated. Six independent simulations do not support high-capacity supervised feature discovery.

**No database or rewrite of existing scientific files.** Different observables (RDF shape, occupancy, residence, SASA burial, diffusion, state populations) retain different estimands.

[Workflows](../workflows.md) · [Configuration](../configuration.md) · [Validation gaps](validation-gaps.md)
