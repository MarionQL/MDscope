# Pipeline graphs

Three different questions need separate diagrams: **execution** (where
frames pass), **recovery** (when saved products are reused or replayed),
and **imports** (Python source dependencies). An import edge does not mean
an analysis opened a trajectory.

Regenerate all diagrams from a checkout:

```bash
python scripts/pipeline_graphs.py --output /tmp/mdscope-graphs
# Graphviz's dot executable is optional for SVG rendering.
# Open execution.svg, recovery.svg and imports_focused.svg (or their DOT files).
```

The script regenerates the existing AST inventory using
`scripts/architecture_map.py`; it also writes `recovery.dot`, `imports_full.dot`, and all four SVGs when
Graphviz is installed. Execution and recovery diagrams are **code-reviewed
conceptual flow maps**, not instrumented traces. Its
conditions must be reviewed when `runner.py`, `streaming.py` or
`fused_all.py` scheduling changes.

## Execution: where trajectories are read

```mermaid
flowchart TD
    CLI["run / replicate CLI"] --> RUN["runner.run_replicate"]
    RUN --> PREF["Validate existing products"]
    PREF --> DEC{"Core streaming or aligned cache required?"}
    DEC -->|Yes| STR["streaming.run_streaming_pass<br/>Primary shared frame loop"]
    DEC -->|No, pending instances| FO["fused_all.run_fused_only<br/>Selected observers and named runs, one frame loop"]
    DEC -->|No work needed| AUD["Per-family output audit"]
    STR --> CORE["Core measurements + checkpoints"]
    STR --> MEM["FusedMembrane"]
    STR --> EXT["FusedExtensions"]
    STR --> DIST["FusedDistances"]
    STR --> MECH["SDC4Mechanisms"]
    FO --> MEM
    FO --> EXT
    FO --> DIST
    FO --> MECH
    CORE --> AUD
    MEM --> AUD
    EXT --> AUD
    DIST --> AUD
    MECH --> AUD
    AUD -.->|Missing outputs after stream| FALL["Selective fused-family fallback readers"]
    AUD --> VER["Artifact validation + legacy stage markers"]
    FALL --> VER
    VER -.->|Enabled, not fused or already satisfied| TASK["Specialized tasks / multi-origin MSD"]
    VER --> SUM["Replica summary"]
    TASK --> SUM
    SUM -.->|Pending companions| LATE["Membrane density / hydration completion"]
    SUM --> FIN["Offline statistics / findings / figures"]
    LATE --> FIN
```

**Reading this graph:** solid arrows show possible execution flow, not
simultaneous calls; dashed arrows show conditions. All fused observer families
consume the same **source-coordinate** frame before alignment in the primary
stream, and can also run under observer-only scheduling. Some standalone or
recovery operations reopen trajectories. The runner validates existing artifacts rather than assuming that files
are valid. In observer-only execution it selects **pending (analysis, run ID)**
pairs; completed named SASA and distance outputs are not recomputed. Replay
chunks are staged separately from committed streaming chunks. Per-family
fallbacks remain available after a completed streaming pass.

## Recovery: reuse, replay, or fail safely

This is **separate from execution** because output completeness and
checkpoint provenance—not just the enabled analysis set—determine whether a
trajectory must be reopened.

```mermaid
flowchart TD
    START["run_replicate"] --> PREF["Preflight: validate existing primary and required companions"]
    PREF --> BAD{"Existing artifact invalid?"}
    BAD -->|Yes| STOP["Fail without overwriting existing results"]
    BAD -->|No| CORE{"Primary streaming required?"}
    CORE -->|Yes| STATE{"Stream checkpoint state?"}
    STATE -->|Completed, valid| REUSE["Validate committed stream inventory; reuse core products"]
    STATE -->|Absent or resumable| STREAM["Resume or run primary stream after checkpoint validation"]
    STATE -->|Invalid or inconsistent| STOP
    REUSE --> FAMILY["Check fused-family output contracts"]
    STREAM --> FAMILY
    FAMILY -->|Missing primary outputs| FALL["Select only missing membrane / extension / distance / SDC4 outputs"]
    FAMILY -->|Complete| VALID["Validate registered artifacts and companions"]
    FALL --> SCRATCH2["Independent selective replay; staging isolated from stream chunks"]
    SCRATCH2 --> VALID
    CORE -->|No| PENDING{"Any pending fused output instances?"}
    PENDING -->|No| SKIP["Reuse saved fused products; do not open a trajectory"]
    PENDING -->|Yes| ONLY["One observer-only source-frame pass<br/>Only requested analysis + run IDs"]
    ONLY --> SCRATCH["Temporary chunks; publish selected results"]
    SCRATCH --> VALID
    SKIP --> VALID
    VALID --> SUMMARY{"Summary reusable and inputs unchanged?"}
    SUMMARY -->|Yes| EXTRA["Optional specialized / secondary passes"]
    SUMMARY -->|No| BUILD["Rebuild summary from saved products"]
    BUILD --> EXTRA
    EXTRA --> DONE["Validated replica; project finalization is offline"]
```

**Recovery contract:** existing invalid primary or companion artifacts fail
preflight; the pipeline must not silently replace them. Absent outputs can
be computed, and an incomplete lipid footprint lacking its commit reference
can be rebuilt. SDC4 residue occlusion has a documented companion-only
completion path. Complete streaming checkpoints must pass their chunk
inventory validation. If a non-streaming configuration has several pending
fused families, they share **one** source-frame traversal. A stream-enabled
configuration may require independent selective fallback passes for outputs
that were added after the stream finished.

Scratch replay chunks are discarded after successful publication; completed
outputs are not regenerated just because another run in the same family is
missing. This diagram is conceptual: see `runner.py`,
`fused_all.py` and each observer's validation for exact failure modes.

## Dependencies: what imports what

This compact graph shows the relevant layers; the generated DOT/SVG files above
contain **actual current import edges**, including conditional imports and
cycles, rather than this intentionally simplified overview.

```mermaid
flowchart TD
    CLI["cli"] -.-> RUN["pipeline.runner"]
    CLI -.-> PLAN["planning"]
    PLAN --> REG["registry + output_paths"]
    PLAN --> FAM["statistics_families"]
    RUN --> STREAM["pipeline.streaming"]
    RUN --> FONLY["pipeline.fused_all"]
    RUN --> TASK["pipeline.replica_tasks"]
    STREAM --> OBS["observables"]
    STREAM --> INT["interactions"]
    STREAM --> MECH["mechanisms"]
    FONLY --> OBS
    FONLY --> MECH
    TASK --> OBS
    TASK --> MECH
    OBS --> CORE["core / config"]
    MECH --> CORE
    INT --> CORE
    RUN --> STATS["statistics"]
    STATS --> CORE
    RUN --> REPORT["reporting"]
```

**Dependency legend:** solid arrows represent architectural dependencies in
this simplified diagram; dashed arrows represent deferred/command-specific
imports. In generated import graphs, a solid arrow means a module-load import
and a dashed arrow means a deferred import. A purple edge (when present)
identifies a multi-module strongly connected component. An import does **not**
prove that an analysis ran or a trajectory was opened.

## What each package owns

| Package | Responsibility | Examples |
|---|---|---|
| `core/` | Shared geometry, coordinate contracts, checkpoints, selection utilities | Minimum-image distances and frame provenance |
| `observables/` | **Measure** physical quantities, often directly from trajectory frames | RMSD/RMSF, SASA, membrane order, lipid enrichment, water coordination |
| `interactions/` | Detect and persist molecular contact/event relationships | H bonds, salt bridges, lipid contacts |
| `mechanisms/` | Generic scientific features and models, including saved-data states | Local environments, motif residence, occlusion, shared conformational states |
| `systems/sdc4/` | SDC4-specific observers and chemistry | SDC4 motif/HS geometry and potassium bridge audit |
| `pipeline/` | Decide **when/how** implementations run and resume; owns the shared frame loops | Streaming, fused observers, fallback dispatch |
| `statistics/` | Compare replicas, estimate uncertainty, and evaluate contrasts from saved products | Hierarchical/bootstrap CIs, RDF contrasts |
| `reporting/` | Present saved results, without reopening trajectories | Figures and findings |

An *observable* is not merely a post-hoc grouping of existing tables: many
perform PBC-aware coordinate calculations on each frame. Offline aggregation
and inference belong in `statistics/`, while display-only grouping belongs
in `reporting/`.

`systems/sdc4/` is intentionally separate from `pipeline/`: SDC4
and HS science is **what** is measured, whereas `pipeline/` owns **when** it
runs. In particular, `mechanisms/shared_states.py` works from saved features
and does not need a trajectory pass; `SDC4Mechanisms` is a per-frame observer.
Moving all mechanisms under `pipeline/` would conflate scheduling with
domain-specific scientific implementations.

**Refactor test:** share PBC geometry, atom selection and frame plumbing when
the numerical definitions agree. Do not collapse physically different
estimators into one generic observer just because they have similar output
tables. Prefer eliminating redundant trajectory passes and duplicate
aggregation contracts before moving package boundaries.


## Boundary diagnostics and safe recovery

These are **failure categories**, not fallback-to-zero behaviors. Errors should
identify the affected system/replica or saved file, while leaving committed
outputs unchanged.

| Boundary | Condition | What to inspect / safe next action |
|---|---|---|
| `validate`, `plan --check-topology`, `replicate` | Unreadable topology or invalid/empty MDAnalysis selection | Check the reported topology path, system and selection; correct the YAML/topology before running |
| `replicate` resume | Corrupt checkpoint or committed chunk inventory | Inspect the reported replica checkpoint and committed artifacts; archive the affected replica before any clean restart rather than deleting chunks to bypass validation |
| `statistics` / `finalize` | Missing, invalid or nonmonotonic saved replica summary/required companion | Finish or inspect the named replica and file; never interpret missing data as zero or silently overwrite earlier results |
| `plotting` | Selected statistical figure lacks completed statistics | Complete statistics first, or select a descriptive-only figure supported by the saved replicas; plotting does not infer missing CIs |
| `report STATISTICS_DIR` | Missing/invalid manifest or contrast-coverage file | Inspect the named statistics directory; complete the statistical export before regenerating findings; preserve existing findings |

The production checks live in `pipeline/preflight.py`,
`pipeline/streaming.py::validate_stream_inventory` and
`pipeline/runner.py::_validate_replica_products` and
`reporting/report.py::generate_report`. The **recovery graph**
above describes their ordering; topology-only preflight does not read an XTC,
and statistics/plotting reuse saved products. Keep this table as the single
developer reference for failure categories; Studio should use these existing
backend failures rather than inventing separate validation semantics.
