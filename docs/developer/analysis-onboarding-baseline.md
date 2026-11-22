# Analysis onboarding baseline — issue #269 (P0)

**Comparison revision:** merged `main` `e436bc5dc8df4a51f4d24d981156e9134cb82ca9` (2026-11-22). **Method:** static source/edit-site walkthrough using *existing* analogous MDscope analyses, plus two read-only plan-contract probes added to `tests/test_registry_execution_boundary.py`. This did **not** implement a synthetic scientific kernel or open a trajectory. Numbers below are candidate code touchpoints, not measured developer time, runtime speedups or a completed feature-branch edit count.

## Reproducible import baseline

The GitHub Action architecture artifact `mdscope-architecture-map` (workflow 45062963237, artifact 14770193512, source hash `47f46051b2f388c0735eaf01e8aae1719acfacbd`) recorded **134 modules, 422 possible directed module import edges**, and **one two-module reporting cycle** (`reporting.figure_inventory ↔ reporting.table_reports`), with deferred imports. The artifact includes conditional/function-local imports; this is not a claim of a runtime cycle. Regenerate graph after the P1 changes with `python scripts/architecture_map.py --output /tmp/mdscope-architecture`; compare exact source hashes before quoting a new count.

Before P1, direct reverse imports included `registry → pipeline.analysis_selection`, `core.alignment → config`, and `contracts.state_feature_spec → output_paths`. The SDC4 science → pipeline edge had already been removed in merged PR #283. The P1 changes remove these remaining three *causes* while preserving consumers' import identity tests and scientific behavior; removing modules changes the graph's total node/edge count, so raw counts should not be equated with quality.

## Probe A: add a source-frame scalar like radius of gyration

Trace the existing `radius_of_gyration` path as a developer would for a new `shape_scalar`. It currently uses:

| Owner | Baseline edit decision | Required for a new analogous scalar? |
|---|---|---|
| `config/defaults.py` | YAML default, selection and enabled status | Yes |
| `config_schema.py` | Studio/editor option semantics for analyses supporting selection | Only if selection/options differ |
| `core/alignment.py` (P1: `config/runs.py`) | whether structural protein selection is required | Conditional on input choices |
| `observables/radius_of_gyration.py` | source-coordinate scientific operation | Yes, for a new estimator |
| `pipeline/streaming.py` | binding to shared pass and per-frame summary | Yes under current implementation |
| `registry.py` | scalar identity, units and summary storage | Yes under current implementation |
| `planning.py` | display purpose, statistics and figure-family eligibility | Yes for cataloguing |

**Observation:** six always-or-usually involved owners plus optional selection/schema choices (seven listed), even though actual numerical work is one tiny function. The generic time-series statistics and plotting engines can often reuse the new column without a new estimator implementation. There is no requirement to create a dedicated `analyses/shape_scalar/` package. **Repeated fact:** an analysis/metric identity is described in more than one of streaming, registry and planning; reducing that coordination is P3, not P1.

The read-only CI probe enables `radius_of_gyration` and checks that its saved registry product is the canonical `summary/timeseries.parquet`, without loading an XTC.

## Probe B: add a named-run profile like residue-to-membrane distance

Trace existing `residue_membrane_distance` for a new `named_surface_profile` with two configured run IDs, only one enabled:

| Owner | Baseline edit decision |
|---|---|
| `config/defaults.py` + `MULTI_RUN_ANALYSES` | default settings and named-run eligibility |
| `config_schema.py` + `config/validators/membrane.py` | supported selections, option schema and physical constraints |
| `observables/residue_membrane_distance.py` | actual numerical kernel and saved profiles |
| `pipeline/fused_distance.py` + `pipeline/runner.py` | shared pass, selective named-run replay, restart/complete checks |
| `contracts/outputs.py` | persisted primary/companion paths and array/column schema |
| `registry.py` | configured profile identities and conditional artifact inventory |
| `planning.py` | analysis metadata, preview, figure patterns, named-run inventory |
| `statistics/profiles.py` / `reporting/plotting.py` | optional aggregation/inference and rendering if generic handling is insufficient |

**Observation:** at least 10 independently owned code files for a fully integrated named profile (and more if it needs custom statistics/figures). The meaningful complexity is *multiple conditionally coordinated production decisions*, not the count of packages. Required run-specific behavior includes no output for a disabled run, canonical named-run paths, and reprocessing only a missing eligible instance. The CI probe checks the first two; #279 owns the selective-recovery implementation and tests.

## Decisions that should remain separate

YAML interpretation, scientific numerical kernels, checkpoint-safe execution, persisted artifact identity, replica-aware statistics, and figures change for different reasons. Do not force them into one mega-descriptor or one giant class merely to lower the number of touched modules.

## Duplication and next measurement gates

1. Before P1, the same enabled-family selection logic was callable through `pipeline.analysis_selection` from producer classes *and* the read-only registry. P1 moves its authoritative predicates to existing `config/runs.py`, removing the execution import from the registry; production uses the identical callable rather than copying logic.
2. `runner.run_replicate` still has four family-specific replay/fallback branches (membrane, extensions, residue-distance, SDC4), on top of the observer-only `fused_all` pass. This is **not proof they are interchangeable**: their checkpoint keys, conditional companions, and replay semantics differ. #279 should derive a single work-decision plan before changing dispatch.
3. Output schema is centralized in `contracts/outputs.py`, but `registry.required_artifacts()` resolves conditional companions, while `planning.ANALYSIS_METADATA` and `analysis_figure_candidates()` advertise optional capabilities. Avoid duplicated facts without pretending production and figure routing are the same responsibility.
4. **Runtime baseline not yet measured**: full XTC reader/pass counts, no-op I/O, wall time, RSS and import latency need representative local/ARC runs. Existing logic intends one shared source-frame streaming pass, or one selected fused-only pass, with additional fallback passes if products remain missing; this is a code-path statement, not measured efficiency.
5. **Next scorecard:** after #279/#269 P3, repeat the same two probes and report code owner touchpoints, duplicate registration sites, unexpected unavailable outputs, exact YAML/plan/config parity, and real trajectory pass counts. Do not attribute P1's layer cleanup a measured speedup.
