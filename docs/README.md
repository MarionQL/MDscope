# MDscope documentation

Start with the [repository README](../README.md) for installation and CLI usage. The guides below describe methods, configuration and scientific contracts.

## Architecture and development

- [Architecture](ARCHITECTURE.md) — current code ownership and validation gates
- [Repository rules](REPOSITORY_RULES.md) — refactoring requirements
- [Registry](ANALYSIS_REGISTRY.md) — analysis configuration and registration
- [Coordinate contracts](COORDINATE_CONTRACTS.md) — PBC, alignment and frame semantics
- [Trajectory passes](TRAJECTORY_PASSES.md) — fused execution and checkpoints

## Scientific analyses

- [Neuronal extensions](NEURONAL_EXTENSION.md) — optional membrane and interface analyses
- [SDC4 mechanisms](SDC4_MECHANISMS.md) — topology and mechanistic observer settings
- [Issue analyses](ISSUE_ANALYSES.md) — opt-in observable definitions
- [Report workflow](REPORT.md) — end-to-end SDC4 example
- [Multi-origin MSD](MULTI_ORIGIN_MSD.md) — diffusion methodology

## Historical numerical audits

These records document earlier implementation choices; they are **not** authoritative module-location guides. Consult the current architecture and source code before following a historical path.

- [Geometry audit](GEOMETRY_AUDIT.md)
- [Geometry consolidation audit](GEOMETRY_CONSOLIDATION_AUDIT.md)
- [Issue 31 numerical primitive audit](issue31-numerical-primitive-audit.md)
