# MDscope documentation

Start with the [repository README](../README.md) for installation, basic usage, and CLI commands. For the production syndecan workflow, see [workflows](workflows.md) and the checked-in `synd4_full_mechanisms.yaml`.

## Using MDscope

- [Configuration and analysis registration](configuration.md) — YAML options, selection groups and registry behavior.
- [Scientific analyses and methods](analyses.md) — membrane, lipid, SDC4, solvation, interaction and diffusion methods; assumptions and limitations.
- [Workflows and reporting](workflows.md) — trajectory streaming, SLURM execution, finalization and report artifacts.

## Developing MDscope

- [Architecture](developer/architecture.md) — package boundaries and installation model.
- [Scientific audit and outstanding validation](developer/scientific-audit.md) — prioritized methodology findings, corrections, benchmarks and production commands.
- [Scientific contracts](developer/scientific-contracts.md) — coordinate conventions, geometry, numerical regression and validation.
- [Numerical reference project and release checklist](developer/release-validation.md) — portable YAML fixture, independent expected outputs, release evidence and citation/license status.
- [Contributing](developer/contributing.md) — refactoring rules and testing expectations.

Documentation is organized by enduring functionality, not GitHub issue numbers. Historical development discussions remain available in Git history and closed issues.

