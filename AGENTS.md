# MDscope agent guide

Read the guides relevant to the change:

| Task | Source of truth |
|---|---|
| Scope, modularity, GitHub Flow, reviews | [Contributing](docs/developer/contributing.md) |
| Ownership, execution, recovery | [Architecture](docs/developer/architecture.md) |
| Architecture decisions and compatibility gates | [ADR-0002](docs/developer/adr/0002-evidence-first-architecture.md) |
| Coordinates, PBC, units, identities, estimators | [Scientific contracts](docs/developer/scientific-contracts.md) |
| Known scientific limits and verification gaps | [Scientific audit](docs/developer/scientific-audit.md) |
| Checks, commands and evidence reuse | [Check policy](docs/developer/check-policy.md) |
| SDC4 numerical references | [SDC4 fixtures](docs/developer/sdc4-fixtures.md) |
| Generalization with SDC4 compatibility | [YAML contract](docs/configuration.md#yaml-contract-for-generalization) |
| User settings and execution | [Configuration](docs/configuration.md) · [Analyses](docs/analyses.md) · [Workflows](docs/workflows.md) |

## Delivery

- Complete the authorized milestone across implementation, integration, tests and issue tracking. Small commits are review units, not stopping points.
- Stop for a real blocker or user decision, not an arbitrary patch size or number of passing tests.
- Use a focused branch and PR; do not commit directly to `main`.
- Report executed evidence, failures, skips and production limits. Numerical agreement does not establish physical interpretation or convergence.
