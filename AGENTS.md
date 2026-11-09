# MDscope agent instructions

Use the existing project documentation as the source of truth. Read the files relevant to the change **before** editing; do not duplicate their contents here.

- [Contributing and engineering rules](docs/developer/contributing.md): ownership, minimal architecture, implementation and review expectations.
- [Architecture](docs/developer/architecture.md): module boundaries, data flow and packaging.
- [Scientific contracts](docs/developer/scientific-contracts.md): coordinate/PBC conventions, estimator semantics, units and reference expectations.
- [Scientific methodology audit](docs/developer/scientific-audit.md): established limitations and unresolved scientific questions; do not turn caveats into validated claims.
- [Validation effectiveness and gaps](docs/developer/validation-gaps.md): independent oracles, mutation evidence, known coverage and production-validation limits.
- [Configuration](docs/configuration.md), [analyses](docs/analyses.md), and [workflows](docs/workflows.md): user-facing YAML, scientific options, execution, finalization and outputs.
- [SDC4 reference fixtures](docs/developer/sdc4-fixtures.md): system-specific test contracts, when touching SDC4 mechanisms.
- [Repository documentation index](docs/README.md): navigation to other maintained material.

## Local code map

Use Graphify for changes that span modules, such as following an analysis from
its runner through output validation and reporting. Query specific symbols with
`graphify explain SYMBOL` or `graphify path FROM TO`, then verify the callers,
YAML, output contracts and tests in the source. Use ordinary source search for
small edits or when Graphify is unavailable.

From the checkout root, build the optional code-only map with:

```sh
uv tool install graphifyy==0.9.104
GRAPHIFY_BACKEND=none graphify extract . --code-only
GRAPHIFY_BACKEND=none graphify cluster-only .
```

The generated JSON, report and HTML live in the ignored `graphify-out/` folder
in each checkout. When using the map, run `GRAPHIFY_BACKEND=none graphify update .`
before querying after a pull or branch change, and after edits that change
callers or module boundaries. Keep the map local to development checkouts.

## Working agreement

1. Inspect existing callers, contracts, tests and relevant docs; state the intended change and its acceptance criteria. Preserve the distinction between **numerically correct** and **physically interpretable**.
2. Prefer the smallest coherent change to existing infrastructure. Do not add compatibility wrappers, duplicated helpers, frameworks or tests without a demonstrated need. Follow the public-versus-internal compatibility policy in contributing.md.
3. Explicitly consider the six failure modes: wrong requirements/assumptions, scientific mathematics, state/data flow, integration contracts, inadequate verification, and failure/recovery handling.
4. For meaningful bug fixes, reproduce the behavior and use a regression test that would fail on the original implementation. Prefer independent analytical or alternative-implementation oracles for numerical claims. Do not weaken a failing assertion merely to make CI pass.
5. Check required output identity, completeness and provenance where the changed workflow relies on persisted files. A file's existence or a passing test is not sufficient evidence of correctness.
6. Run the smallest relevant tests first, then the applicable repository CI or real-trajectory checks. Report what was **actually executed**, what was not, and any unresolved limitations. Do not claim full production validation from short fixtures or synthetic tests.

For current module and scientific specifics, follow the linked documentation rather than expanding this file into a second policy manual.
