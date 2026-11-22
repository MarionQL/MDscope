# MDscope agent checklist

Read only the guides relevant to the change; do not duplicate them here.

| Task | Source of truth |
|---|---|
| Ownership, minimal changes, reviews | [Contributing](docs/developer/contributing.md) |
| Module boundaries, persisted flow | [Architecture](docs/developer/architecture.md) |
| Where moved or new modules belong (target layout) | [ADR-0001](docs/developer/adr/0001-target-package-architecture.md) |
| Coordinate frames, PBC, units, numerical semantics | [Scientific contracts](docs/developer/scientific-contracts.md) |
| Known scientific limitations | [Scientific audit](docs/developer/scientific-audit.md) |
| Which checks to run and when | [Check policy](docs/developer/check-policy.md) |
| Test effectiveness, remaining coverage | [Validation gaps](docs/developer/validation-gaps.md) |
| SDC4-specific reference cases | [SDC4 fixtures](docs/developer/sdc4-fixtures.md) |
| Generalizing analyses while preserving system-specific choices | [YAML generalization contract](docs/configuration.md#yaml-contract-for-generalization) |
| YAML, observable options, workflows | [Configuration](docs/configuration.md) · [Analyses](docs/analyses.md) · [Workflows](docs/workflows.md) |

## Avoid

- Dead code and stale exports; remove obsolete paths only after checking callers and compatibility.
- Swallowed errors or redundant guards; catch expected failures at boundaries with actionable context.
- Duplicate algorithms, unjustified dependencies, and speculative options.
- Redundant prose; retain concise scientific contracts and non-obvious rationale.

Use **GitHub Flow** and the issue title/acceptance conventions in [Contributing](docs/developer/contributing.md). Keep changes on a focused issue branch and PR, not directly on `main`.

## Refactoring and modularity principles

1. **Prefer cohesive modules over maximally small modules.** Group related functions with shared responsibilities; do not split files merely to reduce line counts or lint complexity.
2. **Avoid compatibility adapters for internal imports.** Migrate internal callers to canonical owners and remove obsolete forwarding modules; preserve adapters only for a demonstrated external/public compatibility contract or a temporary, tracked migration.
3. **Separate by reason to change.** Split code when ownership, dependencies, invariants, or change lifecycles differ. Keep functions together when they typically evolve together.
4. **Do not introduce abstractions without a concrete benefit.** Add registries, interfaces, base classes, wrappers, or helper layers only when they remove actual duplication, coupling, or a tested maintenance problem.
5. **Measure complexity across the whole system.** Review import depth, indirection, navigation cost, number of touched modules per feature, duplicate responsibilities, cycles, and testability—not just individual file length or Ruff scores.

Before a substantial architectural refactor, assess whether existing modules can be consolidated instead. Retain boundaries that protect scientific semantics, clear ownership, or dependency direction.

## Before merging

1. Identify intended behavior, inputs/outputs, units, populations, provenance and failure outcomes. Separate **numerical correctness** from **physical interpretation**.
2. Review six risks: assumptions, mathematical errors, state/data flow, integration, insufficient verification, and failure/recovery.
3. Complete the smallest coherent solution that satisfies the authorized milestone, across modules when needed. Do not interpret minimal changes as a reason to stop after a small patch. Remove superseded internal paths within scope and after migrating callers; use adapters or separate tests when they protect a real contract.
4. For substantive bugs, reproduce failure and test the affected entry point. Prefer independent analytical references; never weaken assertions to pass CI.
5. Validate identity, frame/time alignment, PBC, missing-vs-zero and persisted-product completeness. Existing output files are not automatically valid.
6. Follow the [check policy](docs/developer/check-policy.md): targeted local tests, automatic Python 3.11–3.13 CI, dispatched figure review, and local real-data/performance checks (`scripts/run_real_data_checks.sh`) for affected behavior. Report **executed** evidence, failures, skipped checks and production limits. Short fixtures do not prove convergence.
   
## Code Quality & Anti-Bloat Rules
Before opening a PR or declaring a task complete, audit your changes to remove AI clutter:

1. **Run deterministic checks:**
   - Run `ruff check --select F,I,C90,SIM,RET,UP .` to remove unused imports/variables, enforce modern Python syntax, and flag overly complex or nested logic.
   - Run `vulture . --exclude ".venv,venv,.git,build,dist" --min-confidence 80` to detect unused functions, classes, methods, and variables. Delete dead code rather than commenting it out.
2. **Eliminate unnecessary abstractions:**
   - Inline single-use helper functions or wrapper classes unless they isolate complex logic or external I/O.
   - Avoid speculative configuration flags, unused `**kwargs` pass-throughs, or premature base classes.
3. **Keep error handling lean:**
   - Do not wrap standard operations in broad `try/except Exception:` blocks unless handling a specific, expected failure mode at an I/O or boundary layer.
4. **Strip redundant comments & docstrings:**
   - Delete inline comments that merely narrate *what* the next line does (e.g., `# Loop through items`). Keep comments only when explaining *why* a non-obvious decision, edge case, or domain rule is needed.
   - Keep docstrings concise—omit boilerplate `Args:`/`Returns:` sections on trivial internal functions when type hints already make the signature clear.


## Delivery scope

Work through an authorized milestone to completion in one sustained effort when feasible.
Small commits aid review; they are not stopping points. Continue through related
implementation, integration, tests and issue tracking. Stop for a real blocker or
user decision, not an arbitrary file count, patch size or number of passing tests.
