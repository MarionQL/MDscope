# MDscope agent checklist

Read only the guides relevant to the change; do not duplicate them here.

| Task | Source of truth |
|---|---|
| Ownership, minimal changes, reviews | [Contributing](docs/developer/contributing.md) |
| Module boundaries, persisted flow | [Architecture](docs/developer/architecture.md) |
| Coordinate frames, PBC, units, numerical semantics | [Scientific contracts](docs/developer/scientific-contracts.md) |
| Known scientific limitations | [Scientific audit](docs/developer/scientific-audit.md) |
| Which checks to run and when | [Check policy](docs/developer/check-policy.md) |
| Test effectiveness, remaining coverage | [Validation gaps](docs/developer/validation-gaps.md) |
| SDC4-specific reference cases | [SDC4 fixtures](docs/developer/sdc4-fixtures.md) |
| YAML, observable options, workflows | [Configuration](docs/configuration.md) · [Analyses](docs/analyses.md) · [Workflows](docs/workflows.md) |

## Avoid

- Dead code and stale exports; remove obsolete paths only after checking callers and compatibility.
- Swallowed errors or redundant guards; catch expected failures at boundaries with actionable context.
- Duplicate algorithms, unjustified dependencies, and speculative options.
- Redundant prose; retain concise scientific contracts and non-obvious rationale.

Use **GitHub Flow** and the issue title/acceptance conventions in [Contributing](docs/developer/contributing.md). Keep changes on a focused issue branch and PR, not directly on `main`.

## Optional development skills

Use an installed skill **only when its instructions are available in the
current agent environment** and the task calls for it. Names below are
workflow cues, not a replacement for this repository's authoritative docs.
If the named skill cannot be read, follow the equivalent repository process
and say it was not invoked. Do not block work on a missing skill.

| Task | Prefer skill | Repository contract |
|---|---|---|
| Reproduce a failure, find its cause | **Diagnosing Bugs** | [Validation gaps](docs/developer/validation-gaps.md) and affected entry-point tests |
| Change numerical semantics or fix a regression | **TDD** | [Scientific contracts](docs/developer/scientific-contracts.md) and [Check policy](docs/developer/check-policy.md) |
| Complex feature or cross-module change | **To Spec** → **Implement Spec** | [Contributing](docs/developer/contributing.md); explicit acceptance criteria |
| Review dependencies and module ownership | **Improve Codebase Architecture** / **Codebase Design** | [Architecture](docs/developer/architecture.md), generated dependency graphs and measured evidence |
| Review a PR for science, recovery and maintainability | **Code Review** | [Contributing](docs/developer/contributing.md), [Check policy](docs/developer/check-policy.md) |
| Prepare change summary and evidence | **PR** | GitHub Flow, actual checks, and review status |
| Prioritize multiple bugs/issues | **Triage** | Issue scope/acceptance criteria in [Contributing](docs/developer/contributing.md) |
| Challenge an architectural proposal | **Grill with Docs** | Compare against real producers, saved contracts and developer docs |
| Repeat a class of failures | **Retro** | Record one actionable prevention rule in the canonical guide |
| Establish local commit-time checks | **Setup Pre-Commit** | Add only lightweight, reproducible checks; keep CI authoritative |

For code changes, a good default is **Diagnosing Bugs → TDD →
Code Review → PR**, using **To Spec / Implement Spec** for larger work.
Architecture reviews should recommend bounded changes with evidence before
moving modules. Skills never authorize changing a scientific estimator,
discarding saved outputs, weakening tests, or bypassing review gates.

## Before merging

1. Identify intended behavior, inputs/outputs, units, populations, provenance and failure outcomes. Separate **numerical correctness** from **physical interpretation**.
2. Review six risks: assumptions, mathematical errors, state/data flow, integration, insufficient verification, and failure/recovery.
3. Complete the smallest coherent solution that satisfies the authorized milestone, across modules when needed. Do not interpret minimal changes as a reason to stop after a small patch. Remove superseded internal paths within scope and after migrating callers; use adapters or separate tests when they protect a real contract.
4. For substantive bugs, reproduce failure and test the affected entry point. Prefer independent analytical references; never weaken assertions to pass CI.
5. Validate identity, frame/time alignment, PBC, missing-vs-zero and persisted-product completeness. Existing output files are not automatically valid.
6. Follow the [check policy](docs/developer/check-policy.md): targeted local tests, automatic Python 3.11/3.13 CI, and explicitly dispatched visual/real-data/performance checks for affected behavior. Report **executed** evidence, failures, skipped checks and production limits. Short fixtures do not prove convergence.

## Delivery scope

Work through an authorized milestone to completion in one sustained effort when feasible.
Small commits aid review; they are not stopping points. Continue through related
implementation, integration, tests and issue tracking. Stop for a real blocker or
user decision, not an arbitrary file count, patch size or number of passing tests.
