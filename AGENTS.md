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

The [agent skill-routing catalog](docs/developer/agent-skill-routing.md) maps
all 36 user-listed development/writing skills to MDscope tasks, including
**positive triggers, exclusions, multi-phase routes, and fallback behavior**.
Read the relevant row when selecting a workflow; do **not** assume this repo
or the user's installed-skills screenshot makes a skill callable.

**Availability gate:** inspect the *current agent environment*, find the
actual skill and read its complete instructions before invoking it. Skills
visible in another application may not be exposed here. If unavailable,
perform the equivalent work using the authoritative repository documents
above, explicitly state that the skill was **not invoked**, and never block
progress merely because its instructions are missing. Skill names in the
routing catalog are workflow cues, not verified copies of their instructions.

**Default routing, by task (not a mandatory chain):**

- Unexpected behavior, failed checks, regression: **Diagnosing Bugs**;
  add **TDD** for a behavior-changing fix, then **Code Review** and **PR**
  when relevant.
- New cross-module behavior or scientific output: **To Spec** if contracts
  need clarification; **Domain Modeling** if entities/invariants are unclear;
  **Implement Spec** or **Implement**, backed by scientific tests.
- Structural changes: **Improve Codebase Architecture**, optionally
  **Grill with Docs**; use measured graphs and a bounded refactor rather
  than splitting modules based solely on length.
- Uncertain Studio interaction: **Prototype** or **Wizard** before production
  implementation if a disposable UX spike would reduce risk.
- Multiple issues: **Triage**, then **Wayfinder**/**To Tickets** only if
  sequencing and independent deliverables justify them.
- Agent instructions: **Writing for Agents**; recurring process failure:
  **Retro**. Use **Research** for external technical/scientific evidence.

Select one primary skill per phase rather than invoking every possibly
related skill. Python MDscope development does not automatically call
TypeScript, educational, editorial, migration, or Claude-only skills.
No skill can override [scientific contracts](docs/developer/scientific-contracts.md),
[check policy](docs/developer/check-policy.md), saved-artifact compatibility,
review requirements, or the actual capabilities of the active environment.

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
