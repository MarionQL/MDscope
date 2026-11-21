# Agent skill routing for MDscope

This is a **task-to-skill decision guide**, not a declaration that any skill is installed or callable. The 36 names below come from the user-provided installed-skills screenshot (2026-11-20). **Availability is per execution environment**: discover the named skill, read its actual instructions, and only then invoke it. The ChatGPT skill registry inspected on 2026-11-20 exposed **none of these 36 skills as individually readable** in this chat. This observation does not establish whether Claude Code, another session, or a future environment can read them.

## Invariants and selection policy

1. Identify task, expected deliverable, relevant repository contracts, and whether the skill is actually exposed. Never interpret the screenshot, this catalog, or a skill's title as executable instructions.
2. Use **one primary skill per phase** and add a second only when it changes the outcome. Avoid indiscriminate chains; skipping a skill must not block work. If instructions are missing, implement the same task using canonical repository guidance, and report that the named skill was **not invoked**.
3. Scientific contracts, checks, and producer/consumer schemas take precedence over generic skill advice: [scientific-contracts](scientific-contracts.md), [architecture](architecture.md), [contributing](contributing.md), [check-policy](check-policy.md), [validation-gaps](validation-gaps.md). Do not let a skill alter estimators or schemas without explicit scope, or weaken checks.
4. Distinguish **writing a plan**, **performing changes**, **running tests**, and **claiming verification**. A skill invocation does not imply tests ran, CI passed, a PR was created, or human review occurred.
5. Do not automatically use speculative product planning, teaching, content writing, background delegation, or TypeScript-specific skills for ordinary Python/Streamlit maintenance. Choose them only if task and runtime warrant.
6. For larger tasks, track acceptance criteria and implement the full authorized milestone; checkpoints are for coherence, not a reason to stop after an arbitrary number of commits.

## Skill-by-skill decision matrix

**Status:** *Screenshot only* for every entry below in the inspected ChatGPT environment. These are recommended triggers based on the visible names/descriptions, **not claims about the uninspected skill internals**. If the instructions become available, re-evaluate the trigger against them.

| Skill | Trigger in MDscope | Don't invoke just because… / output |
|---|---|---|
| **Triage** | Several open failures/issues need grouping, severity, ownership, or ordered acceptance work | One well-scoped bug already has an owner; output prioritized issues and evidence |
| **Wayfinder** | A large, multi-PR initiative needs a navigation/decision trail and explicit stage gates | Simple patch or single bounded milestone; output phased decisions and dependencies |
| **To Spec** | Behavior spans CLI, pipeline, Studio, persistence, or requires explicit contracts before implementation | Acceptance criteria are already unambiguous; output testable specification |
| **PR** | Code/docs changed and a reviewer-readable PR narrative is needed | It is not a substitute for running checks or opening an actual PR; output summary, risks and verification |
| **Wizard** | User explicitly wants an interactive setup wizard (e.g., Studio first-run inputs) | Ordinary form refinements or analysis options; output user-flow prototype/spec, not a requirement to change frameworks |
| **Teach** | User requests a guided lesson about MD, Python, architecture, statistics or code | Routine implementation or review; output learning progression |
| **Handoff** | Long or interrupted session requires a concise state transfer with branch, changes and blockers | Ordinary progress note; output factual, restartable handoff |
| **Setup TS Deep Modules** | A real TypeScript subproject needs deep module boundary enforcement | Main MDscope package is Python; do not introduce TS merely to use a skill |
| **Writing for Agents** | Editing AGENTS.md or agent-facing instructions, including this routing policy | End-user prose or program docs; output concise, discoverable, nonduplicative agent policy |
| **Retro** | Repeated incidents, flaky checks or workflow friction warrant one preventative correction | Every individual fix; output concrete cause, guardrail, owner |
| **Wait What** | A proposed plan/pitch lacks critical context and needs targeted missing-context questions | Existing issue has adequate acceptance criteria; do not block task on broad questionnaires |
| **To Questionnaire** | User explicitly needs a reusable question document for stakeholders or maintainers | One or two readily answered clarifications |
| **Implement Spec** | Approved spec or tickets exist; work needs production code, integration and tests | No specification and high uncertainty; use To Spec first when needed |
| **Code Review** | Nontrivial diff/PR merits independent correctness, scientific contract, lifecycle and recovery review | No diff exists; prioritize severity-ranked findings and affected files |
| **software factory** | User explicitly asks for an end-to-end structured product-building process | Single MDscope issue, especially when repository process already suffices |
| **Domain Modeling** | Unclear semantics for conditions, replicas, selections, outputs, saved artifacts or state transitions | Renaming a UI label; output glossary/invariants, not a parallel registry |
| **Research** | Need evidence from external scientific literature, algorithms, or competing implementations | Internal source search/bug reproduction; cite sources and distinguish evidence from inference |
| **Writing Shape** | User wants a coherent user-facing article explaining a complex scientific/software topic | Routine docstring or PR narrative |
| **Scaffold Exercises** | User asks for educational coding exercises/tests to learn MDscope internals | Do not generate practice tasks during production implementation |
| **Setup Matt Pocock Skills** | User explicitly configures those skills in a TypeScript-oriented repo | Not relevant to Python MDscope by default |
| **Diagnosing Bugs** | A reproducible exception, regression, CI failure, corrupted persisted output, or unexpected numerical value | Merely implementing a known feature; establish reproduction, root cause, fix and regression evidence |
| **Chief Of Staff** | A long-running, goal-oriented working session needs prioritized milestones and decision tracking | One patch or simple PR summary; output focus and open decisions |
| **Improve Codebase Architecture** | Evidence shows change coupling, bad boundaries, cycles, inconsistent state owners, or mixed UI/backend concerns | File length alone; start with dependency graphs and bounded seam candidates |
| **Ask Matt** | User specifically asks which Matt Pocock skill/workflow suits an unfamiliar task | Default routing already identifies a good fit |
| **To Tickets** | A reviewed spec must be decomposed into independent testable GitHub issues/tasks | A single implementation can be completed coherently; output dependencies and acceptance checks |
| **Prototype** | Uncertain UX or architecture hypothesis benefits from a disposable spike (e.g., selection builder interaction) | Production code is already specified; do not silently merge throwaway work |
| **Grilling** | Early design or technical premise needs systematic stress-testing before commitment | Routine completed PR review; output unresolved risks and decisions |
| **Setup Pre-Commit** | Local developer checks would prevent common, cheap formatting/lint/guard failures | Expensive MD/real-data runs or duplicated CI; keep hooks fast and opt-in where appropriate |
| **Migrate to Shoehorn** | An actual repo migration from test assertions to shoehorn is requested with a verified benefit | Existing pytest suite works; do not change test tooling opportunistically |
| **Codebase Design** | Need vocabulary and design boundaries for deep modules or dependency architecture | It is not itself a mandate to restructure working code |
| **Claude Handoff** | Delegating work to a Claude background agent when such agent capability actually exists | Current chat has no background worker; never promise unattended progress without a real invoked agent |
| **Writing Fragments** | User supplies rough material to organize into reusable writing fragments | Engineering source code or standalone analysis |
| **Git Guardrails for Claude Code** | Claude Code environment needs explicit protection against destructive git commands | Generic GitHub connector work; current repository policies still govern safe operations |
| **Grill Me** | A decision is underdefined and user requests a focused interview to clarify it | Requirements are sufficient; do not stall implementation with questions |
| **Grill with Docs** | Competing architecture claims require critical reading against repo source and actual contracts | General brainstorming without source material; document source-backed findings |
| **Loop Me** | User wants reusable skill/workflow specification for recurring work | Single-use action; do not equate with autonomous recurring scheduling |
| **Writing Beats** | User explicitly builds article/story structure from supplied material | PR body or code review checklist |
| **TDD** | Fixing regression or changing scientific calculations, state transitions, validators or persisted contracts | Docs-only changes; establish failing/characterizing test then make it pass with meaningful assertions |
| **Implement** | Defined issue/spec/ticket is ready for coding, integration, validation and reporting | Research-only or speculative discussion; honor actual scope |
| **skill-creator** | User asks to design/update an installable skill and its instructions | Editing repository AGENTS.md alone does not install or expose a skill |

## Common MDscope routes

| Scenario | Primary route | Additional phase only if useful |
|---|---|---|
| Studio input-path bug or stale topology state | Diagnosing Bugs → TDD → Implement | Code Review → PR |
| New scientific metric or changed statistical inference | To Spec → Domain Modeling (if concepts unclear) → TDD → Implement Spec | Code Review, Research for external validation |
| Architecture audit / tangled dependencies | Improve Codebase Architecture → Grill with Docs | To Spec → To Tickets if authorized refactor is large |
| Large Studio UX interaction | Prototype or Wizard for uncertain user flow → To Spec | Implement → Code Review |
| Multiple issues and release scope | Triage → Wayfinder or Chief Of Staff | To Tickets only for independent work |
| Release/check failure | Diagnosing Bugs → TDD | Retro if failure repeats |
| PR ready for review | Code Review → PR | Handoff only if work continues in another session |
| Scientific literature question | Research | Writing Shape only when producing an explanatory article |
| Agent policy update | Writing for Agents | Retro if addressing observed process failure |

## Execution and reporting contract

At the start of an issue, check the applicable repo docs and discover the capability actually available in that environment. If a relevant skill exists, **read its full instructions** before use; if unavailable, follow normal engineering practice without pretending to have invoked it. When reporting, distinguish: `skill available/read/invoked`, `workflow manually followed`, `tests run`, `CI observed`, `PR updated`. Include actual links/check identifiers when available. Do not copy generic skill bodies or every table row into AGENTS.md; this guide is the single routing catalog.
