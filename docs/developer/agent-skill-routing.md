# Agent skill routing for MDscope

Read this guide when choosing a development skill. It records actual instruction review, task triggers and limits; it is not an installation manifest.

## Availability and evidence

On **2026-11-20 (America/New_York)**, all **40 catalog entries** were readable in this ChatGPT Work session: **39 personal skills plus the platform/plugin skill-creator**. Each main `SKILL.md` was read through the skill reader with no continuation outstanding. This supersedes the earlier screenshot-only / none-readable assessment in this PR. There are **no missing main instruction files among these 40**.

The audit covered main instructions, not every linked template, script or reference, and did not execute all 40 workflows. Read the branch-specific supporting resources when actually running a skill. Readability does not prove a CLI, background worker, tracker operation or external service is usable. Other sessions must discover their own capabilities; repository text does not expose or install a skill.

## Selection and repository integration

1. Match the requested outcome to one primary workflow; load helper references and required dependencies when that workflow calls for them. A workflow can legitimately require several skills.
2. Respect invocation intent: **Auto** means relevant task-driven use is allowed; **Request** reflects `disable-model-invocation: true` and requires the user to request that workflow, by name or clear intent. Suggest a request-only workflow when useful instead of silently entering it. **Managed** is the platform skill manager. Auditing a skill reads its instructions without running its procedure.
3. Use [contributing](contributing.md) for GitHub Flow, issue titles, labels and review, [scientific contracts](scientific-contracts.md) for numerical semantics, and [check policy](check-policy.md) for checks and cadence. Generic skill advice must preserve independent scientific references, estimator distinctions, persisted compatibility and meaningful validation.
4. Resolve tracker operations through MDscope's existing GitHub workflow and the available connector or authenticated CLI. Inspect any `docs/agents/` configuration before using tracker-dependent skills. Missing optional skill setup is not a reason to interrupt a bounded task already supported by repository guidance; report a manual fallback accurately. Adding skill-family setup, new labels or a parallel glossary requires its own justified scope.
5. Apply existing authorization to decisions already settled in the session. Keep genuine user decisions explicit; do factual investigation yourself. A skill's phase boundary is not a reason to stop an authorized implementation milestone. Planning/interview workflows remain planning/interview workflows when those are what the user requested.
6. If required agents, templates or runtime capabilities are unavailable, identify the specific missing capability and follow repository practice where possible. Describe partial/manual work as such. Never claim a full invocation, independent review, test pass, scheduled task or background execution from a read alone.

For TDD and spec synthesis, confirm proposed test interfaces when they have not already been agreed. Use the existing check policy for when to broaden runs; generic requests for frequent/full checks do not add heavy dispatches to docs or small authoring edits.

## Skill-by-skill matrix

Names below are exact callable names; display names in the UI may use title case. All rows have main instructions **readable/read** in the dated audit.

| Skill | Mode | Trigger in MDscope | Actual behavior and limits |
|---|---|---|---|
| `ask-matt` | Request | Ask which workflow fits an unfamiliar task. | Router over this skill family; read the chosen skill before recommending it. Routine routing already covered here needs no interview. |
| `code-review` | Auto | Review a PR or substantive diff against a fixed base. | Two separate reviews: repository standards and originating spec, using parallel fresh agents. Resolve the base and spec first; report axes separately. A manual review is not the full skill. |
| `codebase-design` | Auto | Design a module interface, test seam, or bounded architectural improvement. | Reference vocabulary and deletion test for deep modules; not a mandate to split files or introduce abstractions. |
| `diagnosing-bugs` | Auto | Reported failure, numerical regression, stale Studio state, or measured slowdown. | First run a fast symptom-specific red-capable command; minimise, test ranked hypotheses, fix and recheck the original scenario. Document missing access or an inadequate regression seam. |
| `domain-modeling` | Auto | Resolve overloaded condition/replica/selection/artifact terminology or document a consequential decision. | Actively sharpen terms; update glossary as terms settle. ADRs require a costly reversal, non-obvious rationale and a real tradeoff. Reading a glossary alone is not invocation. |
| `implement` | Request | Implement one defined piece of work or ticket. | Fetch referenced ticket, use TDD at agreed seams where applicable, review, commit. Apply MDscope check policy to validation cadence; this is not the parallel whole-spec orchestrator. |
| `grill-with-docs` | Request | Sharpen an unresolved design with the user while retaining domain decisions. | Combines grilling and domain-modeling. Produces glossary/ADR updates as warranted; not a passive source audit or automatic prerequisite to a known fix. |
| `implement-spec` | Request | Implement an entire spec with an existing dependency graph of tickets. | Fresh implementer agents in isolated worktrees, frontier scheduling, integration branch, merger agents and final two-axis review. Reading a spec is insufficient; verify tracker and agent capabilities. |
| `improve-codebase-architecture` | Request | Survey architectural friction and choose a deepening opportunity. | Inspect hot spots, glossary and ADRs; fresh exploration agent; visual HTML candidate report; user selects a candidate before grilling. File length/graphs alone do not constitute this workflow. |
| `pr` | Auto | Write or revise a PR body. | Smallest useful visual/diff, before/after evidence, reversible-door and blast-radius assessment. Does not open, test or merge a PR by itself. |
| `prototype` | Auto | Answer a specific uncertain logic/state or Studio appearance question with runnable disposable code. | Logic: one interactive HTML demo. UI: distinct visual variants. Capture question/verdict and prototype on a separate branch with an issue pointer; production retains validated decisions. |
| `research` | Auto | Investigate external scientific methods, library behavior or official API facts. | Fresh background researcher, primary sources, cited Markdown findings following repo convention. Requires actual delegation; a direct source read can be useful without claiming this full skill ran. |
| `retro` | Request | Review a coding session to prevent repeated process failures. | Read session evidence, propose environment improvements ranked by severity; mechanical failures belong in deterministic checks, judgement in review standards. Not automatic code cleanup after every patch. |
| `setup-matt-pocock-skills` | Request | Configure this engineering skill family's repository integration. | Configures tracker, triage label mapping and domain-doc layout; not TypeScript-specific and not a skill installer. Reuse existing MDscope conventions and confirm proposed configuration changes. |
| `tdd` | Auto | Build or fix observable behavior test-first, especially scientific or persistence contracts. | Agree public seams, then one meaningful red/green slice at a time with independent expected values. Honor already-agreed scope/seams; docs-only edits do not need behavior tests. |
| `to-spec` | Request | Synthesize a discussed feature into a publishable, testable spec. | No new requirements interview; explore existing code, agree test seams and synthesize the conversation. Publishes to configured tracker; does not implement or silently settle unresolved user decisions. |
| `triage` | Request | Process incoming raw issues or explicitly named external PR requests into agent-ready briefs. | Verify claims and redundancy, recommend category/state, use configured labels and required AI triage disclaimer. Generated implementation tickets are already prepared; ordinary backlog prioritization need not run this state machine. |
| `to-tickets` | Request | Decompose a discussed plan/spec into independently verifiable work. | Vertical slices with explicit blocking edges; confirm breakdown before publishing. Wide mechanical refactors may use expand/contract. Preserve parent issue; reuse existing labels. |
| `wayfinder` | Request | Chart an initiative whose unresolved decisions exceed one session. | Issue-backed decision map; research/prototype/grilling/task tickets, claims and dependencies. Planning by default, one non-research decision ticket per session; hand off to a spec once clear. Not an implementation milestone loop. |
| `wizard` | Auto | Guide steps only the human can perform, such as third-party credentials or manual cutover. | Generate an interactive bash procedure using the supplied template, ordered stages, captured values and hidden secret entry. Confirm stages; not a Studio first-run UI builder or agent-runnable setup substitute. |
| `chief-of-staff` | Request | Coordinate a sustained strategic goal through delegated work and durable knowledge. | All work in fresh subagents; maintain dated, sourced wiki claims and consider useful schedules. Agree a new wiki location. Ordinary milestone tracking does not warrant this whole operating model. |
| `claude-handoff` | Request | Transfer work immediately to an actual Claude background worker. | Writes a redacted summary and launches claude --bg with a descriptive name. Readable here, but execution requires that CLI, authentication and background-worker support; generic Codex agents are not proof it ran. |
| `loop-me` | Request | Design a reusable workflow for a recurring pattern in the user's work. | Stateful grilling producing workflows/*.md and contextual notes. Defines trigger, checkpoints and deliverables where needed; does not install a skill, execute the workflow or create a schedule. |
| `setup-ts-deep-modules` | Request | Enforce private implementation imports in an actual TypeScript subproject. | Installs dependency-cruiser, documents entry points and proves pass/fail/pass on a forbidden import. Inapplicable to Python MDscope by default. |
| `software-factory` | Request | Design an agent system that turns tracker/git events into reviewed code. | Shared vocabulary for tasks, runs, runners and external service slots; explicitly design triggers, deduplication, brakes and human re-entry. A reference/design skill, not an end-to-end build executor. |
| `writing-beats` | Request | Turn fixed raw writing material into a chosen journey of beats. | Set audience prerequisites; offer alternatives and append only the selected beat. Interactive article composition, not an ordinary PR summary. |
| `writing-fragments` | Request | Explore ideas for an article by collecting raw fragments. | Interview and append unstructured fragments to one file; reread before writes and preserve user edits. Intentionally does not organize an article. |
| `writing-shape` | Request | Shape supplied raw material into a separate article block by block. | Read input in full, settle prerequisites and opening, append agreed blocks. Preserve raw source; user chooses completion. Not a generic automatic documentation formatter. |
| `git-guardrails-claude-code` | Auto | Explicitly configure Claude Code hooks against dangerous git operations. | Project/global scope choice, merge hook settings, test blocking. Blocks every git push as well as destructive commands; discuss that consequence. Does not protect GitHub connector writes. |
| `migrate-to-shoehorn` | Auto | Requested TypeScript test-data migration from as assertions. | Use @total-typescript/shoehorn only in TS tests; verify typecheck. Does not apply to pytest or production scientific calculations. |
| `setup-pre-commit` | Auto | Requested Husky/lint-staged/Prettier setup in a Node project. | Adds Node dependencies and commit-time typecheck/test scripts, potentially full suites. Not Python pre-commit guidance; use MDscope-native tooling for Python and keep check-policy cadence. |
| `scaffold-exercises` | Auto | Requested course exercise scaffolding in a compatible ai-hero course workspace. | Specific numbered problem/solution/explainer layout, main.ts convention and pnpm ai-hero-cli linter. Not a generic Python tutorial generator; do not add that toolchain to MDscope. |
| `grill-me` | Request | Stress-test a plan without retaining repository domain docs. | Stateless wrapper around grilling. For repository design that should retain terms/decisions, prefer grill-with-docs when requested. |
| `grilling` | Auto | User asks to stress-test a decision through an interview. | Ask independent frontier questions in rounds with recommendations, wait for decisions and look up facts through agents. Not routine PR review or an instruction to interrogate already-specified work. |
| `handoff` | Request | Move work to another harness, directory, colleague or separate session. | Portable redacted state summary with suggested skills and pointers to existing artifacts. Writes a document; launches no worker. |
| `teach` | Request | User requests sustained learning about an MDscope-related topic. | Stateful teaching workspace with mission, sourced lessons, references and learning records. Keep learning artifacts in the agreed teaching workspace; not automatic explanatory prose during implementation. |
| `to-questionnaire` | Request | A collaborator holds facts needed for a decision. | Clarify recipient and required answers, produce an async discovery document covering each gap. Drafting a questionnaire does not authorize sending it. |
| `wait-what` | Request | User says the previous explanation did not land. | Re-pitch the last explanation with context, simplified technical English and domain vocabulary. Not a questionnaire or an early product-planning gate. |
| `writing-for-agents` | Auto | Edit AGENTS.md, CLAUDE.md, skill instructions or linked agent guidance. | Use branch-specific context pointers, progressive disclosure, clear completion criteria and one authoritative home per rule. Keep detailed skill routing here rather than duplicating it in AGENTS.md. |
| `skill-creator` | Managed | User requests creating, installing, editing or managing a skill, or asks how skills work. | Use the platform-provided skill manager, resolving its active package at runtime. Repository routing docs do not install skills. This name is supplied by the platform/plugin rather than the 39 personal packages. |

## Common task routes

These are choices by task, not a mandatory chain. Request-only workflows retain the invocation rule above.

| Task | Start with | Follow-up when needed |
|---|---|---|
| Failed topology loading, invalid resume, numerical regression | diagnosing-bugs: command that catches the reported symptom | tdd at an agreed public seam; two-axis code-review for substantive changes; pr for the narrative |
| Defined scientific metric or Studio feature | Existing issue and contracts; implement when requested | tdd and code-review; research for unresolved external facts |
| Discussed feature needs a durable spec | to-spec when requested | to-tickets for a multi-ticket build; implement-spec only for an authorized parallel task graph |
| Uncertain Studio logic or appearance | prototype, with one explicit question | Fold the validated decision into the production task; wizard is reserved for human-only procedures |
| Architecture improvement | codebase-design for a selected interface | improve-codebase-architecture when an interactive candidate survey is requested |
| Incoming raw reports | triage when requested | Prepared tickets proceed directly to implementation |
| Huge unresolved initiative | wayfinder when requested | to-spec once decisions are clear; not direct implementation of a decision map |
| Agent policy change | writing-for-agents | Keep detailed routing in this document; AGENTS.md holds a short pointer |
| Unclear explanation | wait-what when requested | Re-explain the same work with necessary context |
| Python developer hooks | Existing Python tooling and check policy | setup-pre-commit applies only to a separately scoped Node/Husky setup |

## Reporting

Keep evidence distinct: **instructions read**, **workflow executed**, **manual adaptation**, **tests run**, **CI observed**, and **PR updated**. Record actual commands, commits or links for verification. Read this catalog on demand; load only the chosen skill and resources needed for the current phase.
