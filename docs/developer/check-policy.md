# Check policy

Run checks for the behavior changed, not every available workflow.

| Change | Local evidence | Before merge |
|---|---|---|
| Docs/comments only | Review links/examples; compile if executable snippets changed | Automatic Python checks; no heavy dispatch |
| Studio/planning/schema metadata | Studio/planner/config tests; verify contract against actual producer | Automatic Python matrix; no real-data run unless execution/preflight semantics changed |
| Scientific kernels, selections, coordinates, streaming or recovery | Affected independent references and entry-point/restart tests | Python matrix + relevant real GROMACS smoke/parity workflow |
| Statistics estimators, selection or aggregation | Reference and family/dispatcher tests; missingness/replica checks | Python matrix; real-data validation if affected data path needs it |
| Renderer, figure routing, layout or labels | Production-renderer synthetic previews; inspect selected figures | Python matrix + dispatched figure review, visually inspect artifact |
| Packaging/dependencies/CLI | Installed import/CLI checks outside checkout | Existing matrix/build/wheel/sdist checks; heavy jobs only if their runtime is affected |
| Performance claim or major hot-path change | Numerical parity first; representative measured comparison | Explicit relevant benchmark dispatch; record CPU/RSS/I/O, not just toy timing |

## Automation

- `MDscope checks` stays automatic: Python 3.11/3.13 regressions and existing
  packaging checks protect cross-module integration. Superseded runs cancel.
- Figure review, real GROMACS smoke, release H-bond parity, release solvation
  benchmark and pipeline/statistics benchmark are **workflow_dispatch** jobs.
  Run them through Actions → select workflow → Run workflow → select PR branch.
- Other existing manual scientific audits remain manual. No scheduled production
  computation is added by this policy.
- Dispatch heavy checks once a coherent change is ready for review. Do not run a
  full trajectory download for every small authoring or documentation commit.

GitHub PR path filters use the cumulative PR diff, not just the last commit.
They cannot reliably express “this push only changed Studio” once the PR has
also touched config/runner files. Explicit dispatch avoids that repeated work.

## Reuse and reporting

- Record checked commit, command/workflow, outcome and limitations in the PR.
- Reuse passed evidence while its relevant implementation, dependencies and test
  inputs remain unchanged. After edits, rerun affected tests; broaden for failures,
  integration changes or unresolved risks. Do not repeat the full local suite
  merely because targeted tests passed; CI already covers the Python matrix.
- A later docs/UI-only commit does not invalidate numerical evidence from the
  same unchanged kernels. A changed numerical path does: rerun its real/reference
  checks. Visual changes similarly invalidate the affected visual evidence.
- A skipped/manual job is not a pass. Explain why a heavy check was not applicable.
  Required checks in the table remain merge obligations even when dispatch is manual.
- Before workflow-filter changes merge, verify repository-required status checks
  do not demand a now-manual job; adjust requirements explicitly if necessary.
  Never use skip markers or weaken assertions to bypass a failing relevant check.
