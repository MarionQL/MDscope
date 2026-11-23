# Check policy

Run checks for changed behavior after a coherent milestone; do not repeat heavy checks for small edits.

| Change | Local evidence | Before merge |
|---|---|---|
| Docs/comments | Links/examples; compile changed executable snippets | Automatic Python checks; no heavy dispatch |
| Studio/planning/schema metadata | Studio/planner/config tests; check against actual producer | Python matrix; real data only if execution/preflight changes |
| Kernels/selections/coordinates/streaming/recovery | Independent references and entry-point/restart tests | Python matrix + relevant local real-data check |
| Statistics | Reference/family/dispatcher tests; missingness/replica checks | Python matrix; real data where the affected path requires it |
| Figures/routing/layout/labels | Inspect production-renderer synthetic previews | Python matrix + dispatched figure review |
| Packaging/dependencies/CLI | Installed import/CLI outside checkout | Matrix, wheel and sdist checks |
| Performance/hot paths | Numerical parity, representative comparison | Relevant benchmark; record CPU, RSS and I/O |

## Development commands

```bash
python -m pip install -e '.[test]'  # add sasa for native FreeSASA
python -m pytest -q tests          # default excludes benchmark tests
ruff check src tests scripts      # CI-blocking E9/F rules
ruff check --select F,I,C90,SIM,RET,UP .
vulture . --exclude ".venv,venv,.git,build,dist" --min-confidence 80
git ls-files -z '*.sh' | xargs -0 shellcheck
pytest -q -s -m benchmark tests/test_performance_benchmarks.py
```

- Review broader Ruff/Vulture findings within scope; remove confirmed dead code after checking callers and public contracts.
- Intentional re-exports: `import name as name`.
- Keep shell/checksum files LF (`.gitattributes`); Slurm rejects CRLF. On Windows, write with `newline="\n"`.

## Automation

- `MDscope checks`: Python 3.11–3.13, 3.11 minimum dependencies, packaging and `bash -n` for tracked shell scripts. Superseded runs cancel.
- The 3.11 job also runs ShellCheck/Ruff and publishes `mdscope-coverage` (no minimum).
- Workflows use `contents: read`; external actions are SHA-pinned and updated by Dependabot. Review release notes before merging bumps.
- Figure review is the only `workflow_dispatch` job: Actions → MDscope figure review → Run workflow → PR branch. Review instructions: [Contributing](contributing.md#visual-review).
- Real-data checks run locally; other scientific audits remain manual. No scheduled production computation.
- PR path filters see the cumulative diff, so a later docs-only push cannot reliably suppress earlier code-related jobs.

## Real GROMACS checks

```bash
python -m pip install -e '.[sasa]'  # Linux/WSL or ARC; no Windows FreeSASA wheel
scripts/run_real_data_checks.sh smoke
# Other checks: hbond solvation performance finalization box-sweep box-clearance all
```

- Full inputs belong in git-ignored `examples/data/`; outputs default to git-ignored `real_test_results/`. Overrides: `MDSCOPE_REAL_DATA`, `MDSCOPE_REAL_OUTPUT`.
- The runner checks full inputs against `tests/data/real_test_data.sha256`.
- [GitHub test-data assets](https://github.com/MarionQL/MDscope/releases/tag/v0.1.0) include restored short XTCs (41/40 frames, 50 ns spacing). They differ from the full-resolution checksum pins and cannot establish full performance/release validation. Do not replace the pins with short-input hashes.
- For derived smoke iteration, `python scripts/stride_real_data.py --stride 5` writes `examples/data/stride5/`, `provenance.json` and `SHA256SUMS`. Use `MDSCOPE_REAL_DATA=examples/data/stride5`; only `smoke`/`hbond` are allowed. Short inputs do not replace full-data merge evidence.
- Record commit, environment, command, result and input qualification in the PR. Full benchmark windows must fit the source trajectory.

## Evidence reuse

- Reuse passed evidence while relevant code, dependencies and inputs are unchanged.
- After edits, rerun affected checks; broaden for failures, integration changes or unresolved risks. CI covers the Python matrix; targeted success alone does not warrant repeating the full local suite.
- Docs/UI-only edits do not invalidate unchanged numerical evidence. Numerical or visual changes invalidate their affected evidence.
- Report failed, skipped, unavailable and resource-blocked checks explicitly; none counts as a pass. Required checks remain merge obligations.
- Before changing workflow filters, check required statuses. Never bypass a relevant failure with skip markers or weaker assertions.
