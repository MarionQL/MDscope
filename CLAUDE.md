# Claude Code

Read [AGENTS.md](AGENTS.md) for task routing and delivery requirements.

- Development commands and test selection: [Check policy](docs/developer/check-policy.md).
- Module ownership and graph generation: [Architecture](docs/developer/architecture.md).
- CLI test subprocesses use the running interpreter's installed console script. For an uninstalled worktree, use that environment's `python -m pytest` with `PYTHONPATH=<worktree>/src` so subprocesses exercise the intended source.
