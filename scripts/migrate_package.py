#!/usr/bin/env python3
"""One-time, idempotent migration from mdscope.analysis to mdscope.

Run from the repository root on issue35/src-layout-packaging.
"""
from __future__ import annotations

import ast
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src/mdscope/analysis"
DEST = ROOT / "src/mdscope"
ALLOWED = {".py", ".toml", ".sh", ".yml", ".yaml", ".md"}
SKIP = {".git", ".venv", "venv", "dist", "build", ".pytest_cache", "__pycache__"}
OLD = "mdscope.analysis"

def migrate() -> None:
    if not (ROOT / "pyproject.toml").is_file():
        raise RuntimeError("Not an MDscope repository")
    if SOURCE.exists():
        for child in sorted(SOURCE.iterdir()):
            if child.name == "__init__.py":
                continue
            target = DEST / child.name
            if target.exists():
                raise RuntimeError(f"Refusing to overwrite existing destination: {target}")
            shutil.move(str(child), str(target))
        old_init = SOURCE / "__init__.py"
        if old_init.exists():
            old_init.unlink()
        SOURCE.rmdir()
        (DEST / "__init__.py").write_text(
            '"""MDscope: reproducible molecular-dynamics analysis workflows."""\n\n'
            '__version__ = "0.12.0"\n\n__all__ = ["__version__"]\n',
            encoding="utf-8",
        )
    elif not (DEST / "observables").is_dir():
        raise RuntimeError("Neither source package nor migrated package found")

    changes = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == Path(__file__).resolve() or path.suffix not in ALLOWED or any(x in SKIP for x in path.relative_to(ROOT).parts):
            continue
        original = path.read_text(encoding="utf-8")
        updated = original.replace("MDscope.analysis", "mdscope").replace(OLD, "mdscope")
        # Repository-relative paths in tooling, CI, and documentation.
        updated = updated.replace("src/mdscope/analysis/", "src/mdscope/")
        if updated != original:
            path.write_text(updated, encoding="utf-8")
            changes.append(str(path.relative_to(ROOT)))
    # Update the packaging assertion introduced by the initial src-layout PR.
    layout_test = ROOT / "tests/test_source_layout.py"
    if layout_test.exists():
        old_assertion = 'files("mdscope") / "analysis" / "cli.py"'
        new_assertion = 'files("mdscope") / "cli.py"'
        content = layout_test.read_text(encoding="utf-8")
        if old_assertion in content:
            layout_test.write_text(content.replace(old_assertion, new_assertion), encoding="utf-8")
            changes.append(str(layout_test.relative_to(ROOT)))
    # Ensure metadata is not accidentally left pointing at the old entry point.
    metadata = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    if 'mdscope = "mdscope.cli:main"' not in metadata:
        raise RuntimeError("CLI entry point was not migrated")
    leftovers = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == Path(__file__).resolve() or path.suffix not in ALLOWED or any(x in SKIP for x in path.relative_to(ROOT).parts):
            continue
        if OLD in path.read_text(encoding="utf-8"):
            leftovers.append(str(path.relative_to(ROOT)))
    if leftovers:
        raise RuntimeError(f"Unmigrated imports: {leftovers}")
    for path in DEST.rglob("*.py"):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    print(f"Migration completed; updated {len(changes)} text files")
    for p in changes:
        print("  ", p)

if __name__ == "__main__":
    migrate()
