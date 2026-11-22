"""Keep runtime environment variables discoverable in the workflow guide."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "workflows.md"

# Only inspect variable *uses*, not comments, YAML examples, or synthetic atom names.
PYTHON_ENV = re.compile(r"""(?:os\.environ\.(?:get|__getitem__)|os\.getenv)\(\s*["'](MDSCOPE_[A-Z0-9_]+)["']""")
SHELL_ENV = re.compile(r"""\$\{?(MDSCOPE_[A-Z0-9_]+)\b|\b(?:export\s+)?(MDSCOPE_[A-Z0-9_]+)\s*=""")
SHELL_PYTHONUNBUFFERED = re.compile(r"\bPYTHONUNBUFFERED\s*=")


def test_documented_runtime_environment_variables():
    found = set()
    for dirname in ("src", "scripts", "slurm"):
        for path in (ROOT / dirname).rglob("*"):
            if path.suffix not in {".py", ".sh"}:
                continue
            source = path.read_text(encoding="utf-8")
            if path.suffix == ".py":
                found.update(PYTHON_ENV.findall(source))
            else:
                found.update(v for match in SHELL_ENV.findall(source)
                             for v in match if v)
                if SHELL_PYTHONUNBUFFERED.search(source):
                    found.add("PYTHONUNBUFFERED")
    documented = DOC.read_text(encoding="utf-8")
    missing = sorted(name for name in found if f"`{name}`" not in documented)
    assert not missing, f"Undocumented runtime environment variables: {missing}"
