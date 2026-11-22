"""Lower packages must not import higher ones (#256).

Intended order: core -> observables/interactions -> mechanisms ->
statistics/reporting -> pipeline -> cli/studio. The `config/` package may
import `core/` and `contracts/`, not execution or presentation layers.
Function-level imports count:
they still couple the lower package to the higher one.
"""
import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "mdscope"
HIGHER = {"mechanisms", "statistics", "reporting", "pipeline", "studio"}
FORBIDDEN = {
    "core": HIGHER | {"observables", "interactions"},
    "observables": HIGHER,
    "interactions": HIGHER,
    "config": HIGHER | {"observables", "interactions", "systems"},
}
# (importing module, imported module): reason. Remove entries as they are fixed.
ALLOWED = {
    ("observables/solvation_postprocess.py", "mechanisms.motif_residence"):
        "generic episode summary still lives in mechanisms; relocate with #258",
}


def _imports(path: Path):
    parts = path.relative_to(PACKAGE).parts
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level:
            base = list(parts[:-1])[:len(parts) - node.level]
            module = ".".join(base + ([node.module] if node.module else []))
        elif (node.module or "").startswith("mdscope."):
            module = node.module.removeprefix("mdscope.")
        else:
            continue
        yield node.lineno, module


def test_lower_packages_do_not_import_higher_ones():
    violations, used = [], set()
    for path in sorted(PACKAGE.rglob("*.py")):
        relative = path.relative_to(PACKAGE)
        layer = relative.parts[0] if len(relative.parts) > 1 else None
        for lineno, module in _imports(path):
            if module.split(".")[0] not in FORBIDDEN.get(layer, ()):
                continue
            key = (relative.as_posix(), module)
            if key in ALLOWED:
                used.add(key)
            else:
                violations.append(f"{relative.as_posix()}:{lineno} imports {module}")
    assert not violations, "\n".join(violations)
    assert used == set(ALLOWED), f"stale ALLOWED entries: {set(ALLOWED) - used}"
