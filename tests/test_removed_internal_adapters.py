"""Internal code must use canonical owners after adapter removal."""
import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "mdscope"
REMOVED = {
    "statistics.state_feature_spec",
    "mechanisms.shared_states",
    "contracts.residue_response",
    "config.validators.shared_states",
    "config.validators.hydration",
    "config.validators.residue_response",
}


def test_no_imports_from_removed_internal_adapters():
    failures = []
    for path in PACKAGE.rglob("*.py"):
        relative = path.relative_to(PACKAGE)
        parts = relative.parts[:-1]
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.level:
                target = ".".join((*parts[:len(parts) - node.level + 1],
                                   *([node.module] if node.module else [])))
            else:
                target = (node.module or "").removeprefix("mdscope.")
            if target in REMOVED:
                failures.append(f"{relative}:{node.lineno} imports {target}")
    assert not failures, "\\n".join(failures)
