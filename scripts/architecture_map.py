#!/usr/bin/env python3
"""Reproducible, approximate static import/call maps for an installed MDscope checkout.

No MDscope import, trajectory I/O, or third-party graph packages are required.
Import edges include conditional and function-local imports. Call edges include
only statically resolvable local/imported function names and self/cls methods.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path
import platform
import subprocess

from importlib import metadata


def module_name(path: Path, root: Path) -> str:
    parts = path.relative_to(root).with_suffix("").parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join((root.name, *parts))


def import_base(source: str, node: ast.ImportFrom, is_package: bool) -> str:
    if node.level == 0:
        return node.module or ""
    namespace = source.split(".")
    if not is_package:
        namespace.pop()
    for _ in range(node.level - 1):
        if namespace:
            namespace.pop()
    return ".".join((*namespace, *(node.module.split(".") if node.module else ())))


def inspect_package(root: Path) -> dict:
    """Return only provable static relationships; never import any target code."""
    root = Path(root).resolve()
    paths = sorted(root.rglob("*.py"))
    modules = {module_name(p, root): p for p in paths}
    trees = {name: ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
             for name, path in modules.items()}
    imports: dict[str, dict[str, str]] = {}
    module_edges: set[tuple[str, str]] = set()
    functions: dict[str, ast.AST] = {}
    for name, tree in trees.items():
        aliases: dict[str, str] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for item in node.names:
                    alias = item.asname or item.name.split(".")[0]
                    aliases[alias] = item.name if item.asname else item.name.split(".")[0]
                    if item.name in modules and item.name != name:
                        module_edges.add((name, item.name))
            elif isinstance(node, ast.ImportFrom):
                base = import_base(name, node, modules[name].name == "__init__.py")
                if base in modules and base != name:
                    module_edges.add((name, base))
                for item in node.names:
                    if item.name == "*":
                        continue
                    target = ".".join(filter(None, (base, item.name)))
                    aliases[item.asname or item.name] = target
                    if target in modules and target != name:
                        module_edges.add((name, target))
        imports[name] = aliases
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                functions[f"{name}.{node.name}"] = node
            elif isinstance(node, ast.ClassDef):
                for method in node.body:
                    if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        functions[f"{name}.{node.name}.{method.name}"] = method

    calls: set[tuple[str, str]] = set()
    for func_name, func_node in functions.items():
        # An __init__.py package name prefixes child modules: choose the
        # most specific module, not its parent package.
        module = max((m for m in modules if func_name.startswith(m + ".")),
                     key=len)
        local_class = func_name.removeprefix(module + ".").rsplit(".", 1)[0]
        if local_class == func_name.removeprefix(module + "."):
            local_class = None
        aliases = imports[module]
        for node in ast.walk(func_node):
            if not isinstance(node, ast.Call):
                continue
            parts = []
            current = node.func
            while isinstance(current, ast.Attribute):
                parts.append(current.attr)
                current = current.value
            if isinstance(current, ast.Name):
                parts.append(current.id)
            else:
                continue
            parts.reverse()
            if not parts:
                continue
            if parts[0] in ("self", "cls") and local_class and len(parts) == 2:
                candidate = f"{module}.{local_class}.{parts[1]}"
            elif parts[0] in aliases:
                candidate = ".".join((aliases[parts[0]], *parts[1:]))
            else:
                candidate = ".".join((module, *parts))
            if candidate in functions:
                calls.add((func_name, candidate))

    return {
        "modules": sorted(modules),
        "imports": sorted([list(edge) for edge in module_edges]),
        "functions": sorted(functions),
        "calls": sorted([list(edge) for edge in calls]),
        "cycles": cycles(set(modules), module_edges),
        "limitations": [
            "Static source graph, not measured runtime execution.",
            "Conditional and function-local imports are included as potential dependencies.",
            "Call edges are approximate: dynamic dispatch, decorators, reexports, shadowed names, "
            "higher-order functions and imported C extensions are not resolved.",
            "Only module-level functions and immediate class methods are call nodes.",
        ],
    }


def cycles(nodes: set[str], edges: set[tuple[str, str]]) -> list[list[str]]:
    """Return strongly connected components with more than one module."""
    adjacency = {node: [] for node in nodes}
    for source, target in edges:
        adjacency[source].append(target)
    index = 0
    stack = []
    active = set()
    positions = {}
    low = {}
    components = []

    def visit(node):
        nonlocal index
        positions[node] = low[node] = index
        index += 1
        stack.append(node)
        active.add(node)
        for other in sorted(adjacency[node]):
            if other not in positions:
                visit(other)
                low[node] = min(low[node], low[other])
            elif other in active:
                low[node] = min(low[node], positions[other])
        if low[node] == positions[node]:
            part = []
            while stack:
                member = stack.pop()
                active.remove(member)
                part.append(member)
                if member == node:
                    break
            if len(part) > 1:
                components.append(sorted(part))

    for node in sorted(nodes):
        if node not in positions:
            visit(node)
    return sorted(components)


def write_dot(path: Path, name: str, nodes: list[str], edges: list[list[str]]):
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"digraph {name} {{" + "\n  rankdir=LR;\n  overlap=false;\n")
        for node in nodes:
            handle.write("  " + json.dumps(node) + ";\n")
        for source, target in edges:
            handle.write("  " + json.dumps(source) + " -> " +
                         json.dumps(target) + ";\n")
        handle.write("}\n")


def git_revision(repo_root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=False,
        )
    except OSError:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path,
                        default=Path(__file__).resolve().parents[1] / "src/mdscope")
    parser.add_argument("--output", type=Path, required=True,
                        help="Directory for DOT graphs and machine-readable inventory")
    args = parser.parse_args()
    root = args.package_root.resolve()
    if not (root / "__init__.py").is_file():
        parser.error(f"Not a Python package: {root}")
    graph = inspect_package(root)
    try:
        installed_version = metadata.version(root.name)
    except metadata.PackageNotFoundError:
        installed_version = None
    graph["environment"] = {
        "git_sha": git_revision(root.parents[1]),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "installed_distribution_version": installed_version,
        "source": root.as_posix(),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    write_dot(args.output / "modules.dot", "module_imports",
              graph["modules"], graph["imports"])
    write_dot(args.output / "calls.dot", "statically_resolved_calls",
              graph["functions"], graph["calls"])
    (args.output / "architecture.json").write_text(
        json.dumps(graph, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Modules: {len(graph['modules'])}; import edges: {len(graph['imports'])}; "
          f"resolved call edges: {len(graph['calls'])}; "
          f"import SCCs: {len(graph['cycles'])}")
    print(f"Revision: {graph['environment']['git_sha'] or 'unknown'}")
    print(f"Graph files: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
