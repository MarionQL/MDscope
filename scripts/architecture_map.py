#!/usr/bin/env python3
"""Reproducible, approximate static import/call maps for an MDscope checkout.

``--output DIR`` writes the JSON inventory and DOT graphs; ``--write-docs``
regenerates the package-dependency diagram in docs/developer/architecture.md
and ``--check-docs`` fails if it is stale (a test runs the check).

No MDscope import, trajectory I/O, or third-party graph packages are required.
Import edges include conditional and function-local imports. Call edges include
only statically resolvable local/imported function names and self/cls methods.
"""
from __future__ import annotations

import argparse
import ast
import json
import platform
import subprocess
from importlib import metadata
from pathlib import Path


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
    import_sites: set[tuple[str, str, int, str]] = set()
    functions: dict[str, ast.AST] = {}
    for name, tree in trees.items():
        aliases: dict[str, str] = {}
        parents = {child: parent for parent in ast.walk(tree)
                   for child in ast.iter_child_nodes(parent)}

        def record_import(target: str, node: ast.AST) -> None:
            if target not in modules or target == name:
                return
            module_edges.add((name, target))
            ancestor = node
            phase = "module_load"
            while ancestor in parents:
                ancestor = parents[ancestor]
                if isinstance(ancestor, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                    phase = "deferred"
                    break
            import_sites.add((name, target, node.lineno, phase))

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for item in node.names:
                    alias = item.asname or item.name.split(".")[0]
                    aliases[alias] = item.name if item.asname else item.name.split(".")[0]
                    record_import(item.name, node)
            elif isinstance(node, ast.ImportFrom):
                base = import_base(name, node, modules[name].name == "__init__.py")
                record_import(base, node)
                for item in node.names:
                    if item.name == "*":
                        continue
                    target = ".".join(filter(None, (base, item.name)))
                    aliases[item.asname or item.name] = target
                    record_import(target, node)
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
        "import_sites": [
            {"source": source, "target": target, "line": line, "phase": phase}
            for source, target, line, phase in sorted(import_sites)
        ],
        "functions": sorted(functions),
        "calls": sorted([list(edge) for edge in calls]),
        "cycles": cycles(set(modules), module_edges),
        "limitations": [
            "Static source graph, not measured runtime execution.",
            "Conditional and function-local imports are included as potential dependencies.",
            "Import-site phase is module_load or deferred; module_load may still be conditional.",
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


REPO = Path(__file__).resolve().parents[1]
ARCHITECTURE_DOC = REPO / "docs/developer/architecture.md"
LAYERING_TEST = REPO / "tests/test_import_layering.py"
BEGIN = ("<!-- BEGIN GENERATED package dependencies; "
         "regenerate with: python scripts/architecture_map.py --write-docs -->")
END = "<!-- END GENERATED package dependencies -->"


def layering_exceptions(test_path: Path = LAYERING_TEST) -> set[tuple[str, str]]:
    """Module imports the layering test allows, as (importer, imported) names."""
    tree = ast.parse(test_path.read_text(encoding="utf-8"))
    allowed = next(ast.literal_eval(node.value) for node in tree.body
                   if isinstance(node, ast.Assign)
                   and any(getattr(t, "id", None) == "ALLOWED" for t in node.targets))
    return {("mdscope." + path.removesuffix(".py").replace("/", "."), "mdscope." + module)
            for path, module in allowed}


def _unit(module: str) -> str:
    """The top-level subpackage or module that contains *module*."""
    parts = module.split(".")
    return parts[1] if len(parts) > 1 else parts[0]


def package_dependency_mermaid(graph: dict, exceptions: set[tuple[str, str]]) -> str:
    """Deterministic Mermaid graph of top-level units, transitively reduced.

    Units that import each other are merged, an edge implied by a longer path
    is dropped, and allowed layering exceptions are listed instead of drawn
    so they do not distort the layer order.
    """
    modules = set(graph["modules"])
    packages = {_unit(m) for m in modules if m.count(".") > 1}
    units = {_unit(m) for m in modules} - {"mdscope"}
    # An excepted import statement also records its parent packages: skip the line.
    excepted = {(s["source"], s["line"]) for s in graph["import_sites"]
                if (s["source"], s["target"]) in exceptions}
    phases: dict[tuple[str, str], set[str]] = {}
    for site in graph["import_sites"]:
        source, target = _unit(site["source"]), _unit(site["target"])
        if source != target and "mdscope" not in (source, target) and \
                (site["source"], site["line"]) not in excepted:
            phases.setdefault((source, target), set()).add(site["phase"])
    merged = {m: tuple(group) for group in cycles(units, set(phases)) for m in group}
    node = {u: merged.get(u, (u,)) for u in units}
    edges: dict[tuple[tuple, tuple], set[str]] = {}
    for (source, target), phase in phases.items():
        if node[source] != node[target]:
            edges.setdefault((node[source], node[target]), set()).update(phase)
    adjacency: dict[tuple, set[tuple]] = {}
    for source, target in edges:
        adjacency.setdefault(source, set()).add(target)

    def implied(source, target):
        stack, seen = [n for n in adjacency[source] if n != target], set()
        while stack:
            current = stack.pop()
            if current == target:
                return True
            if current not in seen:
                seen.add(current)
                stack.extend(adjacency.get(current, ()))
        return False

    def ident(group):
        return "_".join(group)

    def label(group):
        return " + ".join(f"{u}/" if u in packages else f"{u}.py" for u in group)

    lines = ["```mermaid", "flowchart TD"]
    lines += [f'    {ident(g)}["{label(g)}"]' for g in sorted(set(node.values()))]
    for (source, target), phase in sorted(edges.items()):
        if not implied(source, target):
            arrow = "-->" if "module_load" in phase else "-.->"
            lines.append(f"    {ident(source)} {arrow} {ident(target)}")
    lines.append("```")
    notes = [f"- `{label(g)}` import each other and are drawn as one node."
             for g in sorted(set(merged.values()))]
    notes += [f"- Allowed layering exception, not drawn: `{a.removeprefix('mdscope.')}` → "
              f"`{b.removeprefix('mdscope.')}`." for a, b in sorted(exceptions)]
    return "\n".join(lines + [""] + notes)


def update_docs(graph: dict, path: Path = ARCHITECTURE_DOC, *, check: bool = False) -> bool:
    """Rewrite the generated block in *path*; return whether it was out of date."""
    raw = path.read_bytes().decode("utf-8")
    newline = "\r\n" if "\r\n" in raw else "\n"
    text = raw.replace("\r\n", "\n")
    start, end = text.index(BEGIN), text.index(END) + len(END)
    block = "\n".join((BEGIN, package_dependency_mermaid(graph, layering_exceptions()), END))
    updated = text[:start] + block + text[end:]
    if updated != text and not check:
        path.write_bytes(updated.replace("\n", newline).encode("utf-8"))
    return updated != text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=REPO / "src/mdscope")
    parser.add_argument("--output", type=Path,
                        help="Directory for DOT graphs and machine-readable inventory")
    docs = parser.add_mutually_exclusive_group()
    docs.add_argument("--write-docs", action="store_true",
                      help=f"Regenerate the package diagram in {ARCHITECTURE_DOC.name}")
    docs.add_argument("--check-docs", action="store_true",
                      help="Exit 1 if the package diagram in the docs is out of date")
    args = parser.parse_args()
    root = args.package_root.resolve()
    if not (root / "__init__.py").is_file():
        parser.error(f"Not a Python package: {root}")
    if not (args.output or args.write_docs or args.check_docs):
        parser.error("Give --output, --write-docs or --check-docs")
    graph = inspect_package(root)
    if args.write_docs or args.check_docs:
        stale = update_docs(graph, check=args.check_docs)
        print(f"{ARCHITECTURE_DOC.name}: {'out of date' if stale else 'current'}"
              + (" (rewritten)" if stale and args.write_docs else ""))
        if stale and args.check_docs:
            return 1
        if args.output is None:
            return 0
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
