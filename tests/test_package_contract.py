"""Stable import and packaging contract."""
from importlib import metadata
from pathlib import Path
import subprocess
import sys

def test_version_and_public_entrypoints():
    import mdscope
    from mdscope.cli import main
    from mdscope.pipeline.runner import run_all
    from mdscope.observables.msd import fit_lateral_diffusion
    assert mdscope.__version__
    assert callable(main) and callable(run_all) and callable(fit_lateral_diffusion)

def test_console_entrypoint_metadata():
    eps = metadata.distribution("mdscope").entry_points
    assert any(ep.group == "console_scripts" and ep.name == "mdscope"
               and ep.value == "mdscope.cli:main" for ep in eps)

def test_module_help():
    result = subprocess.run([sys.executable, "-m", "mdscope.cli", "--help"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "validate" in result.stdout and "replicate" in result.stdout


"""Validate the packaged module structure, not a working-directory alias."""
from importlib import import_module
from importlib.resources import files
from pathlib import Path

def test_source_layout_and_subpackage_discovery():
    root = Path(__file__).resolve().parents[1]
    assert (root / "src/mdscope/cli.py").is_file()
    assert not (root / "analysis").exists()
    for sub in ("core", "interactions", "mechanisms", "observables",
                "pipeline", "reporting", "statistics"):
        assert import_module(f"mdscope.{sub}").__file__
    assert (files("mdscope") / "cli.py").is_file()



def _declared_requirements():
    import re
    import tomllib
    root = Path(__file__).resolve().parents[1]
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    def names(requirements):
        return {re.split(r"[<>=!~;\[ ]", r, maxsplit=1)[0].lower().replace("_", "-") for r in requirements}
    return project, names(project["dependencies"]), {
        extra: names(reqs) for extra, reqs in project["optional-dependencies"].items()}


def test_single_version_source_matches_metadata_and_citation():
    import mdscope
    project, _, _ = _declared_requirements()
    assert "version" not in project and "version" in project["dynamic"]
    citation = (Path(__file__).resolve().parents[1] / "CITATION.cff").read_text(encoding="utf-8")
    assert f"\nversion: {mdscope.__version__}\n" in citation
    assert metadata.version("mdscope") == mdscope.__version__


def test_cli_reports_version():
    import mdscope
    result = subprocess.run([sys.executable, "-m", "mdscope.cli", "--version"],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"mdscope {mdscope.__version__}"


def test_every_third_party_import_is_a_declared_dependency():
    """Transitive availability (e.g. SciPy via MDAnalysis) is not a declaration."""
    import ast
    _, required, extras = _declared_requirements()
    optional = {"freesasa": "sasa", "streamlit": "studio", "tqdm": None}  # tqdm: guarded fallback
    distributions = metadata.packages_distributions()
    src = Path(__file__).resolve().parents[1] / "src" / "mdscope"
    undeclared = {}
    for path in src.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and not node.level and node.module:
                modules = [node.module]
            else:
                continue
            for module in modules:
                top = module.split(".")[0]
                if top in sys.stdlib_module_names or top == "mdscope":
                    continue
                if top in optional:
                    extra = optional[top]
                    assert extra is None or top in extras[extra], (top, extra)
                    continue
                dists = {d.lower().replace("_", "-") for d in distributions.get(top, [top])}
                if not dists & required:
                    undeclared.setdefault(top, set()).add(str(path.relative_to(src)))
    assert not undeclared, f"Imported but not declared in [project].dependencies: {undeclared}"
