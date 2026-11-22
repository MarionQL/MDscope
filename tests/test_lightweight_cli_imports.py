"""Read-only commands must not eagerly load execution or statistics machinery."""

import subprocess
import sys

from mdscope.statistics.families import STATISTICS_FAMILIES


def test_cli_import_keeps_execution_modules_unloaded():
    code = (
        "import sys\n"
        "import mdscope.cli\n"
        "for module in ('mdscope.pipeline.runner', 'mdscope.pipeline.streaming',\n"
        "               'mdscope.reporting.report', 'mdscope.statistics.dispatch',\n"
        "               'mdscope.pipeline.preflight'):\n"
        "    assert module not in sys.modules, module\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, timeout=60,
                            check=False)
    assert result.returncode == 0, result.stderr


def test_planner_import_does_not_load_statistical_estimators():
    code = (
        "import sys\n"
        "from mdscope.planning import build_plan\n"
        "from mdscope.statistics.families import STATISTICS_FAMILIES\n"
        "assert STATISTICS_FAMILIES\n"
        "for module in ('mdscope.statistics.dispatch', 'mdscope.statistics.profiles',\n"
        "               'mdscope.statistics.time_series', 'mdscope.pipeline.runner',\n"
        "               'mdscope.pipeline.fused_membrane'):\n"
        "    assert module not in sys.modules, module\n"
    )
    result = subprocess.run([sys.executable, "-c", code],
                            capture_output=True, text=True, timeout=60,
                            check=False)
    assert result.returncode == 0, result.stderr


def test_statistics_dispatcher_uses_same_family_identifiers():
    from mdscope.statistics import STATISTICS_FAMILIES as EXECUTION_FAMILIES

    assert EXECUTION_FAMILIES is STATISTICS_FAMILIES
    assert STATISTICS_FAMILIES == (
        "time_series", "profiles", "interactions", "solvation_rdfs",
        "diffusion", "coupling", "lipid_contacts",
    )
