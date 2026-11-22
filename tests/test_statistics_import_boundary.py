"""Planning and lightweight statistics contracts do not eagerly import estimators."""
import subprocess
import sys


def test_statistics_package_and_planning_are_import_light():
    code = """
import sys
import mdscope.statistics
assert 'mdscope.statistics.time_series' not in sys.modules
assert 'mdscope.statistics.profiles' not in sys.modules
assert 'mdscope.statistics.interactions' not in sys.modules
import mdscope.planning
for module in ('mdscope.statistics.time_series', 'mdscope.statistics.profiles',
               'mdscope.statistics.interactions', 'mdscope.statistics.solvation'):
    assert module not in sys.modules, module
from mdscope.statistics import STATISTICS_FAMILIES
from mdscope.statistics.families import STATISTICS_FAMILIES as canonical
assert STATISTICS_FAMILIES is canonical
import importlib.util
for old in ('statistics_contracts', 'statistics_families', 'statistics_manifest',
            'statistics_selection', 'permutation_resolution', 'state_feature_spec'):
    assert importlib.util.find_spec('mdscope.' + old) is None, old
"""
    p = subprocess.run([sys.executable, "-c", code], check=False,
                       capture_output=True, text=True, timeout=45)
    assert p.returncode == 0, p.stderr


def test_dispatcher_stays_explicit():
    from mdscope.statistics.dispatch import analyze_project
    from mdscope.statistics.families import STATISTICS_FAMILIES
    assert callable(analyze_project)
    assert "time_series" in STATISTICS_FAMILIES
