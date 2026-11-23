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
from mdscope.contracts.statistics import STATISTICS_FAMILIES
assert STATISTICS_FAMILIES
import importlib.util
for old in ('statistics_contracts', 'statistics_families', 'statistics_manifest',
            'statistics_selection', 'permutation_resolution', 'state_feature_spec',
            'output_paths', 'biophysical_common', 'statistics.contracts', 'contracts.metric_selection'):
    assert importlib.util.find_spec('mdscope.' + old) is None, old
"""
    p = subprocess.run([sys.executable, "-c", code], check=False,
                       capture_output=True, text=True, timeout=45)
    assert p.returncode == 0, p.stderr


def test_dispatcher_stays_explicit():
    from mdscope.contracts.statistics import STATISTICS_FAMILIES
    from mdscope.statistics.dispatch import analyze_project
    assert callable(analyze_project)
    assert "time_series" in STATISTICS_FAMILIES


def test_dispatch_loads_only_enabled_estimator_families(tmp_path):
    code = '''
import sys
from copy import deepcopy
from pathlib import Path
from mdscope.config import DEFAULTS
from mdscope.contracts.statistics import STATISTICS_FAMILIES
from mdscope.statistics.dispatch import analyze_project
modules = ('time_series', 'profiles', 'interactions', 'solvation', 'diffusion', 'coupling', 'lipid_contacts')
assert not any('mdscope.statistics.' + name in sys.modules for name in modules)
cfg = deepcopy(DEFAULTS)
cfg['progress']['enabled'] = False
cfg['statistics']['families'] = {name: name == 'time_series' for name in STATISTICS_FAMILIES}
analyze_project(cfg, Path(sys.argv[1]), {})
assert 'mdscope.statistics.time_series' in sys.modules
assert not any('mdscope.statistics.' + name in sys.modules for name in modules[1:])
'''
    result = subprocess.run([sys.executable, '-c', code, str(tmp_path)],
                            capture_output=True, text=True, timeout=45, check=False)
    assert result.returncode == 0, result.stderr
