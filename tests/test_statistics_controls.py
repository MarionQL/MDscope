"""Family isolation, gapped-time safeguards and finite-replica reference answers."""
from copy import deepcopy
from itertools import product
import numpy as np
import pandas as pd
import pytest
from mdscope.statistics import analyze_project
from mdscope.statistics.common import _block_bootstrap, _permutation
from mdscope.statistics.time_series import analyze_time_series
from mdscope.config import DEFAULTS, validate_config, statistics_family_config


def test_nonconstant_moving_blocks_enumerated_distribution():
    class Enumeration:
        def integers(self, low, high, size, dtype):
            assert (low, high, size) == (0, 2, (4, 2))
            return np.array(list(product([0, 1], repeat=2)), dtype=dtype)
    # Four equally likely concatenations of [0,2] or [2,8], truncated to 3.
    actual = _block_bootstrap([0., 2., 8.], 2, 4, Enumeration())
    np.testing.assert_allclose(actual, [2/3, 4/3, 10/3, 4.])


def test_three_vs_three_independent_permutation_resolution():
    cfg = {"statistics": {"hypothesis_tests": {"enabled": True}}}
    # Perfect separation: two of C(6,3)=20 allocations are as extreme.
    assert _permutation(dict(enumerate([0., 1., 2.])), dict(enumerate([8., 9., 10.])), cfg, np.random.default_rng(4)) == .1


def test_disabled_statistics_families_never_enter_kernels(tmp_path, monkeypatch):
    import mdscope.statistics as stats
    import mdscope.statistics.time_series as ts
    import mdscope.statistics.profiles as profiles
    import mdscope.statistics.interactions as interactions
    import mdscope.statistics.solvation as solvation
    def forbidden(*a, **k):
        raise AssertionError("disabled family entered data collection")
    monkeypatch.setattr(ts, 'collect_time_series', forbidden)
    monkeypatch.setattr(profiles, 'collect_temporal_profiles', forbidden)
    monkeypatch.setattr(interactions, 'interaction_tables', forbidden)
    monkeypatch.setattr(solvation.pd, 'read_parquet', forbidden)
    cfg = {'statistics': {'families': {k: False for k in ('time_series','profiles','interactions','solvation_rdfs')}}}
    analyze_project(cfg, tmp_path, {'x': [tmp_path/'rep1']})
    import json
    manifest = json.loads((tmp_path/'project_summary/statistics/manifest.json').read_text())
    assert set(manifest['families'].values()) == {'disabled'}


def test_family_overrides_do_not_mutate_or_change_default_config():
    cfg = {'statistics': {'bootstrap_draws': 5000, 'families': {'interactions': {'bootstrap_draws': 100, 'hypothesis_tests': False}}}}
    before = deepcopy(cfg)
    assert statistics_family_config(cfg, 'time_series') is cfg
    selected = statistics_family_config(cfg, 'interactions')
    assert selected['statistics']['bootstrap_draws'] == 100
    assert selected['statistics']['hypothesis_tests']['enabled'] is False
    assert cfg == before


def test_global_bootstrap_disable_can_be_overridden_by_one_family():
    cfg = {"statistics": {"bootstrap":False,"bootstrap_draws":100,
                          "families":{"profiles":{"bootstrap":True,"hierarchy":{"enabled":False}}}}}
    before = deepcopy(cfg)
    assert statistics_family_config(cfg,"time_series")["statistics"]["bootstrap_draws"] == 0
    selected = statistics_family_config(cfg,"profiles")["statistics"]
    assert selected["bootstrap_draws"] == 100 and selected["bootstrap"]
    assert selected["hierarchy"]["enabled"] is False
    assert cfg == before


@pytest.mark.parametrize('bad', [{'unknown': False}, {'profiles': {'enabled': 'no'}}, {'time_series': {'bootstrap_draws': 0}}, {'profiles': {'ci': 2}}, {'profiles': {'typo': 1}}])
def test_invalid_statistics_controls_fail_before_processing(tmp_path, bad):
    cfg = deepcopy(DEFAULTS)
    top = tmp_path/'a.pdb'; top.touch()
    traj = tmp_path/'a.xtc'; traj.touch()
    cfg['systems'] = {'x': {'topology': str(top), 'replicates': [str(traj)]}}
    cfg['statistics']['families'] = bad
    with pytest.raises(ValueError, match='statistics.families'):
        validate_config(cfg)


def test_gapped_scalar_series_keeps_means_but_not_temporal_intervals(tmp_path):
    paths = {'x': []}
    for i in range(2):
        rep = tmp_path/'x'/f'rep{i}'
        (rep/'summary').mkdir(parents=True)
        pd.DataFrame({'time_ps': np.arange(6)*10., 'value_A': [0., 1., np.nan, 3., 4., 5.]}).to_parquet(rep/'summary/timeseries.parquet')
        paths['x'].append(rep)
    cfg = deepcopy(DEFAULTS); cfg['statistics']['bootstrap_draws'] = 100
    analyze_time_series(cfg, tmp_path, paths)
    rep = pd.read_parquet(tmp_path/'project_summary/statistics/time_series/replicate_statistics.parquet')
    system = pd.read_parquet(tmp_path/'project_summary/statistics/time_series/system_statistics.parquet')
    assert rep['mean'].eq(2.6).all()
    assert rep['effective_n'].isna().all()
    assert rep.mean_ci_low.isna().all() and system.mean_ci_low.isna().all()


@pytest.mark.parametrize('field,value,diagnostic', [
    (('trajectory','stride'),0,'stride'), (('trajectory','chunk_frames'),0,'chunk_frames'),
    (('time','dt_ps'),0.,'dt_ps'), (('time','dt_ps'),np.nan,'dt_ps'),
    (('time','origin_ps'),np.inf,'origin_ps'),
    (('statistics','ci'),1.2,'ci'), (('statistics','bootstrap_draws'),0,'bootstrap_draws'),
])
def test_invalid_numeric_configuration_fails_before_trajectory_processing(tmp_path,field,value,diagnostic):
    cfg=deepcopy(DEFAULTS)
    top=tmp_path/'a.pdb';top.touch();trj=tmp_path/'a.xtc';trj.touch()
    cfg['systems']={'x':{'topology':str(top),'replicates':[str(trj)]}}
    cfg[field[0]][field[1]]=value
    with pytest.raises(ValueError,match=diagnostic):validate_config(cfg)


def test_empty_forward_window_is_rejected_early(tmp_path):
    cfg=deepcopy(DEFAULTS);cfg['systems']={'x':{}}
    cfg['trajectory'].update(start=10,stop=10)
    with pytest.raises(ValueError,match='forward frame window'):validate_config(cfg)


@pytest.mark.parametrize('design', ['independent', 'paired'])
def test_one_invalid_replica_cannot_be_dropped_from_contrast_intervals(tmp_path, design):
    paths = {'a': [], 'b': []}
    for system in paths:
        for i in range(3):
            rep = tmp_path/system/f'rep{i}'
            (rep/'summary').mkdir(parents=True)
            values = np.arange(6, dtype=float) + i + (10 if system == 'b' else 0)
            if system == 'b' and i == 1:
                values[2] = np.nan
            pd.DataFrame({'time_ps': np.arange(6)*10., 'value_A': values}).to_parquet(rep/'summary/timeseries.parquet')
            paths[system].append(rep)
    cfg = deepcopy(DEFAULTS)
    cfg['statistics']['bootstrap_draws'] = 100
    cfg['statistics']['hypothesis_tests']['design'] = design
    analyze_time_series(cfg, tmp_path, paths)
    output = pd.read_parquet(tmp_path/'project_summary/statistics/time_series/system_comparisons.parquet')
    assert len(output) == 1
    for prefix in ('mean', 'median'):
        assert output[f'{prefix}_difference_ci_low'].isna().all()
        assert output[f'{prefix}_difference_ci_high'].isna().all()
        assert output[f'{prefix}_probability_difference_gt_0'].isna().all()
    # Finite-observation descriptive means retain the gapped replica equally.
    assert output.mean_difference_system2_minus_system1.item() == pytest.approx(10 + .1/3)
