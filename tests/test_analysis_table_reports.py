"""Saved-data reports retain populations, exact selections and publication integrity."""
import json

import numpy as np
import pandas as pd
import pytest
from pypdf import PdfReader

from mdscope.reporting.table_reports import table_report_session


def fixture(root):
    folder = root / 'project_summary/statistics/profiles'
    folder.mkdir(parents=True)
    (folder.parent / 'manifest.json').write_text(json.dumps({
        'status': 'complete', 'methods': {}, 'families': {'profiles': 'enabled'}}))
    rows, conditions = [], []
    for species, tails in [('POPC', ['sn1', 'sn2']), ('PSM', ['n_acyl', 'sphingosine'])]:
        for tail in tails:
            for leaf in ['', '__lower', '__upper']:
                profile = f'lipid_order__{species}__{tail}{leaf}'
                for position in range(45):
                    for a, b in [('A', 'B'), ('A', 'C')]:
                        rows.append(dict(profile=profile, position=position, position_label=str(position),
                                         system_1=a, system_2=b, difference_system2_minus_system1=.03,
                                         difference_ci_low=1e-12 if position else 0.,
                                         difference_ci_high=.1, probability_difference_gt_0=.98,
                                         p_value=.01, q_value=np.nan if position == 2 else .04))
                    for system, n in [('A', 3), ('B', 2), ('C', 1)]:
                        conditions.append(dict(profile=profile, position=position, position_label=str(position),
                                               system=system, mean=.2, mean_ci_low=.1,
                                               mean_ci_high=.3, n_replicates=n))
    pd.DataFrame(rows).to_parquet(folder / 'system_comparisons.parquet')
    pd.DataFrame(conditions).to_parquet(folder / 'system_profiles.parquet')
    cfg = {'plots': {'enabled': True, 'table_reports_pdf': True, 'bootstrap_tables': False},
           'statistics': {'enabled': True}, 'analyses': {}}
    return cfg, folder


def test_consolidated_report_preserves_all_rows_and_full_precision(tmp_path):
    from mdscope.reporting.statistical_tables import render_statistical_reports
    cfg, source = fixture(tmp_path)
    original = {p: p.read_bytes() for p in source.iterdir()}
    with table_report_session(tmp_path, enabled=True, cfg=cfg):
        render_statistical_reports(cfg, tmp_path)
    reports = list((tmp_path / 'project_summary/figures').rglob('statistics_tables.pdf'))
    assert len(reports) == 1
    report = reports[0]
    saved = pd.read_csv(report.with_suffix('.csv'))
    assert len(saved) == 1080
    assert set(saved['tail']) == {'sn1', 'sn2', 'n_acyl', 'sphingosine'}
    assert set(saved.leaflet) == {'lower', 'upper', 'combined'}
    assert saved['condition_1_mean'].eq(.2).all()
    assert set(saved.condition_2_n_replicates) == {1, 2}
    assert saved.loc[saved.position == 1, 'difference_ci_low'].eq(1e-12).all()
    assert not saved.loc[saved.position == 0, 'ci_excludes_zero'].any()
    assert saved.loc[saved.position == 1, 'ci_excludes_zero'].all()
    assert not saved.loc[saved.position == 2, 'q_below_alpha'].any()
    pages = PdfReader(report).pages
    assert all(float(page.mediabox.width) == pytest.approx(16.54 * 72) for page in pages)
    assert all(float(page.mediabox.height) == pytest.approx(11.69 * 72) for page in pages)
    text = '\n'.join(page.extract_text() for page in pages)
    assert 'Lower sn1' in text and 'Combined sn2' in text
    assert 'sphingosine' in text and 'n_acyl' in text
    assert '44' in text
    assert {p: p.read_bytes() for p in source.iterdir()} == original
    assert not list((tmp_path / 'project_summary/figures').rglob('*.png'))


def test_exact_source_selection_does_not_expand_to_other_metrics(tmp_path):
    from mdscope.reporting.statistical_tables import render_statistical_reports
    cfg, _ = fixture(tmp_path)
    chosen = 'profiles:lipid_order__POPC__sn1__lower:comparison_table'
    cfg['plots']['include'] = [chosen]
    with table_report_session(tmp_path, enabled=True, cfg=cfg):
        render_statistical_reports(cfg, tmp_path)
    report = next((tmp_path / 'project_summary/figures').rglob('statistics_tables.pdf'))
    saved = pd.read_csv(report.with_suffix('.csv'))
    assert saved.profile.unique().tolist() == ['lipid_order__POPC__sn1__lower']
    assert len(saved) == 90
    manifest = json.loads(report.with_name('statistics_tables.manifest.json').read_text())
    assert manifest['source_ids'] == [chosen]


def small_fixture(root):
    folder = root / 'project_summary/statistics/time_series'
    folder.mkdir(parents=True)
    (folder.parent / 'manifest.json').write_text(json.dumps({'status': 'complete', 'methods': {}}))
    rows = []
    for group in ['chain_A', 'chain_B']:
        rows.append(dict(metric='rmsd_A', group=group, system_1='A', system_2='B',
                         mean_difference_system2_minus_system1=.1, mean_difference_ci_low=0.,
                         mean_difference_ci_high=.2, median_difference_ci_low=-1e-12,
                         median_difference_ci_high=-1e-14, mean_q_value=np.nan, median_q_value=.01))
    pd.DataFrame(rows).to_parquet(folder / 'system_comparisons.parquet')
    pd.DataFrame([dict(metric='rmsd_A', group=g, system=s, mean=v, mean_ci_low=v-.1,
                       mean_ci_high=v+.1, n_replicates=n)
                  for g, v in [('chain_A', 1.), ('chain_B', 2.)]
                  for s, n in [('A', 3), ('B', 2)]]).to_parquet(folder / 'system_statistics.parquet')
    cfg = {'plots': {'enabled': True, 'bootstrap_tables': False, 'table_reports_pdf': True}}
    return cfg, folder


def test_named_groups_median_estimand_and_unknown_variables(tmp_path):
    from mdscope.reporting.statistical_tables import collect_statistical_reports
    from mdscope.registry import statistical_variable_owner
    cfg, _ = small_fixture(tmp_path)
    rows = collect_statistical_reports(cfg, tmp_path)
    assert rows.condition_1_mean.tolist() == [1., 2.]
    assert rows.analysis.tolist() == ['rmsd', 'rmsd']
    assert not rows.ci_excludes_zero.any()
    assert rows.median_difference_ci_excludes_zero.all()
    assert rows.median_q_value_below_alpha.all()
    assert statistical_variable_owner(cfg, 'time_series', {'metric': 'surprise', 'group': 'solvation/water/protein/all'})[0] == 'unmapped'
    assert statistical_variable_owner(cfg, 'interactions', {'metric': 'surprise', 'interaction_type': 'salt_bridges'})[0] == 'unmapped'
    runs = {'analyses': {'sasa': {'enabled': True, 'runs': {'environment': {'enabled': True}}}}}
    assert statistical_variable_owner(runs, 'profiles', {'profile': 'sasa_mean__environment'}) == ('sasa', 'nm^2')


def test_consolidated_csv_failure_rolls_back_and_freshness_is_checked(tmp_path, monkeypatch):
    import pytest
    from mdscope.reporting import table_reports
    from mdscope.reporting.figure_inventory import figure_status
    from mdscope.reporting.statistical_tables import render_statistical_reports
    cfg, folder = small_fixture(tmp_path)
    cfg['plots']['include'] = ['time_series:rmsd_A__group_chain_A:comparison_table']
    with table_report_session(tmp_path, enabled=True, cfg=cfg):
        render_statistical_reports(cfg, tmp_path)
    report = next((tmp_path / 'project_summary/figures').rglob('statistics_tables.pdf'))
    assert figure_status(report, cfg, tmp_path, {}) == 'legacy_unverified'
    files = [report, report.with_suffix('.csv'), report.with_name('statistics_tables.manifest.json')]
    before = {p: p.read_bytes() for p in files}
    replace = table_reports.os.replace
    def fail_csv(source, destination):
        if str(source).endswith('.pending.csv'):
            raise OSError('CSV publication failure')
        return replace(source, destination)
    with monkeypatch.context() as patch:
        patch.setattr(table_reports.os, 'replace', fail_csv)
        with pytest.raises(OSError, match='CSV publication'):
            with table_report_session(tmp_path, enabled=True, cfg=cfg):
                render_statistical_reports(cfg, tmp_path)
    assert {p: p.read_bytes() for p in files} == before
    cfg['plots']['include'].append('time_series:rmsd_A__group_chain_B:comparison_table')
    assert figure_status(report, cfg, tmp_path, {}) == 'stale_configuration'
    cfg['plots']['include'].pop()
    saved = pd.read_parquet(folder / 'system_statistics.parquet')
    saved.loc[0, 'mean'] = 9.
    saved.to_parquet(folder / 'system_statistics.parquet')
    assert figure_status(report, cfg, tmp_path, {}) == 'stale_source'
    (folder / 'system_statistics.parquet').unlink()
    assert figure_status(report, cfg, tmp_path, {}) == 'missing_source'


def test_malformed_selected_statistics_preserve_reports(tmp_path):
    import pytest
    from mdscope.reporting.statistical_tables import render_statistical_reports
    cfg, folder = small_fixture(tmp_path)
    with table_report_session(tmp_path, enabled=True, cfg=cfg):
        render_statistical_reports(cfg, tmp_path)
    report = next((tmp_path / 'project_summary/figures').rglob('statistics_tables.pdf'))
    before = {p: p.read_bytes() for p in report.parent.iterdir()}
    pd.DataFrame({'garbage': [1]}).to_parquet(folder / 'system_comparisons.parquet')
    with pytest.raises(ValueError, match='missing identity columns'):
        with table_report_session(tmp_path, enabled=True, cfg=cfg):
            render_statistical_reports(cfg, tmp_path)
    assert {p: p.read_bytes() for p in report.parent.iterdir()} == before


def test_interval_classification_uses_finite_ordered_bounds_and_configured_alpha(tmp_path):
    from mdscope.reporting.statistical_tables import collect_statistical_reports
    cfg, folder = small_fixture(tmp_path)
    cfg['statistics'] = {'hypothesis_tests': {'alpha': .02}}
    base = pd.read_parquet(folder / 'system_comparisons.parquet').iloc[0].to_dict()
    bounds = [(-1., -.5), (0., .1), (1e-13, .1), (np.nan, 1.), (1., np.inf), (2., 1.), (-.1, 0.)]
    rows = [{**base, 'group': f'case_{i}', 'mean_difference_ci_low': lo,
             'mean_difference_ci_high': hi, 'mean_q_value': .03}
            for i, (lo, hi) in enumerate(bounds)]
    pd.DataFrame(rows).to_parquet(folder / 'system_comparisons.parquet')
    data = collect_statistical_reports(cfg, tmp_path)
    assert data.ci_excludes_zero.tolist() == [True, False, True, False, False, False, False]
    assert not data.q_below_alpha.any()
    assert data.condition_1_mean.isna().all()  # unsupported is never filled with zero
