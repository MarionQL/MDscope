"""Output migration, scientific display, selection and transactional publication."""
import json
from copy import deepcopy
from pathlib import Path
from runpy import run_path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from PIL import Image

from mdscope.config import DEFAULTS, figure_enabled, merge
from mdscope.config.validators.plots import validate_plots
from mdscope.contracts.plot_validation import preset_plot_ids
from mdscope.registry import figure_descriptor
from mdscope.reporting.figure_inventory import figure_generation, figure_inventory
from mdscope.reporting.overview import (
    draw_effect_intervals,
    draw_replica_distributions,
    plot_overviews,
)
from mdscope.reporting.plot_ids import figure_path, id_from_figure_path


def fixture(root):
    return run_path(str(Path(__file__).parents[1] / 'scripts/plot_preview.py'))['create_fixture'](root)


@pytest.mark.parametrize('identifier,owner,tier', [
    ('time_series:rmsd_A:timeseries', 'rmsd', 'secondary'),
    ('time_series:rmsd_A:autocorrelation', 'rmsd', 'tertiary'),
    ('profiles:lipid_order__POPC__sn1:profile', 'membrane_order', 'primary'),
    ('profiles:mystery:profile', 'unmapped', 'primary'),
    ('profiles:residue_response:heatmap', 'integrative', 'primary'),
    ('sdc4:condition:nxip_height_occlusion_joint', 'sdc4_mechanisms', 'secondary'),
    ('shared_states:state_sampling_diagnostics', 'shared_states', 'tertiary'),
    ('diffusion:diffusion_coefficients', 'membrane', 'primary'),
])
def test_catalog_paths_roundtrip_without_filename_guesses(tmp_path, identifier, owner, tier):
    meta = figure_descriptor(DEFAULTS, identifier)
    assert (meta['analysis'], meta['tier']) == (owner, tier)
    path = tmp_path / meta['destination']
    assert path.parent.parts[-2:] == (owner, tier)
    assert id_from_figure_path(path) == identifier
    assert meta['source_candidates'] == meta['prerequisites']


def test_exact_safe_identity_and_page_roundtrip(tmp_path):
    ids = ['overview:name@with/slash:effects', 'overview:name_with_slash:effects']
    paths = [figure_path({}, tmp_path/'project_summary/figures/temp.png', identifier=x) for x in ids]
    assert paths[0] != paths[1]
    assert all(path.is_relative_to(tmp_path/'project_summary/figures') for path in paths)
    for identifier, path in zip(ids, paths):
        assert id_from_figure_path(path) == identifier
        assert id_from_figure_path(path.with_stem(path.stem+'@page_002')) == identifier


def test_tiers_opt_in_exact_and_empty_authoritative():
    cfg = deepcopy(DEFAULTS)
    cfg['plots'].update(enabled=True, presets=['primary'])
    assert figure_enabled(cfg, 'overview:rmsd:distributions')
    assert not figure_enabled(cfg, 'time_series:rmsd_A:timeseries')
    cfg['plots']['include'] = ['time_series:rmsd_A:timeseries']
    assert figure_enabled(cfg, 'time_series:rmsd_A:timeseries')
    assert not figure_enabled(cfg, 'overview:rmsd:distributions')
    cfg['plots']['include'] = []
    assert not figure_enabled(cfg, 'overview:rmsd:distributions')
    ids = ['time_series:rmsd_A:autocorrelation', 'profiles:rmsf:profile']
    assert preset_plot_ids(ids[::-1], ['primary']) == ['profiles:rmsf:profile']


@pytest.mark.parametrize('presets', [['bogus'], ['primary', 'primary'], 'primary', [None]])
def test_invalid_presets(presets):
    with pytest.raises(ValueError, match='plots.presets'):
        validate_plots({'plots': {'presets': presets}})


def test_display_preserves_zero_missing_negative_and_invalid_interval():
    fig, ax = plt.subplots()
    rows = [dict(label=str(i), effect=value, low=lo, high=hi) for i, (value,lo,hi) in enumerate([
        (0., 0., 0.), (-2., -3., -1.), (1., 3., 2.), (np.nan, np.nan, np.nan)])]
    draw_effect_intervals(ax, rows, unit='nm²')
    assert len(ax.collections) == 3  # every finite effect, including zero
    assert sum(text.get_text() == 'CI unavailable' for text in ax.texts) == 2
    assert len(ax.get_yticklabels()) == 4
    assert 'nm²' in ax.get_xlabel()
    plt.close(fig)
    fig, ax = plt.subplots()
    draw_replica_distributions(ax, [('measured zero', [0.,0.]), ('missing', []), ('negative', [-1.,-2.])], unit='Å')
    assert len(ax.get_xticklabels()) == 3
    assert any(text.get_text() == 'Missing' for text in ax.texts)
    plt.close(fig)


def test_publication_failure_retry_and_unselected_files(tmp_path):
    cfg = merge(DEFAULTS, {'plots': {'enabled': True, 'include': ['time_series:rmsd_A:timeseries']}})
    old = tmp_path / 'project_summary/figures/time_series/rmsd_A__timeseries.png'
    old.parent.mkdir(parents=True)
    Image.new('RGB', (12, 12), 'red').save(old)
    unrelated = old.parent / 'user_notes.txt'
    unrelated.write_text('keep me')
    original = old.read_bytes()
    with pytest.raises(RuntimeError, match='injected'):
        with figure_generation(cfg, tmp_path, {}) as staging:
            target = figure_path(cfg, staging/'project_summary/figures/time_series/rmsd_A__timeseries.png')
            Image.new('RGB', (12,12), 'blue').save(target)
            raise RuntimeError('injected')
    assert old.read_bytes() == original
    assert unrelated.read_text() == 'keep me'
    with figure_generation(cfg, tmp_path, {}) as staging:
        target = figure_path(cfg, staging/'project_summary/figures/time_series/rmsd_A__timeseries.png')
        Image.new('RGB', (12,12), 'blue').save(target)
    canonical = tmp_path / figure_descriptor(cfg, 'time_series:rmsd_A:timeseries')['destination']
    assert canonical.is_file() and not old.exists()
    assert unrelated.read_text() == 'keep me'
    records = json.loads((tmp_path/'project_summary/figures/figure_provenance.json').read_text())
    assert records[canonical.relative_to(tmp_path).as_posix()]['plot_id'] == 'time_series:rmsd_A:timeseries'
    inventory = figure_inventory(tmp_path, {}, cfg, include_rejected=True)
    assert inventory.iloc[0]['analysis'] == 'rmsd'
    assert inventory.iloc[0]['tier'] == 'secondary'


def test_saved_overview_uneven_replicas_mixed_units_and_no_numeric_writes(tmp_path, monkeypatch):
    cfg, reps = fixture(tmp_path)
    cfg = merge(DEFAULTS, cfg)
    cfg['plots'].update(include=['overview:rmsd:distributions', 'overview:membrane:effects'])
    cfg['analyses']['membrane']['enabled'] = True
    comparisons = tmp_path/'project_summary/statistics/time_series/system_comparisons.parquet'
    pd.DataFrame([dict(metric='membrane_thickness_A', group=None, system_1='condition_A', system_2='condition_B',
                       mean_difference_system2_minus_system1=-1., mean_difference_ci_low=-2., mean_difference_ci_high=0.)]).to_parquet(comparisons)
    paths = list(tmp_path.rglob('*.parquet'))
    before = {path: path.read_bytes() for path in paths}
    import MDAnalysis as mda
    monkeypatch.setattr(mda, 'Universe', lambda *a, **k: pytest.fail('No topology/trajectory reads'))
    plot_overviews(cfg, tmp_path, reps)
    figures = list((tmp_path/'project_summary/figures').rglob('*.png'))
    assert {id_from_figure_path(path) for path in figures} == set(cfg['plots']['include'])
    assert all(path.read_bytes() == before[path] for path in paths)
    cfg['plots']['include'] = []
    monkeypatch.setattr(pd, 'read_parquet', lambda *a, **k: pytest.fail('Empty selection must read nothing'))
    plot_overviews(cfg, tmp_path, reps)


def test_large_production_heatmap_pages_cover_all_species(tmp_path):
    from mdscope.reporting.lipid_contact_bootstrap import _draw
    from mdscope.reporting.templates import heatmap_pages
    labels = [f'species_{i:03d}' for i in range(61)]
    values = np.arange(61*61, dtype=float).reshape(61,61)
    values[0,0], values[0,1] = np.nan, 0.
    coverage = np.zeros(values.shape, dtype=int)
    for fig, ax, rs, cs, _ in heatmap_pages(values, labels, labels, title='Full coverage', units='ratio', vmin=0, vmax=4000):
        coverage[rs,cs] += 1
        assert ax.images[0].get_clim() == (0,4000)
        plt.close(fig)
    assert np.all(coverage == 1)
    table = pd.DataFrame([dict(species_a=a, species_b=b, mean=0.) for a in labels for b in labels])
    identifier = 'lipid_contacts:lipid_contact_matrix_generic_upper_enrichment'
    cfg = {'plots': {'enabled': True, 'include': [identifier], 'dpi': 50}}
    _draw(cfg, tmp_path/'project_summary/figures/lipid_contacts', table, ordered=labels,
          value='mean', supported=False, title='Generic saved matrix', stem=identifier.split(':')[1],
          metric='enrichment', limit=1., comparison=False)
    outputs = list((tmp_path/'project_summary/figures').rglob('*.png'))
    assert len(outputs) == 6
    assert all(id_from_figure_path(path) == identifier for path in outputs)


def test_shrinking_selected_overview_retires_only_its_old_pages(tmp_path):
    groups = {f'chain_{i:02d}': {'selection': 'protein', 'label': f'Protein {i}'} for i in range(12)}
    cfg = merge(DEFAULTS, {'selections': {'groups': groups},
                          'plots': {'enabled': True, 'include': ['overview:rmsd:distributions'], 'dpi': 50}})
    rep = tmp_path/'generic/rep1'
    summary = rep/'summary'
    summary.mkdir(parents=True)
    pd.DataFrame({'time_ps': [0.,1.,2.], 'rmsd_A': [0.,1.,2.]}).to_parquet(summary/'timeseries.parquet')
    pd.DataFrame([dict(group=group, time_ps=float(i), rmsd_A=float(i)) for group in groups for i in range(3)]).to_parquet(summary/'group_timeseries.parquet')
    with figure_generation(cfg, tmp_path, {'generic': [rep]}) as staging:
        plot_overviews(cfg, staging, {'generic': [rep]})
    folder = tmp_path/'project_summary/figures/rmsd/primary'
    assert len(list(folder.glob('*.png'))) == 4
    unrelated = folder/'user-note.txt'
    unrelated.write_text('keep')
    cfg['selections']['groups'] = {'chain_00': groups['chain_00']}
    with figure_generation(cfg, tmp_path, {'generic': [rep]}) as staging:
        plot_overviews(cfg, staging, {'generic': [rep]})
    assert len(list(folder.glob('*.png'))) == 1
    assert unrelated.read_text() == 'keep'
    assert len(list((tmp_path/'.archives/figures').rglob('*@page_004.png'))) == 1


def test_canonical_interaction_inventory_and_residue_settings(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    cfg = merge(DEFAULTS, {'plots': {'enabled': True}})
    identifier = 'interactions:salt_bridges:occupancy_heatmap'
    path = tmp_path/figure_descriptor(cfg, identifier)['destination']
    path.parent.mkdir(parents=True)
    Image.new('RGB',(10,10)).save(path)
    assert figure_status(path, cfg, tmp_path, {}) == 'legacy_unverified'
    cfg['statistics']['residue_response']['enabled'] = False
    path = tmp_path/figure_descriptor(cfg, 'profiles:residue_response:heatmap')['destination']
    assert figure_status(path, cfg, tmp_path, {}) == 'disabled'


def test_presentation_text_roundtrips_and_never_changes_values():
    from mdscope.contracts.plot_validation import validate_figure_options
    from mdscope.reporting.plot_ids import apply_figure_options
    options = {'title': 'Generic protein response', 'annotations': ['Independent replicas; saved intervals.'], 'dpi': 100}
    validate_figure_options({'overview:rmsd:effects': options})
    fig, ax = plt.subplots()
    data = np.array([-2., 0., 3.])
    ax.plot(data)
    apply_figure_options(fig, options)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), data)
    assert fig._suptitle.get_text() == options['title']
    assert fig.texts[-1].get_text() == options['annotations'][0]
    plt.close(fig)


def test_preset_expansion_is_deterministic_and_explains_unavailable_statistics(tmp_path):
    from mdscope.planning import resolve_figure_presets
    cfg, reps = fixture(tmp_path)
    cfg = merge(DEFAULTS, cfg)
    cfg['analyses']['rmsd']['enabled'] = True
    cfg['plots'].update(presets=['primary'])
    for spec in cfg['systems'].values():
        spec['topology'] = 'unopened.pdb'
        for replica in spec['replicates']:
            replica['trajectory'] = 'unopened.xtc'
    one, omitted = resolve_figure_presets(cfg, tmp_path, tmp_path/'config.yaml')
    two, _ = resolve_figure_presets(cfg, tmp_path, tmp_path/'config.yaml')
    assert one['plots']['include'] == sorted(one['plots']['include']) == two['plots']['include']
    assert 'overview:rmsd:distributions' in one['plots']['include']
    assert 'time_series:rmsd_A:timeseries' not in one['plots']['include']
    assert 'overview:rmsd:effects' in one['plots']['include']
    cfg['statistics']['enabled'] = False
    disabled, omitted = resolve_figure_presets(cfg, tmp_path, tmp_path/'config.yaml')
    assert 'overview:rmsd:distributions' in disabled['plots']['include']
    assert not any(figure_descriptor(cfg, identifier)['requires_statistics'] for identifier in disabled['plots']['include'])
    cfg['plots']['include'] = []
    assert resolve_figure_presets(cfg, tmp_path, tmp_path/'config.yaml') == (cfg, [])


def test_canonical_sources_and_presentation_options_invalidate_provenance(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    cfg = merge(DEFAULTS, {'plots': {'enabled': True, 'include': ['overview:rmsd:distributions']}})
    rep = tmp_path/'generic/rep1'
    source = rep/'summary/timeseries.parquet'
    source.parent.mkdir(parents=True)
    pd.DataFrame({'time_ps': [0.,1.], 'rmsd_A': [0.,1.]}).to_parquet(source)
    with figure_generation(cfg, tmp_path, {'generic': [rep]}) as staging:
        plot_overviews(cfg, staging, {'generic': [rep]})
    image = next((tmp_path/'project_summary/figures').rglob('*.png'))
    assert figure_status(image, cfg, tmp_path, {'generic': [rep]}) == 'valid'
    cfg['plots']['options'] = {'overview:rmsd:distributions': {'title': 'Changed'}}
    assert figure_status(image, cfg, tmp_path, {'generic': [rep]}) == 'stale_configuration'
    cfg['plots']['options'] = {}
    pd.DataFrame({'time_ps': [0.,1.], 'rmsd_A': [2.,3.]}).to_parquet(source)
    assert figure_status(image, cfg, tmp_path, {'generic': [rep]}) == 'stale_source'


def test_namespace_unit_notes_do_not_mix_incompatible_effects(tmp_path, monkeypatch):
    from mdscope.reporting import overview
    folder = tmp_path/'project_summary/statistics/time_series'
    folder.mkdir(parents=True)
    pd.DataFrame([dict(metric=metric, system_1='A', system_2='B',
                       mean_difference_system2_minus_system1=1., mean_difference_ci_low=-1.,
                       mean_difference_ci_high=2., ci_level=.95)
                  for metric in ('sdc4_tm_crossing_angle_deg', 'sdc4_local_thickness_A',
                                 'sdc4_local_upper_curvature_Ainv')]).to_parquet(folder/'system_comparisons.parquet')
    cfg = merge(DEFAULTS, {'analyses': {'sdc4_mechanisms': {'enabled': True}},
                          'plots': {'enabled': True, 'include': ['overview:sdc4_mechanisms:effects']}})
    captured = []
    original = overview.draw_effect_intervals
    def inspect(ax, rows, *, unit):
        captured.append((len(rows), unit))
        return original(ax, rows, unit=unit)
    monkeypatch.setattr(overview, 'draw_effect_intervals', inspect)
    overview.plot_overviews(cfg, tmp_path, {})
    assert sorted(captured) == [(1, '1/Å'), (1, 'degrees'), (1, 'Å')]
    assert len(list((tmp_path/'project_summary/figures').rglob('*.png'))) == 3


def test_coverage_counts_consolidated_pdf_sources_as_rendered(tmp_path):
    from mdscope.reporting.plotting import plot_all
    cfg, reps = fixture(tmp_path)
    identifier = 'profiles:lipid_order__POPC__sn1:comparison_table'
    cfg['plots'].update(include=[identifier], table_reports_pdf=True)
    with figure_generation(cfg, tmp_path, reps) as staging:
        plot_all(cfg, staging, reps)
    coverage = json.loads((tmp_path/'project_summary/figures/figure_coverage.json').read_text())
    assert coverage['figures'][0]['id'] == identifier
    assert coverage['figures'][0]['status'] == 'rendered'
    assert coverage['figures'][0]['outputs'] == ['project_summary/figures/membrane_order/secondary/statistics_tables.pdf']


def test_missing_effect_group_is_labeled_global_across_pandas_versions(tmp_path, monkeypatch):
    from mdscope.reporting import overview
    folder = tmp_path/'project_summary/statistics/time_series'
    folder.mkdir(parents=True)
    pd.DataFrame([dict(metric='rmsd_A', group=np.nan, system_1='A', system_2='B',
                       mean_difference_system2_minus_system1=0., mean_difference_ci_low=-1.,
                       mean_difference_ci_high=1.)]).to_parquet(folder/'system_comparisons.parquet')
    cfg = merge(DEFAULTS, {'plots': {'enabled': True, 'include': ['overview:rmsd:effects']}})
    labels = []
    original = overview.draw_effect_intervals
    def inspect(ax, rows, *, unit):
        labels.extend(row['label'] for row in rows)
        return original(ax, rows, unit=unit)
    monkeypatch.setattr(overview, 'draw_effect_intervals', inspect)
    overview.plot_overviews(cfg, tmp_path, {})
    assert labels == ['rmsd_A · global']
