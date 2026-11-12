"""Offline family selection and adversarial stale-figure integrity checks."""
from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
from PIL import Image
from mdscope.reporting.figure_inventory import (
    record_generated_figures,figure_inventory,figure_status)
from mdscope.reporting.plotting import plot_interactions, _metric_meta


def _fixture(tmp_path):
    rep=tmp_path/'condition/rep1'; (rep/'summary').mkdir(parents=True)
    pd.DataFrame({'time_ps':[0.,10.,20.],'rmsd_A':[1.,2.,3.]}).to_parquet(rep/'summary/timeseries.parquet')
    fig=tmp_path/'project_summary/figures/time_series/rmsd_A__timeseries.png'
    fig.parent.mkdir(parents=True);Image.new('RGB',(50,50),'blue').save(fig)
    cfg={'project':{'name':'qa'},'systems':{'condition':{'replicates':[{'label':'rep1'}]}},
         'plots':{'enabled':True},'statistics':{'enabled':True}}
    return cfg,{'condition':[rep]},fig


def test_readable_older_figures_survive_render_settings_changes(tmp_path):
    cfg, paths, fig = _fixture(tmp_path)
    assert figure_status(fig, cfg, tmp_path, paths) == 'legacy_unverified'
    record_generated_figures(cfg, tmp_path, paths, {})
    assert figure_status(fig, cfg, tmp_path, paths) == 'valid'
    changed = deepcopy(cfg)
    changed['trajectory'] = {'stride': 4}
    assert figure_status(fig, changed, tmp_path, paths) == 'valid'
    changed['plots']['families'] = {'time_series': False}
    assert figure_status(fig, changed, tmp_path, paths) == 'disabled'
    fig.write_bytes(b'invalid PNG')
    assert figure_status(fig, cfg, tmp_path, paths) == 'invalid_figure'


def test_source_rewrite_keeps_figure_eligible_but_missing_source_does_not(tmp_path):
    cfg, paths, fig = _fixture(tmp_path)
    record_generated_figures(cfg, tmp_path, paths, {})
    source = paths['condition'][0] / 'summary/timeseries.parquet'
    pd.DataFrame({'time_ps': [0., 10.], 'rmsd_A': [7., 9.]}).to_parquet(source)
    assert figure_status(fig, cfg, tmp_path, paths) == 'valid'
    source.unlink()
    assert figure_status(fig, cfg, tmp_path, paths) == 'missing_source'


def test_disabled_interaction_family_never_reads_high_cardinality_data(tmp_path,monkeypatch):
    import mdscope.reporting.plotting as mod
    def forbidden(*a,**k):raise AssertionError('disabled plot read raw identities')
    monkeypatch.setattr(mod,'interaction_tables',forbidden)
    plot_interactions({'plots':{'families':{'interactions':False}}},tmp_path,{'x':[tmp_path/'rep1']})


def test_excluded_interaction_kind_is_filtered_before_read(tmp_path,monkeypatch):
    import mdscope.reporting.plotting as mod
    seen=[]
    def tables(*a,**k):
        seen.append(k['exclude_kinds']);return {}
    monkeypatch.setattr(mod,'interaction_tables',tables)
    plot_interactions({'plots':{'interaction_exclude_kinds':['lipid_hydrogen_bonds']}},tmp_path,{'x':[tmp_path/'rep1']})
    assert seen==[{'lipid_hydrogen_bonds'}]


def test_scientific_labels_preserve_raw_identifier_in_data_only():
    raw='sdc4_motif_acidic_patch_PROB_1_steric_occlusion_fraction'
    label,units=_metric_meta(raw)
    assert '_' not in label and 'chain B' in label and 'proxy' in label
    assert units=='fraction'
