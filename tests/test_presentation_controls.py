"""Offline family selection and adversarial stale-figure integrity checks."""
import json
import zipfile
from copy import deepcopy
import numpy as np
import pandas as pd
import pytest
from PIL import Image
from mdscope.reporting.figure_inventory import (
    figure_snapshot,record_generated_figures,figure_inventory,figure_status)
from mdscope.reporting.plotting import plot_interactions, _metric_meta
from mdscope.reporting.report import generate_report


def _fixture(tmp_path):
    rep=tmp_path/'condition/rep1'; (rep/'summary').mkdir(parents=True)
    pd.DataFrame({'time_ps':[0.,10.,20.],'rmsd_A':[1.,2.,3.]}).to_parquet(rep/'summary/timeseries.parquet')
    fig=tmp_path/'project_summary/figures/time_series/rmsd_A__timeseries.png'
    fig.parent.mkdir(parents=True);Image.new('RGB',(50,50),'blue').save(fig)
    cfg={'project':{'name':'qa'},'systems':{'condition':{'replicates':[{'label':'rep1'}]}},
         'plots':{'enabled':True},'report':{'enabled':True},'statistics':{'enabled':True}}
    return cfg,{'condition':[rep]},fig


def test_unverified_stale_disabled_and_changed_figures_are_distinct(tmp_path):
    cfg,paths,fig=_fixture(tmp_path)
    assert figure_status(fig,cfg,tmp_path,paths)=='unverified'
    record_generated_figures(cfg,tmp_path,paths,{})
    assert figure_status(fig,cfg,tmp_path,paths)=='valid'
    changed=deepcopy(cfg);changed['plots']['families']={'time_series':False}
    assert figure_status(fig,changed,tmp_path,paths)=='disabled'
    changed=deepcopy(cfg);changed['trajectory']={'stride':4}
    assert figure_status(fig,changed,tmp_path,paths)=='stale_config'
    for key, value in [('time', {'dt_ps': 20.}), ('selections', {'protein': 'segid B'})]:
        changed = deepcopy(cfg); changed[key] = value
        assert figure_status(fig,changed,tmp_path,paths)=='stale_config'
    Image.new('RGB',(50,50),'red').save(fig)
    assert figure_status(fig,cfg,tmp_path,paths)=='changed_figure'
    before=figure_snapshot(tmp_path,paths)
    record_generated_figures(cfg,tmp_path,paths,before)
    assert figure_status(fig,cfg,tmp_path,paths)=='changed_figure' # unchanged PNG is not blessed


def test_source_change_and_missing_source_invalidate_figure(tmp_path):
    cfg,paths,fig=_fixture(tmp_path);record_generated_figures(cfg,tmp_path,paths,{})
    source=paths['condition'][0]/'summary/timeseries.parquet'
    pd.DataFrame({'time_ps':[0.,10.],'rmsd_A':[7.,9.]}).to_parquet(source)
    assert figure_status(fig,cfg,tmp_path,paths)=='stale_source'
    source.unlink()
    assert figure_status(fig,cfg,tmp_path,paths)=='missing_source'


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


def test_pdf_bundle_excludes_unverified_png_and_preserves_auditable_rejection(tmp_path):
    cfg,paths,fig=_fixture(tmp_path)
    result=generate_report(cfg,tmp_path,paths)
    with zipfile.ZipFile(result['bundle']) as z:
        assert str(fig.relative_to(tmp_path)) not in z.namelist()
        manifest=pd.read_csv(z.open('project_summary/report/figure_manifest.csv'))
        assert manifest.status.tolist()==['unverified']
    record_generated_figures(cfg,tmp_path,paths,{})
    result=generate_report(cfg,tmp_path,paths)
    with zipfile.ZipFile(result['bundle']) as z:
        assert str(fig.relative_to(tmp_path)) in z.namelist()
    assert json.loads((tmp_path/'project_summary/report/report_provenance.json').read_text())['sampling'][0]['frames']==3


def test_scientific_labels_preserve_raw_identifier_in_data_only():
    raw='sdc4_motif_acidic_patch_PROB_1_steric_occlusion_fraction'
    label,units=_metric_meta(raw)
    assert '_' not in label and 'chain B' in label and 'proxy' in label
    assert units=='fraction'


def test_pdf_section_controls_and_q_not_interval_sign(tmp_path):
    from pypdf import PdfReader
    from mdscope.reporting.report_inventory import inference_status
    cfg,paths,fig=_fixture(tmp_path)
    cfg['report']['sections']={'overview':False,'solvation':False,'atlas':False}
    result=generate_report(cfg,tmp_path,paths)
    text=' '.join(page.extract_text() for page in PdfReader(result['pdf']).pages)
    assert 'Systems and analyzed sampling' not in text
    assert 'Scientific figure atlas and coverage audit' not in text
    assert 'Condition-level scalar results' in text
    assert inference_status(.99,.1)=='Exploratory'
    assert inference_status(.04,.01)=='Adjusted support'
    assert inference_status(np.nan,np.nan)=='Descriptive / unadjusted'


def test_disabled_inline_solvation_figures_are_not_constructed(tmp_path,monkeypatch):
    import mdscope.reporting.report as report
    cfg,paths,fig=_fixture(tmp_path)
    cfg['plots']['families']={'solvation':False}
    folder=paths['condition'][0]/'extensions/solvation';folder.mkdir(parents=True)
    pd.DataFrame({'time_ps':[0.], 'solvent':['water'], 'target_class':['protein'],
                  'site_class':['acidic'], 'mean_n_solvent':[1.]}).to_parquet(folder/'species_timeseries.parquet')
    def forbidden(*a,**k):raise AssertionError('disabled inline plot was constructed')
    for name in ('solvation_rdf_figure','solvation_timeseries_figure','solvation_survival_figure'):
        monkeypatch.setattr(report,name,forbidden)
    generate_report(cfg,tmp_path,paths)


def test_bundle_coverage_uses_the_active_family_controls(tmp_path):
    from mdscope.reporting.figure_inventory import report_coverage
    cfg, paths, _ = _fixture(tmp_path)
    cfg["plots"]["families"] = {"time_series": False}
    expected = report_coverage(tmp_path, paths, cfg)
    result = generate_report(cfg, tmp_path, paths)
    with zipfile.ZipFile(result["bundle"]) as archive:
        actual = pd.read_csv(archive.open("project_summary/report/figure_coverage.csv"))
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    assert actual.loc[actual.analysis.eq("Time series"), "status"].item() == "disabled"
