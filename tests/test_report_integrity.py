"""Offline export failure, stale numerical inputs, and auditable bundle contracts."""
from copy import deepcopy
import hashlib
import json
import subprocess
import sys
import zipfile
import pandas as pd
import pytest
from pypdf import PdfReader
from mdscope.reporting.report import generate_report
from mdscope.reporting.figure_inventory import statistics_status
from mdscope.statistics import analyze_project
from test_presentation_controls import _fixture


def _statistical_fixture(root):
    cfg, paths, fig = _fixture(root)
    cfg['statistics'].update(bootstrap_draws=20, families={
        'profiles': False, 'interactions': False, 'solvation_rdfs': False})
    analyze_project(cfg, root, paths)
    return cfg, paths


@pytest.mark.parametrize('damage,expected', [('running', 'incomplete'), ('malformed', 'invalid_manifest')])
def test_partial_statistics_never_enter_pdf_or_current_bundle(tmp_path, damage, expected):
    cfg, paths = _statistical_fixture(tmp_path)
    manifest = tmp_path / 'project_summary/statistics/manifest.json'
    if damage == 'running':
        data = json.loads(manifest.read_text())
        data['status'] = 'running'
        manifest.write_text(json.dumps(data))
    else:
        manifest.write_text('broken json')
    assert statistics_status(cfg, tmp_path, 'time_series') == expected
    result = generate_report(cfg, tmp_path, paths)
    text = ' '.join(page.extract_text() for page in PdfReader(result['pdf']).pages)
    assert 'Statistical inputs excluded: ' + expected in text
    with zipfile.ZipFile(result['bundle']) as archive:
        names = archive.namelist()
        assert len(names) == len(set(names))
        assert not any('/time_series/system_statistics.' in name for name in names)
        audit = pd.read_csv(archive.open('project_summary/report/table_manifest.csv'))
        assert expected in audit.status.tolist()


def test_old_statistics_use_recorded_methods_after_config_changes(tmp_path):
    cfg, paths = _statistical_fixture(tmp_path)
    manifest = tmp_path / 'project_summary/statistics/manifest.json'
    data = json.loads(manifest.read_text())
    data.update(statistics_hash='old', analysis_config_hash='old', input_fingerprint='old')
    manifest.write_text(json.dumps(data))
    cfg['statistics']['ci'] = .90
    cfg['trajectory'] = {'stride': 7}
    assert statistics_status(cfg, tmp_path, 'time_series') == 'valid'
    result = generate_report(cfg, tmp_path, paths)
    text = ' '.join(page.extract_text() for page in PdfReader(result['pdf']).pages)
    assert '95% CI' in text and '90% CI' not in text


@pytest.mark.parametrize('stage', ['pdf', 'bundle'])
def test_failed_export_preserves_previous_pdf_bundle_and_manifests(tmp_path,monkeypatch,stage):
    import mdscope.reporting.report as report
    cfg,paths,_ = _fixture(tmp_path)
    result = generate_report(cfg,tmp_path,paths)
    out = tmp_path/'project_summary/report'
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    def fail(*args,**kwargs): raise RuntimeError('injected export failure')
    monkeypatch.setattr(report, '_build_pdf' if stage=='pdf' else '_make_bundle', fail)
    with pytest.raises(RuntimeError,match='injected'): generate_report(cfg,tmp_path,paths)
    after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    assert before==after
    assert not any(p.is_dir() for p in out.iterdir())


def test_current_bundle_preserves_full_parquet_and_current_sidecars(tmp_path):
    cfg,paths = _statistical_fixture(tmp_path)
    first = generate_report(cfg,tmp_path,paths)
    cfg['report']['max_metrics_per_section']=10
    result=generate_report(cfg,tmp_path,paths)
    with zipfile.ZipFile(result['bundle']) as archive:
        original=tmp_path/'project_summary/statistics/time_series/system_statistics.parquet'
        assert archive.read(str(original.relative_to(tmp_path)))==original.read_bytes()
        actual=json.loads(archive.read('project_summary/report/report_provenance.json'))
        assert actual['report_options']==cfg['report']
        assert len(archive.namelist())==len(set(archive.namelist()))
    text=' '.join(page.extract_text() for page in PdfReader(result['pdf']).pages)
    assert 'legacy/unverified' not in text


def test_report_components_import_without_pdf_plotting_or_trajectory_dependencies():
    code='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname.split('.')[0] in {'reportlab','matplotlib','MDAnalysis'}:
            raise AssertionError('unexpected dependency: '+fullname)
sys.meta_path.insert(0,Block())
import mdscope.reporting.report_inventory
import mdscope.reporting.bundle
import mdscope.reporting.report
'''
    subprocess.run([sys.executable,'-c',code],check=True,capture_output=True,text=True)


def test_duplicate_pngs_are_audited_once_without_duplicate_pdf_panels(tmp_path):
    from mdscope.reporting.figure_inventory import record_generated_figures, figure_inventory
    cfg,paths,fig=_fixture(tmp_path)
    duplicate=fig.with_name('duplicate__timeseries.png')
    duplicate.write_bytes(fig.read_bytes())
    record_generated_figures(cfg,tmp_path,paths,{})
    inventory=figure_inventory(tmp_path,paths,cfg,include_rejected=True)
    assert sorted(inventory.status)==['duplicate_figure','valid']
    result=generate_report(cfg,tmp_path,paths)
    content=json.loads((tmp_path/'project_summary/report/report_content.json').read_text())
    assert len(content['figures'])==1
    with zipfile.ZipFile(result['bundle']) as archive:
        assert sum(name.endswith('.png') for name in archive.namelist())==1
        source=tmp_path/'project_summary/report/figure_provenance.json'
        assert archive.read('project_summary/report/figure_provenance.json')==source.read_bytes()


def test_human_labels_ci_level_and_raw_identifier_companions(tmp_path):
    cfg,paths,_=_fixture(tmp_path)
    metric='sdc4_motif_acidic_patch_PROB_1_protein_neighbor_residues'
    source=paths['condition'][0]/'summary/timeseries.parquet'
    pd.DataFrame({'time_ps':[0.,10.,20.],metric:[1.,2.,3.]}).to_parquet(source)
    cfg['statistics'].update(bootstrap_draws=20,ci=.9,families={
        'profiles':False,'interactions':False,'solvation_rdfs':False})
    analyze_project(cfg,tmp_path,paths)
    result=generate_report(cfg,tmp_path,paths)
    text=' '.join(page.extract_text() for page in PdfReader(result['pdf']).pages)
    assert '90% CI' in text and '95% CI' not in text
    assert metric not in text
    assert 'neighbor residues' in text and 'Units' in text
    with zipfile.ZipFile(result['bundle']) as archive:
        data=pd.read_csv(archive.open('project_summary/statistics/time_series/system_statistics.csv'))
        assert data.metric.tolist()==[metric]


def test_inputs_changed_during_statistics_leave_an_incomplete_manifest(tmp_path, monkeypatch):
    import mdscope.statistics as statistics
    cfg,paths,_=_fixture(tmp_path)
    cfg['statistics']['bootstrap_draws']=20
    source=paths['condition'][0]/'summary/timeseries.parquet'
    def change_input(*args):
        data=pd.read_parquet(source);data['rmsd_A']+=20.;data.to_parquet(source)
    monkeypatch.setattr(statistics,'analyze_profiles',change_input)
    with pytest.raises(RuntimeError,match='inputs changed'):
        analyze_project(cfg,tmp_path,paths)
    assert statistics_status(cfg,tmp_path,'time_series')=='incomplete'


def test_disabling_bundle_archives_the_previous_export(tmp_path):
    cfg,paths,_=_fixture(tmp_path)
    result=generate_report(cfg,tmp_path,paths)
    old_bundle=__import__('pathlib').Path(result['bundle'])
    content=old_bundle.read_bytes()
    cfg['report']['bundle_enabled']=False
    result=generate_report(cfg,tmp_path,paths)
    assert result['bundle'] is None and not old_bundle.exists()
    archived=list((tmp_path/'.archives/report/bundles').glob('*.zip'))
    assert len(archived)==1 and archived[0].read_bytes()==content
