"""Regression tests for PDF table publication boundaries (issue #215)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from mdscope.reporting.table_reports import table_report_session, add_table_page


def _page(path: Path):
    fig, ax = plt.subplots()
    try:
        ax.text(0.5, 0.5, path.stem)
        assert add_table_page(fig, path)
    finally:
        plt.close(fig)


def test_empty_success_archives_old_pdf(tmp_path):
    selected = tmp_path / "project_summary" / "figures" / "time_series"
    selected.mkdir(parents=True)
    report = selected / "statistics_tables.pdf"
    report.write_bytes(b"old report")
    with table_report_session(tmp_path, enabled=True):
        pass
    assert not report.exists()
    archives = list((tmp_path / ".archives" / "figure_tables" / "pdf_reports").rglob("*.pdf"))
    assert len(archives) == 1
    assert archives[0].read_bytes() == b"old report"


def test_new_pdf_is_published_and_old_png_archived(tmp_path):
    directory = tmp_path / "project_summary" / "figures" / "profiles"
    directory.mkdir(parents=True)
    png = directory / "metric__bootstrap_table.png"
    png.write_bytes(b"legacy png")
    with table_report_session(tmp_path, enabled=True):
        _page(directory / "metric__bootstrap_table.png")
    pdf = directory / "statistics_tables.pdf"
    assert pdf.read_bytes().startswith(b"%PDF-")
    assert b"%%EOF" in pdf.read_bytes()[-1024:]
    assert not png.exists()
    assert list((tmp_path / ".archives" / "figure_tables").rglob("metric__bootstrap_table_*.png"))


def test_failed_session_preserves_prior_report_and_removes_temp(tmp_path):
    directory = tmp_path / "project_summary" / "figures" / "profiles"
    directory.mkdir(parents=True)
    report = directory / "statistics_tables.pdf"
    report.write_bytes(b"previous complete PDF")
    with pytest.raises(RuntimeError, match="renderer failed"):
        with table_report_session(tmp_path, enabled=True):
            _page(directory / "metric__comparison_table.png")
            raise RuntimeError("renderer failed")
    assert report.read_bytes() == b"previous complete PDF"
    assert not list(directory.glob("*.pending.pdf"))


def test_unique_pending_paths_for_distinct_sessions(tmp_path):
    directory = tmp_path / "project_summary" / "figures" / "profiles"
    directory.mkdir(parents=True)
    paths = []
    for _ in range(2):
        with pytest.raises(RuntimeError):
            with table_report_session(tmp_path, enabled=True):
                _page(directory / "metric__comparison_table.png")
                paths.extend(directory.glob("*.pending.pdf"))
                raise RuntimeError("rollback")
    assert len(paths) == 2
    assert paths[0] != paths[1]
    assert not list(directory.glob("*.pending.pdf"))


def test_archives_only_reports_not_regenerated(tmp_path):
    old = tmp_path / "project_summary" / "figures" / "profiles" / "statistics_tables.pdf"
    old.parent.mkdir(parents=True)
    old.write_bytes(b"obsolete")
    active = tmp_path / "project_summary" / "figures" / "time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(active / "metric__comparison_table.png")
    assert not old.exists()
    assert (active / "statistics_tables.pdf").exists()


def test_pdf_inventory_and_snapshot_include_published_report(tmp_path):
    from mdscope.reporting.figure_inventory import figure_inventory, figure_snapshot
    directory = tmp_path / "project_summary" / "figures" / "time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(directory / "density__comparison_table.png")
    snapshot = figure_snapshot(tmp_path, {})
    relative = "project_summary/figures/time_series/statistics_tables.pdf"
    assert relative in snapshot
    inventory = figure_inventory(tmp_path, {})
    row = inventory.loc[inventory["figure_id"] == relative].iloc[0]
    assert row["format"] == "pdf"
    assert row["kind"] == "table_report"
    assert "time_series:density:comparison_table" in row["page_ids"]


def test_disabled_pdf_mode_archives_old_reports(tmp_path):
    directory = tmp_path / "project_summary" / "figures" / "profiles"
    directory.mkdir(parents=True)
    report = directory / "statistics_tables.pdf"
    report.write_bytes(b"obsolete")
    with table_report_session(tmp_path, enabled=False):
        pass
    assert not report.exists()


def test_pdf_inventory_rejects_corruption(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    directory = tmp_path / "project_summary" / "figures" / "time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(directory / "density__comparison_table.png")
    report = directory / "statistics_tables.pdf"
    report.write_bytes(b"%PDF-1.4\ntruncated")
    cfg = {"plots": {"table_reports_pdf": True}}
    assert figure_status(report, cfg, tmp_path, {}) == "invalid_figure"


def test_pdf_inventory_honors_explicit_empty_include(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    directory = tmp_path / "project_summary" / "figures" / "time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(directory / "density__comparison_table.png")
    report = directory / "statistics_tables.pdf"
    cfg = {"plots": {"table_reports_pdf": True, "include": []}}
    assert figure_status(report, cfg, tmp_path, {}) == "disabled"


def test_second_report_replace_failure_rolls_back_first(tmp_path, monkeypatch):
    import mdscope.reporting.table_reports as publication
    one = tmp_path / "project_summary/figures/time_series"
    two = tmp_path / "project_summary/figures/profiles"
    for directory in (one, two):
        directory.mkdir(parents=True)
        (directory / "statistics_tables.pdf").write_bytes(
            f"previous {directory.name}".encode())
    real_replace = publication.os.replace
    failed = False

    def fail_second(source, destination):
        nonlocal failed
        if (not failed and Path(destination) == two / "statistics_tables.pdf"
                and str(source).endswith(".pending.pdf")):
            failed = True
            raise OSError("simulated second report publication failure")
        return real_replace(source, destination)

    monkeypatch.setattr(publication.os, "replace", fail_second)
    with pytest.raises(OSError, match="second report"):
        with table_report_session(tmp_path, enabled=True):
            _page(one / "density__comparison_table.png")
            _page(two / "height__bootstrap_table.png")
    assert (one / "statistics_tables.pdf").read_bytes() == b"previous time_series"
    assert (two / "statistics_tables.pdf").read_bytes() == b"previous profiles"
    assert not list(tmp_path.rglob("*.pending.pdf"))
    assert not list(tmp_path.rglob("*.backup.pdf"))


def test_production_provenance_records_pdf_page_sources_and_sampling(tmp_path):
    import pandas as pd
    import json
    from mdscope.reporting.figure_inventory import (
        figure_snapshot, figure_inventory, record_generated_figures,
    )
    rep = tmp_path / "control/rep1"
    (rep / "summary").mkdir(parents=True)
    pd.DataFrame({"time_ps": [0.0, 10.0]}).to_parquet(
        rep / "summary/timeseries.parquet")
    stats = tmp_path / "project_summary/statistics/time_series"
    stats.mkdir(parents=True)
    pd.DataFrame({"metric": ["density"]}).to_parquet(
        stats / "replicate_statistics.parquet")
    reps = {"control": [rep]}
    before = figure_snapshot(tmp_path, reps)
    folder = tmp_path / "project_summary/figures/time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(folder / "density__comparison_table.png")
    cfg = {"plots": {"enabled": True, "table_reports_pdf": True}}
    record_generated_figures(cfg, tmp_path, reps, before)
    inventory = figure_inventory(tmp_path, reps, cfg=cfg)
    report = inventory.loc[inventory["format"] == "pdf"].iloc[0]
    assert report["status"] == "valid"
    assert "replicate_statistics.parquet" in report["source_files"]
    assert "rep1" in report["sampling"]
    assert "time_series:density:comparison_table" in report["page_ids"]
    provenance = json.loads((
        tmp_path / "project_summary/figures/figure_provenance.json").read_text())
    assert report["figure_id"] in provenance


def test_metadata_staging_failure_preserves_committed_pdf_and_page_identity(tmp_path, monkeypatch):
    from mdscope import config
    folder = tmp_path/'project_summary/figures/time_series'
    with table_report_session(tmp_path, enabled=True):
        _page(folder/'old__comparison_table.png')
    before = {p:p.read_bytes() for p in (tmp_path/'project_summary/figures').rglob('*') if p.is_file()}
    def fail(*args, **kwargs):
        raise OSError('metadata disk full')
    monkeypatch.setattr(config,'atomic_json',fail)
    with pytest.raises(OSError,match='metadata disk full'):
        with table_report_session(tmp_path, enabled=True):
            _page(folder/'new__comparison_table.png')
    assert {p:p.read_bytes() for p in (tmp_path/'project_summary/figures').rglob('*') if p.is_file()} == before


def test_pdf_with_signature_but_no_page_structure_is_rejected(tmp_path):
    from mdscope.reporting.figure_inventory import figure_status
    report=tmp_path/'project_summary/figures/time_series/statistics_tables.pdf'
    report.parent.mkdir(parents=True)
    report.write_bytes(b'%PDF-1.4\nTHIS IS NOT A PDF\n%%EOF\n')
    assert figure_status(report,{'plots':{'table_reports_pdf':True}},tmp_path,{}) == 'invalid_figure'


def test_page_manifest_replace_failure_restores_publication_and_selected_png(tmp_path, monkeypatch):
    import mdscope.reporting.table_reports as publication
    folder=tmp_path/'project_summary/figures/profiles'
    with table_report_session(tmp_path,enabled=True):
        _page(folder/'old__bootstrap_table.png')
    png=folder/'new__bootstrap_table.png'; png.write_bytes(b'old numerical table image')
    before={p:p.read_bytes() for p in folder.iterdir()}
    replace=publication.os.replace
    def fail(source,target):
        if Path(target).name == publication.MANIFEST_NAME:
            if Path(source).name.startswith('.statistics_tables.'):
                raise OSError('manifest replacement failure')
        return replace(source,target)
    monkeypatch.setattr(publication.os,'replace',fail)
    with pytest.raises(OSError,match='manifest replacement'):
        with table_report_session(tmp_path,enabled=True):
            _page(png)
    assert {p:p.read_bytes() for p in folder.iterdir()} == before


def test_overlapping_project_writers_fail_without_corrupting_first_session(tmp_path):
    with table_report_session(tmp_path,enabled=True):
        with pytest.raises(RuntimeError,match='already active'):
            with table_report_session(tmp_path,enabled=True):
                pass
        _page(tmp_path/'project_summary/figures/profiles/chosen__bootstrap_table.png')
    assert not (tmp_path/'.table_reports.lock').exists()


def test_disabling_bootstrap_pages_keeps_comparison_pages_and_png_only_inventory(tmp_path):
    import pandas as pd
    from mdscope.reporting.plotting import _table_figure
    from mdscope.reporting.figure_inventory import figure_inventory
    cfg={'plots':{'enabled':True,'bootstrap_tables':False,'table_reports_pdf':True}}
    folder=tmp_path/'project_summary/figures/time_series'
    with table_report_session(tmp_path,enabled=True,cfg=cfg):
        for variant in ('bootstrap_table','comparison_table'):
            _table_figure(variant,pd.DataFrame({'mean':[1.]}),folder/f'rmsd_A__{variant}.png',cfg)
    result=figure_inventory(tmp_path,{},cfg)
    assert result.iloc[0]['page_ids'] == '["time_series:rmsd_A:comparison_table"]'
    assert result.iloc[0]['page_count'] == 1
    assert figure_inventory(tmp_path,{},cfg,formats=('png',)).empty


def test_bad_pdf_metadata_does_not_break_png_inventory(tmp_path):
    from PIL import Image
    from mdscope.reporting.figure_inventory import figure_inventory
    folder=tmp_path/'project_summary/figures/time_series'
    with table_report_session(tmp_path,enabled=True):
        _page(folder/'rmsd_A__comparison_table.png')
    (folder/'statistics_tables.manifest.json').write_text('{broken')
    Image.new('RGB',(10,10),'red').save(folder/'rmsd_A__timeseries.png')
    cfg={'plots':{'enabled':True,'table_reports_pdf':True}}
    result=figure_inventory(tmp_path,{},cfg)
    assert result['format'].tolist() == ['png']
    rejected=figure_inventory(tmp_path,{},cfg,include_rejected=True)
    assert 'invalid_manifest' in set(rejected.status)


def test_pdf_table_formats_numeric_cells_in_mixed_columns_without_altering_data(tmp_path):
    import pandas as pd
    from pypdf import PdfReader
    from mdscope.reporting.plotting import _table_figure
    frame=pd.DataFrame({'level':['condition_A/rep1','condition_A (system)'],
                        'tau(fr)':[11.9818735105,'replicate-specific']})
    original=frame.copy(deep=True)
    folder=tmp_path/'project_summary/figures/time_series'
    with table_report_session(tmp_path,enabled=True):
        _table_figure('Synthetic mixed table',frame,folder/'rmsd_A__bootstrap_table.png',
                      {'plots':{'enabled':True,'table_reports_pdf':True}})
    text=PdfReader(folder/'statistics_tables.pdf').pages[0].extract_text()
    assert '11.9818735105' not in text and '12' in text
    assert 'replicate-specific' in text
    pd.testing.assert_frame_equal(frame,original)
