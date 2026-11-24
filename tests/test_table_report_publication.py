"""Regression tests for PDF table publication boundaries (issue #215)."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from mdscope.reporting.table_reports import add_table_page, table_report_session


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
    import json

    import pandas as pd

    from mdscope.reporting.figure_inventory import (
        figure_inventory,
        figure_snapshot,
        record_generated_figures,
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
    from mdscope.core import persistence
    folder = tmp_path/'project_summary/figures/time_series'
    with table_report_session(tmp_path, enabled=True):
        _page(folder/'old__comparison_table.png')
    before = {p:p.read_bytes() for p in (tmp_path/'project_summary/figures').rglob('*') if p.is_file()}
    def fail(*args, **kwargs):
        raise OSError('metadata disk full')
    monkeypatch.setattr(persistence,'atomic_json',fail)
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

    from mdscope.reporting.figure_inventory import figure_inventory
    from mdscope.reporting.plotting import _table_figure
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


def _manifest(folder):
    import json
    return json.loads((folder / "statistics_tables.manifest.json").read_text())


def test_superseded_generations_archive_pdf_with_its_own_manifest(tmp_path):
    import json

    from mdscope.reporting.table_reports import pdf_digest
    folder = tmp_path / "project_summary/figures/time_series"
    published = []
    for metric in ("first", "second", "third"):
        with table_report_session(tmp_path, enabled=True):
            _page(folder / f"{metric}__comparison_table.png")
        published.append(((folder / "statistics_tables.pdf").read_bytes(), _manifest(folder)))
    base = tmp_path / ".archives/figure_tables/pdf_reports/project_summary/figures/time_series"
    generations = sorted(base.iterdir())
    assert len(generations) == 2
    for generation, (data, manifest) in zip(generations, published[:2]):
        assert sorted(p.name for p in generation.iterdir()) == [
            "statistics_tables.manifest.json", "statistics_tables.pdf"]
        assert (generation / "statistics_tables.pdf").read_bytes() == data
        archived = json.loads((generation / "statistics_tables.manifest.json").read_text())
        assert archived == manifest
        assert archived["sha256"] == pdf_digest(generation / "statistics_tables.pdf")
    assert _manifest(folder)["pages"][0]["figure_path"] == (
        "project_summary/figures/time_series/third__comparison_table.png")


def test_repeated_unselected_session_keeps_status_identity_and_single_archive(tmp_path):
    from mdscope.reporting.figure_inventory import figure_inventory
    folder = tmp_path / "project_summary/figures/time_series"
    selected = {"plots": {"enabled": True, "table_reports_pdf": True}}
    with table_report_session(tmp_path, enabled=True, cfg=selected):
        _page(folder / "density__comparison_table.png")
    original = (folder / "statistics_tables.pdf").read_bytes()
    pages = _manifest(folder)["pages"]
    other = {"plots": {"enabled": True, "table_reports_pdf": True,
                       "include": ["profiles:lipid_order__POPC__sn1:profile"]}}
    archives = tmp_path / ".archives/figure_tables/pdf_reports"
    records = []
    for _ in range(2):
        with table_report_session(tmp_path, enabled=True, cfg=other):
            pass
        records.append(_manifest(folder))
        assert not (folder / "statistics_tables.pdf").exists()
        assert len([p for p in archives.rglob("*") if p.is_file()]) == 2
    assert records[0] == records[1]
    assert records[0]["status"] == "unselected"
    assert records[0]["superseded_pages"] == pages
    generation = tmp_path / records[0]["archive"]
    assert (generation / "statistics_tables.pdf").read_bytes() == original
    rows = figure_inventory(tmp_path, {}, other, include_rejected=True)
    row = rows.loc[rows["format"] == "pdf"].iloc[0]
    assert (row["status"], row["kind"], row["archive"]) == (
        "unselected", "table_report", records[0]["archive"])
    assert figure_inventory(tmp_path, {}, other).empty
    # Restoring the selection republishes without disturbing the archive.
    with table_report_session(tmp_path, enabled=True, cfg=selected):
        _page(folder / "density__comparison_table.png")
    assert _manifest(folder)["status"] == "complete"
    assert len([p for p in archives.rglob("*") if p.is_file()]) == 2


def test_empty_session_with_corrupt_prior_manifest_fails_and_preserves_outputs(tmp_path):
    folder = tmp_path / "project_summary/figures/time_series"
    with table_report_session(tmp_path, enabled=True):
        _page(folder / "density__comparison_table.png")
    (folder / "statistics_tables.manifest.json").write_text("{broken")
    before = {p: p.read_bytes() for p in folder.iterdir()}
    with pytest.raises(ValueError, match="Invalid PDF page manifest.*statistics_tables"):
        with table_report_session(tmp_path, enabled=True):
            pass
    assert {p: p.read_bytes() for p in folder.iterdir()} == before
    assert not (tmp_path / ".archives").exists()


def _production_fixture(tmp_path):
    from runpy import run_path

    import pandas as pd
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    cfg, repouts = run_path(str(scripts / "plot_preview.py"))["create_fixture"](tmp_path)
    cfg["plots"]["table_reports_pdf"] = True
    # Saved between-condition summaries; rendering only reads them.
    pd.DataFrame([{
        "metric": metric, "group": None, "system_1": "condition_A", "system_2": "condition_B",
        "mean_difference_system2_minus_system1": delta, "mean_difference_ci_low": delta - .1,
        "mean_difference_ci_high": delta + .1, "mean_probability_difference_gt_0": .9,
        "mean_p_value": .04, "mean_q_value": .08,
    } for metric, delta in (("rmsd_A", .2), ("radius_of_gyration_A", .4))]).to_parquet(
        tmp_path / "project_summary/statistics/time_series/system_comparisons.parquet")
    return cfg, repouts


def _figure_bytes(root):
    return {p: p.read_bytes() for p in (root / "project_summary/figures").rglob("*") if p.is_file()}


def test_production_plotting_inventories_each_pdf_once_with_page_identity(tmp_path, monkeypatch):
    import json
    from copy import deepcopy

    from mdscope.reporting import plotting
    from mdscope.reporting.figure_inventory import (
        figure_inventory,
        figure_snapshot,
        record_generated_figures,
    )
    cfg, repouts = _production_fixture(tmp_path)
    before = figure_snapshot(tmp_path, repouts)
    plotting.plot_all(cfg, tmp_path, repouts)
    record_generated_figures(cfg, tmp_path, repouts, before)
    series = "project_summary/figures/statistical_reports/rmsd/statistics_tables.pdf"
    profiles = "project_summary/figures/statistical_reports/membrane_order/statistics_tables.pdf"
    inventory = figure_inventory(tmp_path, repouts, cfg)
    reports = inventory.loc[inventory["format"] == "pdf"].set_index("figure_id")
    assert sorted(reports.index) == [profiles,
        "project_summary/figures/statistical_reports/radius_of_gyration/statistics_tables.pdf", series]
    assert set(reports["status"]) == {"valid"}
    assert set(reports["kind"]) == {"table_report"}
    assert json.loads(reports.loc[series, "page_ids"]) == [
        "time_series:rmsd_A:comparison_table"]
    assert json.loads(reports.loc[profiles, "page_ids"]) == [
        "profiles:lipid_order__POPC__sn1:comparison_table"]
    assert "project_summary/statistics/time_series/system_comparisons.parquet" in json.loads(
        reports.loc[series, "source_files"])
    provenance = json.loads((tmp_path / "project_summary/figures/figure_provenance.json").read_text())
    assert provenance[series]["format"] == "pdf"
    assert provenance[series]["source_metric"] is None
    assert provenance[series]["page_ids"] == json.loads(reports.loc[series, "page_ids"])
    png_only = figure_inventory(tmp_path, repouts, cfg, formats=("png",))
    assert not png_only.empty and set(png_only["format"]) == {"png"}
    assert not any(name.endswith("_table.png") for name in png_only["figure"])

    # One selected table page: no other table pages, and the unselected
    # profiles report is archived and reported rather than left as current.
    single = deepcopy(cfg)
    single["plots"]["include"] = ["time_series:rmsd_A:comparison_table"]
    before = figure_snapshot(tmp_path, repouts)
    plotting.plot_all(single, tmp_path, repouts)
    record_generated_figures(single, tmp_path, repouts, before)
    current = figure_inventory(tmp_path, repouts, single)
    assert current["figure_id"].tolist() == [series]
    assert json.loads(current.iloc[0]["page_ids"]) == ["time_series:rmsd_A:comparison_table"]
    assert current.iloc[0]["status"] == "valid"
    rejected = figure_inventory(tmp_path, repouts, single, include_rejected=True).set_index("figure_id")
    assert rejected.loc[profiles, "status"] == "unselected"
    assert (tmp_path / rejected.loc[profiles, "archive"] / "statistics_tables.pdf").is_file()

    # A failing renderer leaves every committed figure, PDF and sidecar unchanged.
    committed = _figure_bytes(tmp_path)

    def fail(*args, **kwargs):
        raise RuntimeError("profile renderer failed")
    monkeypatch.setattr(plotting, "plot_profiles", fail)
    with pytest.raises(RuntimeError, match="profile renderer failed"):
        plotting.plot_all(cfg, tmp_path, repouts)
    assert _figure_bytes(tmp_path) == committed
    assert not list(tmp_path.rglob("*.pending.pdf"))
    assert not (tmp_path / ".table_reports.lock").exists()
