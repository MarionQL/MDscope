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
        ax.text(0.5, 0.5, "saved statistics")
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
