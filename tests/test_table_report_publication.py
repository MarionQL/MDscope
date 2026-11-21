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
