"""Per-folder PDF table reports and compact lipid display-label regressions."""

import pandas as pd
import pytest

from mdscope.reporting.lipid_contact_bootstrap import lipid_display_label
from mdscope.reporting.plotting import _table_figure
from mdscope.reporting.table_reports import table_report_session


def _cfg(*, include=None):
    return {"plots": {"enabled": True, "table_reports_pdf": True,
                      "include": include, "dpi": 72}}


def test_pdf_replaces_table_pngs_and_groups_multiple_tables(tmp_path):
    figures = tmp_path / "project_summary/figures"
    folder = figures / "time_series"
    folder.mkdir(parents=True)
    stale = folder / "rmsd_A__bootstrap_table.png"
    stale.write_bytes(b"old PNG artifact")

    with table_report_session(tmp_path, enabled=True):
        _table_figure("RMSD bootstrap", pd.DataFrame({"mean": [1.0, 2.0]}),
                      stale, _cfg())
        _table_figure(
            "RMSD comparisons",
            pd.DataFrame({"delta": [0.1], "CI": ["-0.2 to 0.3"]}),
            folder / "rmsd_A__comparison_table.png", _cfg(),
        )

    report = folder / "statistics_tables.pdf"
    assert report.read_bytes().startswith(b"%PDF-")
    assert report.stat().st_size > 1000
    assert not list(folder.glob("*__bootstrap_table.png"))
    assert not list(folder.glob("*__comparison_table.png"))
    archived = list((tmp_path / ".archives/figure_tables"
                     / "project_summary/figures/time_series").glob("rmsd_A__bootstrap_table_*.png"))
    assert len(archived) == 1
    assert archived[0].read_bytes() == b"old PNG artifact"


def test_pdf_reports_split_by_folder_and_skip_empty_data(tmp_path):
    with table_report_session(tmp_path, enabled=True):
        for folder in ("profiles", "interactions/hydrogen_bonds"):
            out = tmp_path / "project_summary/figures" / folder
            _table_figure("Complete", pd.DataFrame({"value": [0.]}),
                          out / "metric__bootstrap_table.png", _cfg())
            _table_figure("Empty", pd.DataFrame({"value": []}),
                          out / "empty__bootstrap_table.png", _cfg())
    reports = sorted((tmp_path / "project_summary/figures").rglob("statistics_tables.pdf"))
    assert len(reports) == 2
    assert all(p.stat().st_size > 1000 for p in reports)


def test_failed_report_does_not_replace_existing_pdf_or_png(tmp_path):
    folder = tmp_path / "project_summary/figures/profiles"
    folder.mkdir(parents=True)
    report = folder / "statistics_tables.pdf"
    report.write_bytes(b"existing report")
    legacy = folder / "metric__bootstrap_table.png"
    legacy.write_bytes(b"existing PNG")
    with pytest.raises(RuntimeError, match="abort"):
        with table_report_session(tmp_path, enabled=True):
            _table_figure("New", pd.DataFrame({"value": [2.]}),
                          legacy, _cfg())
            raise RuntimeError("abort")
    assert report.read_bytes() == b"existing report"
    assert legacy.read_bytes() == b"existing PNG"
    assert not (folder / ".statistics_tables.pending.pdf").exists()


def test_figure_selection_keeps_unselected_pages_out(tmp_path):
    folder = tmp_path / "project_summary/figures/time_series"
    with table_report_session(tmp_path, enabled=True):
        _table_figure("Do not include", pd.DataFrame({"value": [1.]}),
                      folder / "unused__bootstrap_table.png",
                      _cfg(include=["time_series:chosen:bootstrap_table"]))
    assert not (folder / "statistics_tables.pdf").exists()


def test_legacy_png_mode_preserved(tmp_path):
    folder = tmp_path / "project_summary/figures/time_series"
    with table_report_session(tmp_path, enabled=False):
        _table_figure("Legacy", pd.DataFrame({"value": [1.]}),
                      folder / "rmsd_A__bootstrap_table.png", _cfg())
    assert (tmp_path / "project_summary/figures/rmsd/secondary/time_series@rmsd_A@bootstrap_table.png").is_file()
    assert not (folder / "statistics_tables.pdf").exists()


def test_lipid_labels_use_residue_codes_or_short_class_names():
    assert lipid_display_label("POPC", level="species") == "POPC"
    assert lipid_display_label("GM1a", level="species") == "GM1a"
    assert lipid_display_label("phosphatidylcholine", level="class") == "PC"
    assert lipid_display_label("phosphatidylethanolamine", level="class") == "POPE"
    assert lipid_display_label("sphingomyelin", level="class") == "PSM"
    assert lipid_display_label("unknown_class", level="class") == "unknown_class"
