"""Collect rendered statistics tables into one multipage PDF per figure folder.

The numerical Parquet exports remain untouched. Temporary PDFs are published
only after all rendering succeeds; older PNG tables are archived, not deleted.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path
import os

from matplotlib.backends.backend_pdf import PdfPages


_ACTIVE = ContextVar("mdscope_table_reports", default=None)


class TableReports:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.open_reports = {}
        self.pending = {}

    def add(self, fig, figure_path: Path):
        directory = Path(figure_path).parent
        report = directory / "statistics_tables.pdf"
        if report not in self.open_reports:
            directory.mkdir(parents=True, exist_ok=True)
            pending = directory / ".statistics_tables.pending.pdf"
            self.open_reports[report] = PdfPages(pending, metadata={
                "Title": f"MDscope statistics tables - {directory.name}",
                "Creator": "MDscope",
            })
            self.pending[report] = pending
        self.open_reports[report].savefig(fig, bbox_inches="tight")

    def close(self, success: bool):
        try:
            for writer in self.open_reports.values():
                writer.close()
            if success:
                for report, pending in self.pending.items():
                    os.replace(pending, report)
                    self._archive_old_pngs(report.parent)
        finally:
            if not success:
                for pending in self.pending.values():
                    pending.unlink(missing_ok=True)

    def _archive_old_pngs(self, directory: Path):
        old = sorted(set(directory.glob("*__bootstrap_table.png"))
                     | set(directory.glob("*__comparison_table.png")))
        if not old:
            return
        relative = directory.relative_to(self.root)
        dest = self.root / ".archives" / "figure_tables" / relative
        dest.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        for source in old:
            source.replace(dest / f"{source.stem}_{stamp}{source.suffix}")


@contextmanager
def table_report_session(root: Path, *, enabled: bool):
    if not enabled:
        yield
        return
    session = TableReports(root)
    token = _ACTIVE.set(session)
    try:
        yield
    except BaseException:
        session.close(success=False)
        raise
    else:
        session.close(success=True)
    finally:
        _ACTIVE.reset(token)


def add_table_page(fig, path: Path):
    """Return True when a table is routed to an active PDF session."""
    session = _ACTIVE.get()
    if session is None:
        return False
    session.add(fig, path)
    return True
