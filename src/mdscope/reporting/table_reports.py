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
import tempfile

from matplotlib.backends.backend_pdf import PdfPages


_ACTIVE = ContextVar("mdscope_table_reports", default=None)


class TableReports:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.open_reports = {}
        self.pending = {}
        self.pages = {}

    def add(self, fig, figure_path: Path):
        directory = Path(figure_path).parent
        report = directory / "statistics_tables.pdf"
        if report not in self.open_reports:
            directory.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(prefix=".statistics_tables.", suffix=".pending.pdf",
                                         dir=directory, delete=False) as temporary:
                pending = Path(temporary.name)
            self.open_reports[report] = PdfPages(pending, metadata={
                "Title": f"MDscope statistics tables - {directory.name}",
                "Creator": "MDscope",
            })
            self.pending[report] = pending
            self.pages[report] = []
        self.open_reports[report].savefig(fig, bbox_inches="tight")
        from .plot_ids import id_from_figure_path
        self.pages[report].append(id_from_figure_path(figure_path) or str(figure_path.name))

    def close(self, success: bool):
        """Publish complete PDFs only; archive reports with no pages in this run.

        A reporting session is the sole writer for its output root. Unique
        temporary paths isolate independent sessions; they do not grant
        concurrent sessions permission to publish the same destination.
        """
        try:
            for writer in self.open_reports.values():
                writer.close()
            if success:
                # Closing all writers first avoids replacing a previous report
                # when a later PdfPages writer cannot finish.
                for pending in self.pending.values():
                    data = pending.read_bytes()
                    if not (data.startswith(b"%PDF-") and
                            b"%%EOF" in data[-1024:]):
                        raise ValueError(f"Invalid generated PDF: {pending}")
                for report, pending in self.pending.items():
                    os.replace(pending, report)
                    self._archive_old_pngs(report.parent)
                # Selection and available data may change between runs. A
                # report absent from a completed session must not remain active.
                active = set(self.pending)
                for report in sorted(self.root.rglob("statistics_tables.pdf")):
                    if ".archives" in report.parts or report in active:
                        continue
                    self._archive_obsolete_pdf(report)
                self._record_published_pages()
        finally:
            for pending in self.pending.values():
                pending.unlink(missing_ok=True)

    def _record_published_pages(self):
        """Keep report page identities associated with the current PDF set."""
        from ..config import atomic_json
        dest = self.root / "project_summary/figures/table_pdf_pages.json"
        if not self.pending and not dest.exists():
            return
        pages = {
            str(report.relative_to(self.root)): self.pages[report]
            for report in sorted(self.pending)
        }
        atomic_json(dest, pages)

    def _archive_obsolete_pdf(self, report: Path):
        relative = report.relative_to(self.root)
        dest = self.root / ".archives" / "figure_tables" / "pdf_reports" / relative.parent
        dest.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        # The random suffix prevents a same-microsecond archive collision.
        with tempfile.NamedTemporaryFile(prefix=report.stem + "_" + stamp + "_",
                                         suffix=report.suffix, dir=dest,
                                         delete=False) as temporary:
            archive = Path(temporary.name)
        try:
            os.replace(report, archive)
        finally:
            if report.exists():
                archive.unlink(missing_ok=True)

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
        # A successful plot-only rerun with PDF export disabled supersedes
        # previously published table PDFs, without touching numerical data.
        yield
        TableReports(root).close(success=True)
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
