"""Single-writer, saved-data PDF table publication with page identity.

A successful session replaces its reports and archives empty/unselected reports.
Rendering failures preserve committed PDFs. Unique staging files plus rollback
protect publication exceptions; process/power-loss recovery is not transactional.
"""
from __future__ import annotations
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
import hashlib
import json
import os
import tempfile

from matplotlib.backends.backend_pdf import PdfPages
from ..core.persistence import archive_output

_ACTIVE = ContextVar('mdscope_table_reports', default=None)
REPORT_NAME = 'statistics_tables.pdf'
MANIFEST_NAME = 'statistics_tables.manifest.json'


def pdf_digest(path):
    """Validate readable PDF pages and return a streamed integrity digest."""
    path = Path(path)
    with path.open('rb') as handle:
        if handle.read(5) != b'%PDF-':
            raise ValueError(f'Invalid PDF header: {path}')
        handle.seek(max(0, path.stat().st_size - 1024))
        if not handle.read().rstrip().endswith(b'%%EOF'):
            raise ValueError(f'Invalid PDF trailer: {path}')
        handle.seek(0)
        from pypdf import PdfReader
        from pypdf.errors import PyPdfError
        try:
            reader = PdfReader(handle, strict=True)
            if reader.is_encrypted or not len(reader.pages):
                raise ValueError('Expected unencrypted PDF pages')
            for page in reader.pages:
                if len(page.mediabox) != 4:
                    raise ValueError('Invalid page media box')
                content = page.get_contents()
                if content is not None:
                    content.get_data()
        except (PyPdfError, KeyError, TypeError, IndexError, RecursionError) as exc:
            raise ValueError(f'Invalid PDF structure: {path}: {exc}') from exc
        handle.seek(0)
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def report_manifest(path):
    """Read the publication-owned page manifest; missing legacy files stay unknown."""
    path = Path(path).with_name(MANIFEST_NAME)
    if not path.exists():
        return None
    value = json.loads(path.read_text())
    if (not isinstance(value, dict) or not isinstance(value.get('pages'), list) or
            any(not isinstance(page, dict) for page in value['pages'])):
        raise ValueError(f'Invalid PDF page manifest: {path}')
    return value


def table_page_enabled(cfg, identifier):
    """Bootstrap display switch does not disable comparison-table pages."""
    from ..config import figure_enabled
    return (figure_enabled(cfg, identifier) and
            (not identifier.endswith(':bootstrap_table') or
             cfg.get('plots', {}).get('bootstrap_tables', True)))


class TableReports:
    def __init__(self, root: Path, *, enabled=True, cfg=None):
        self.root = Path(root).resolve()
        self.enabled = enabled
        self.cfg = cfg
        self.open_reports = {}
        self.pending = {}
        self.pages = {}

    def add(self, fig, figure_path: Path, *, identifier=None):
        from .plot_ids import id_from_figure_path
        figure_path = Path(figure_path).resolve()
        figure_path.relative_to(self.root / 'project_summary/figures')
        report = figure_path.parent / REPORT_NAME
        if report not in self.open_reports:
            report.parent.mkdir(parents=True, exist_ok=True)
            fd, temporary = tempfile.mkstemp(prefix='.statistics_tables.', suffix='.pending.pdf', dir=report.parent)
            os.close(fd)
            self.pending[report] = Path(temporary)
            self.open_reports[report] = PdfPages(temporary, metadata={
                'Title': f'MDscope statistics tables - {report.parent.name}', 'Creator':'MDscope'})
            self.pages[report] = []
        self.open_reports[report].savefig(fig, bbox_inches='tight')
        self.pages[report].append({'page':len(self.pages[report])+1,
                                  'figure_id':identifier or id_from_figure_path(figure_path),
                                  'figure_path':str(figure_path.relative_to(self.root))})

    def _empty_status(self, report):
        from ..config import figure_enabled, plot_enabled
        from .figure_inventory import figure_family
        if not self.enabled or (self.cfg is not None and not plot_enabled(
                {**self.cfg, 'plots':{**self.cfg.get('plots', {}), 'include':None}}, figure_family(report))):
            return 'disabled'
        if self.cfg is not None:
            prior = report_manifest(report)
            pages = prior['pages'] if prior else []
            if pages and not any(page.get('figure_id') and table_page_enabled(self.cfg, page['figure_id']) for page in pages):
                return ('unselected' if not any(page.get('figure_id') and figure_enabled(self.cfg, page['figure_id']) for page in pages)
                        else 'disabled')
            if not pages and self.cfg.get('plots', {}).get('include') == []:
                return 'unselected'
        return 'empty'

    def close(self, success: bool):
        staged = []
        try:
            # Close every writer even if one close fails.
            failure = None
            for writer in self.open_reports.values():
                try:
                    writer.close()
                except Exception as exc:
                    failure = failure or exc
            if failure is not None:
                raise failure
            if not success:
                return
            existing = set((self.root/'project_summary/figures').rglob(REPORT_NAME))
            existing.update(path.with_name(REPORT_NAME) for path in
                            (self.root/'project_summary/figures').rglob(MANIFEST_NAME))
            manifests = {}
            for report in existing | set(self.pending):
                pages = self.pages.get(report, [])
                manifests[report] = {'status':'complete' if pages else self._empty_status(report),
                                     'format':'pdf', 'pages':pages, 'page_count':len(pages),
                                     'sha256':pdf_digest(self.pending[report]) if pages else None}
            # Stage every sidecar before replacing any committed report.
            sidecars = {}
            for report, value in manifests.items():
                fd, temporary = tempfile.mkstemp(prefix='.statistics_tables.', suffix='.json', dir=report.parent)
                staged.append(Path(temporary))
                os.close(fd)
                from ..config import atomic_json
                atomic_json(temporary, value)
                sidecars[report.with_name(MANIFEST_NAME)] = Path(temporary)
            replacements = {**self.pending, **sidecars}
            targets = set(replacements) | existing
            # Selected legacy table PNGs move only after PDF rendering succeeds.
            old_pngs = {self.root / page['figure_path'] for pages in self.pages.values()
                        for page in pages if (self.root / page['figure_path']).is_file()}
            targets |= old_pngs
            backups = {}
            published = []
            try:
                for target in sorted(targets | old_pngs):
                    if target.exists():
                        category = 'figure_tables' if target.suffix == '.png' else 'figure_tables/pdf_reports'
                        archive = self.root / '.archives' / category / target.parent.relative_to(self.root)
                        label = target.name
                        if target.suffix == '.png':
                            from datetime import datetime, timezone
                            label = f'{target.stem}_{datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")}.png'
                        backups[target] = archive_output(target, archive, label)
                for target, temporary in replacements.items():
                    os.replace(temporary, target)
                    published.append(target)
            except BaseException:
                for target in published:
                    target.unlink(missing_ok=True)
                for target, backup in backups.items():
                    os.replace(backup, target)
                raise
        finally:
            for path in [*self.pending.values(), *staged]:
                path.unlink(missing_ok=True)


@contextmanager
def table_report_session(root: Path, *, enabled: bool, cfg=None):
    # Only one writer per project is supported. Reject overlap rather than
    # merging pages from independently configured runs.
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock = root / '.table_reports.lock'
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f'A table-report writer is already active for {root}; inspect {lock} if an earlier process stopped') from exc
    os.close(fd)
    session = TableReports(root, enabled=enabled, cfg=cfg)
    token = _ACTIVE.set(session if enabled else None)
    try:
        try:
            yield
        except BaseException:
            session.close(success=False)
            raise
        else:
            session.close(success=True)
    finally:
        _ACTIVE.reset(token)
        lock.unlink(missing_ok=True)


def add_table_page(fig, path: Path, *, identifier=None):
    """Return True when the selected table was routed to the PDF writer."""
    session = _ACTIVE.get()
    if session is None:
        return False
    session.add(fig, path, identifier=identifier)
    return True
