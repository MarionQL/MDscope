"""Single-writer, saved-data PDF table publication with page identity.

A successful session replaces its reports and archives empty/unselected reports;
each archived PDF keeps its own manifest in one timestamped generation directory.
Rendering failures preserve committed PDFs. Unique staging files plus rollback
protect publication exceptions; process/power-loss recovery is not transactional.
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import UTC
from pathlib import Path

from matplotlib.backends.backend_pdf import PdfPages

from ..core.persistence import archive_output
from .plot_ids import figure_family

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
    try:
        value = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f'Invalid PDF page manifest: {path}: {exc}') from exc
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
        self.companions = {}
        self.metadata = {}

    def add(self, fig, figure_path: Path, *, identifier=None, source_ids=None):
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
        self.open_reports[report].savefig(fig, bbox_inches=None if source_ids else 'tight')
        self.pages[report].append({'page':len(self.pages[report])+1,
                                  'figure_id':identifier or id_from_figure_path(figure_path),
                                  'figure_path':figure_path.relative_to(self.root).as_posix(),
                                  **({'source_ids':source_ids} if source_ids else {})})

    def add_data(self, report, dataframe, metadata):
        report = Path(report)
        fd, temporary = tempfile.mkstemp(prefix='.statistics_tables.', suffix='.pending.csv', dir=report.parent)
        os.close(fd)
        self.pending[report.with_suffix('.csv')] = Path(temporary)
        dataframe.to_csv(temporary, index=False, na_rep='N/A')
        self.companions[report] = report.with_suffix('.csv')
        self.metadata[report] = metadata

    def _empty_status(self, report, superseded):
        """Classify a report without pages from the pages it last published."""
        from ..config import figure_enabled, plot_enabled
        if not self.enabled or (self.cfg is not None and not plot_enabled(
                {**self.cfg, 'plots':{**self.cfg.get('plots', {}), 'include':None}}, figure_family(report))):
            return 'disabled'
        if self.cfg is not None:
            ids = report_source_ids({'pages':superseded})
            if ids and not any(table_page_enabled(self.cfg, identifier) for identifier in ids):
                return ('unselected' if not any(figure_enabled(self.cfg, identifier) for identifier in ids)
                        else 'disabled')
            if not ids and self.cfg.get('plots', {}).get('include') == []:
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
            from datetime import datetime
            stamp = datetime.now(UTC).strftime('%Y%m%dT%H%M%S%fZ')
            figures = self.root / 'project_summary/figures'
            existing = set(figures.rglob(REPORT_NAME))
            existing.update(path.with_name(REPORT_NAME) for path in figures.rglob(MANIFEST_NAME))
            manifests = {}
            # Existing target -> archive directory, or None for a status-only
            # sidecar that is set aside only until publication succeeds.
            retired = {}
            for report in existing | set(self.pages):
                pages = self.pages.get(report, [])
                sidecar = report.with_name(MANIFEST_NAME)
                generation = None
                if report.exists():
                    # A published PDF and its sidecar are one archived generation.
                    generation = (self.root / '.archives/figure_tables/pdf_reports' /
                                  report.parent.relative_to(self.root) / stamp)
                    retired[report] = generation
                if sidecar.exists():
                    retired[sidecar] = generation
                csv = report.with_suffix('.csv')
                if csv.exists():
                    retired[csv] = generation
                if pages:
                    manifests[report] = {'status':'complete', 'format':'pdf', 'pages':pages,
                                         'page_count':len(pages), 'sha256':pdf_digest(self.pending[report]),
                                         **self.metadata.get(report, {})}
                    if report in self.companions:
                        csv = self.companions[report]
                        manifests[report]['csv'] = csv.relative_to(self.root).as_posix()
                        manifests[report]['csv_sha256'] = file_digest(self.pending[csv])
                    continue
                prior = report_manifest(report) or {}
                superseded = prior.get('pages') or prior.get('superseded_pages', [])
                manifests[report] = {
                    'status':self._empty_status(report, superseded), 'format':'pdf', 'pages':[],
                    'page_count':0, 'sha256':None, 'superseded_pages':superseded,
                    'archive':(generation.relative_to(self.root).as_posix() if generation
                               else prior.get('archive'))}
            # Stage every sidecar before replacing any committed report.
            sidecars = {}
            from ..core.persistence import atomic_json
            for report, value in manifests.items():
                fd, temporary = tempfile.mkstemp(prefix='.statistics_tables.', suffix='.json', dir=report.parent)
                staged.append(Path(temporary))
                os.close(fd)
                atomic_json(temporary, value)
                sidecars[report.with_name(MANIFEST_NAME)] = Path(temporary)
            replacements = {**self.pending, **sidecars}
            # Selected legacy table PNGs move only after PDF rendering succeeds.
            for pages in self.pages.values():
                for page in pages:
                    identifiers = page.get('source_ids') or [page.get('figure_id')]
                    for identifier in identifiers:
                        from ..registry import figure_descriptor
                        descriptor = figure_descriptor(self.cfg or {}, identifier)
                        for relative in (descriptor['destination'], descriptor['legacy_destination']):
                            png = self.root / relative
                            if png.is_file():
                                retired[png] = self.root / '.archives/figure_tables' / png.parent.relative_to(self.root)
            backups = {}
            published = []
            try:
                for target, archive in sorted(retired.items()):
                    if archive is None:
                        fd, aside = tempfile.mkstemp(prefix='.statistics_tables.', suffix='.previous.json',
                                                     dir=target.parent)
                        os.close(fd)
                        staged.append(Path(aside))
                        os.replace(target, aside)
                        backups[target] = Path(aside)
                    else:
                        label = f'{target.stem}_{stamp}.png' if target.suffix == '.png' else target.name
                        backups[target] = archive_output(target, archive, label)
                for target, temporary in replacements.items():
                    os.replace(temporary, target)
                    published.append(target)
            except BaseException:
                for target in published:
                    target.unlink(missing_ok=True)
                for target, backup in backups.items():
                    os.replace(backup, target)
                for archive in {archive for archive in retired.values() if archive is not None}:
                    if archive.is_dir() and not any(archive.iterdir()):
                        archive.rmdir()
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


def add_table_page(fig, path: Path, *, identifier=None, source_ids=None):
    """Return True when the selected table was routed to the PDF writer."""
    session = _ACTIVE.get()
    if session is None:
        return False
    session.add(fig, path, identifier=identifier, source_ids=source_ids)
    return True


def file_digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def report_source_ids(report):
    """Legacy single-source pages and consolidated pages share one reader."""
    return list(dict.fromkeys(identifier for page in report.get('pages', [])
                              for identifier in (page.get('source_ids') or [page.get('figure_id')])
                              if identifier))


def add_report_data(report, dataframe, metadata):
    session = _ACTIVE.get()
    if session is None:
        raise RuntimeError('Consolidated reports require an active publication session')
    session.add_data(report, dataframe, metadata)


def table_selection(cfg):
    """Settings that determine the exact population of a consolidated report."""
    plots = cfg.get('plots', {})
    return {key: (sorted(plots[key]) if key == 'include' and plots.get(key) is not None else plots.get(key, default))
            for key, default in (('include', None), ('bootstrap_tables', True), ('families', {}),
                                 ('time_series', True), ('profiles', True), ('interactions', True),
                                 ('interaction_exclude_kinds', []))}
