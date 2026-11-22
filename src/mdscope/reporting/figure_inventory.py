"""Inventory and provenance for manually generated scientific figures."""
from __future__ import annotations

from datetime import UTC
from pathlib import Path

import pandas as pd

PRIMARY_SUFFIXES = (
    "__timeseries.png", "__profile.png", "__system_difference.png",
    "occupancy_heatmap.png", "lifetime_distribution.png",
    "distance_heatmap.png", "lipid_enrichment_heatmap.png",
    "state_populations.png", "state_sampling_diagnostics.png",
    "lipid_density_maps.png", "lipid_enrichment_maps.png",
    "lipid_region_mean_contact_fraction.png",
    "lipid_region_mean_matched_local_enrichment.png",
    "lipid_region_mean_matched_local_contact_fraction.png",
    "lipid_region_moiety_contact_fraction.png",
)
SUPPLEMENT_SUFFIXES = (
    "lipid_region_mean_local_enrichment.png",
    "__distribution.png", "__autocorrelation.png",
    "__bootstrap_table.png", "__comparison_table.png",
)


def _relative(path, root):
    """Persisted product identities are portable POSIX paths below the project."""
    return Path(path).relative_to(root).as_posix()


def figure_inventory(root, system_repouts, cfg=None, include_rejected=False, *, formats=("png", "pdf")):
    """List reporting products; formats=('png',) retains PNG-only consumers."""
    from filecmp import cmp
    root=Path(root)
    paths = [path for path in _figure_paths(root, system_repouts) if path.suffix[1:] in formats]
    import json
    provenance = Path(root)/'project_summary/figures/figure_provenance.json'
    records = json.loads(provenance.read_text()) if cfg is not None and provenance.exists() else {}
    rows=[]
    seen_png = {}
    for path in sorted(set(paths)):
        name=path.name
        kind=("primary" if name.startswith("lipid_contact_matrix_") or name.endswith(PRIMARY_SUFFIXES) else
              "table_report" if path.suffix == '.pdf' else "supplement" if name.endswith(SUPPLEMENT_SUFFIXES) else "other")
        status = figure_status(path,cfg,root,system_repouts,records=records) if cfg is not None else 'unverified'
        if path.suffix == '.png' and status in {'valid', 'legacy_unverified'}:
            peers = seen_png.setdefault(path.stat().st_size, [])
            if any(cmp(path, peer, shallow=False) for peer in peers):
                status = 'duplicate_figure'
            else:
                peers.append(path)
        if cfg is not None and status not in {'valid', 'legacy_unverified'} and not include_rejected:
            continue
        record = records.get(_relative(path, root),{})
        from .table_reports import report_manifest
        try:
            report = report_manifest(path) if path.suffix == '.pdf' else None
        except (OSError, ValueError, TypeError):
            report = None
        pages = report.get('pages', []) if report else []
        rows.append({"format":path.suffix[1:], "page_count":report.get('page_count') if report else None,
                     "page_ids":json.dumps([page.get('figure_id') for page in pages if isinstance(page, dict)]),
                     "figure_id":_relative(path, root),
                     "display_label":record.get('display_label', path.stem.replace('_',' ')),
                     "source_metric":record.get('source_metric'),"units":record.get('units'),
                     "generated_at":record.get('generated_at'),
                     "source_files":json.dumps(record.get('source_files',[])),
                     "sampling":json.dumps(record.get('sampling',[])),
                     "status":status,"family":figure_family(path),"figure":_relative(path, root),
                     "section":_relative(path.parent, root),
                     "kind":kind,"bytes":path.stat().st_size,"archive":None})
    if include_rejected and 'pdf' in formats:
        from .table_reports import MANIFEST_NAME, REPORT_NAME, report_manifest
        for sidecar in (root/'project_summary/figures').rglob(MANIFEST_NAME):
            path = sidecar.with_name(REPORT_NAME)
            if path.exists():
                continue
            try:
                report = report_manifest(path)
                status, archive = report.get('status', 'invalid_manifest'), report.get('archive')
            except (OSError, ValueError, TypeError):
                status, archive = 'invalid_manifest', None
            rows.append(dict(figure=_relative(path, root), figure_id=_relative(path, root),
                             format='pdf', page_count=0, page_ids='[]', bytes=0,
                             family=figure_family(path), status=status, kind='table_report',
                             section=_relative(path.parent, root), archive=archive))
    return pd.DataFrame(rows,columns=["figure","figure_id","display_label","section","kind","bytes","family","status","source_metric","units","generated_at","source_files","sampling","format","page_count","page_ids","archive"])


def figure_family(path):
    parts = Path(path).parts
    for token, family in (('shared_states','shared_states'),('sdc4','sdc4'),('sdc4_mechanisms','sdc4'),
                          ('solvation','solvation'),('diffusion','diffusion'),
                          ('lipid_contacts','lipid_contacts'),('membrane','lipid_enrichment'),
                          ('interactions','interactions'),('profiles','profiles'),
                          ('time_series','time_series')):
        if token in parts:
            return family
    return 'other'


def _figure_paths(root, system_repouts):
    paths = list((Path(root)/'project_summary').rglob('*.png'))
    for reps in system_repouts.values():
        for rep in reps:
            paths += list((Path(rep)/'extensions').rglob('*.png'))
            paths += list((Path(rep)/'figures').rglob('*.png'))
    paths += list((Path(root)/'project_summary').rglob('*.pdf'))
    return sorted({path for path in paths if not path.name.startswith('.')})


def figure_snapshot(root, system_repouts):
    return {_relative(p, root): (p.stat().st_mtime_ns,p.stat().st_size)
            for p in _figure_paths(root, system_repouts)}


def _source_files(root, system_repouts, family):
    root = Path(root)
    folders = {'time_series':'statistics/time_series','profiles':'statistics/profiles',
               'interactions':'statistics/interactions','solvation':'statistics/solvation',
               'shared_states':'shared_states','diffusion':'diffusion', 'lipid_contacts':'statistics/lipid_contacts'}
    paths = list((root/'project_summary'/folders[family]).glob('*.parquet')) if family in folders else []
    for reps in system_repouts.values():
        for rep in reps:
            rep = Path(rep)
            if family == 'interactions':
                paths += list((rep/'interactions').glob('*/occupancy.parquet'))
                paths += list((rep/'interactions').glob('*/presence.parquet'))
                paths += list((rep/'interactions').glob('*/lifetimes.parquet'))
            else:
                folder = {'sdc4':'sdc4_mechanisms','lipid_contacts':'lipid_contacts',
                          'solvation':'solvation','lipid_enrichment':'lipid_enrichment'}.get(family)
                if folder:
                    paths += list((rep/'extensions'/folder).glob('*.parquet'))
                    paths += list((rep/'extensions'/folder).glob('*.npz'))
                if family == 'profiles':
                    paths += list((rep/'cache').glob('*.npz'))
                if family in ('time_series','other'):
                    paths += list((rep/'summary').glob('*timeseries.parquet'))
    return sorted(set(paths))


def record_generated_figures(cfg, root, system_repouts, before):
    """Record source paths and sampling for figures written in this run."""
    import json
    from datetime import datetime

    from ..config import atomic_json
    root = Path(root)
    dest = root/'project_summary/figures/figure_provenance.json'
    records = json.loads(dest.read_text()) if dest.exists() else {}
    sampling = []
    for system, reps in system_repouts.items():
        for rep in reps:
            path = Path(rep) / "summary/timeseries.parquet"
            if path.is_file():
                times = pd.read_parquet(path, columns=["time_ps"])["time_ps"]
                sampling.append(dict(system=system, replicate=Path(rep).name,
                                     frames=len(times), first_ps=float(times.min()),
                                     last_ps=float(times.max())))
    sources = {}
    for path in _figure_paths(root,system_repouts):
        name = _relative(path, root)
        if before.get(name) == (path.stat().st_mtime_ns,path.stat().st_size):
            continue
        family = figure_family(path)
        if family not in sources:
            sources[family] = _source_files(root,system_repouts,family)
        files = sources[family]
        from .labels import _metric_meta
        if path.suffix == '.pdf':
            # A report spans metrics; each page carries its own figure identity.
            metric, units, label = None, None, f'Statistics tables ({path.parent.name})'
        else:
            metric = path.stem.split('__')[0]
            label, units = _metric_meta(metric)
            if units == metric:
                units = 'see figure axes/source schema'
        records[name] = dict(source_metric=metric,units=units,display_label=label,
            family=family, source_files=[_relative(p, root) for p in files],
            generated_at=datetime.now(UTC).isoformat(), sampling=sampling,
            format=path.suffix[1:])
        if path.suffix == '.pdf':
            from .table_reports import report_manifest
            report = report_manifest(path)
            if report is not None:
                records[name].update(page_count=report['page_count'], pages=report['pages'],
                                     sha256=report['sha256'],
                                     page_ids=[page['figure_id'] for page in report['pages']])
    atomic_json(dest,records)


def residue_response_status(cfg, root):
    """Keep previously rendered exploratory regions tied to effective settings."""
    import json

    from ..config import residue_response_settings
    settings = residue_response_settings(cfg)
    if (not settings["enabled"] or settings["profiles"] == [] or settings["chains"] == []):
        return "disabled"
    path = Path(root) / "project_summary/statistics/profiles/residue_response_manifest.json"
    if not path.exists():
        return "legacy_unverified"
    try:
        record = json.loads(path.read_text())
        if not isinstance(record, dict):
            return "invalid_manifest"
        if record.get("status") != "complete":
            return "incomplete"
        return "valid" if record.get("settings") == settings else "stale_configuration"
    except (OSError, ValueError, TypeError):
        return "invalid_manifest"


def figure_status(path, cfg, root, system_repouts, records=None):
    import json

    from PIL import Image

    from ..config import figure_enabled, plot_enabled
    from .plot_ids import id_from_figure_path
    root, path = Path(root), Path(path)
    family = figure_family(path)
    if path.parent.name == "residue_response" and "profiles" in path.parts:
        status = residue_response_status(cfg, root)
        if status not in {"valid", "legacy_unverified"}:
            return status
    # Old SDC4 diagnostics lived beside each replica's Parquet files.
    # The current figure renderer writes condition-level PNGs instead. Do not
    # advertise stale legacy PNGs as active, and never delete user artifacts.
    if family == "sdc4" and "extensions" in path.parts and "sdc4_mechanisms" in path.parts:
        return 'disabled'
    if not plot_enabled(cfg, family):
        return 'disabled'
    if path.suffix == '.pdf':
        from .table_reports import pdf_digest, report_manifest, table_page_enabled
        if not cfg.get('plots', {}).get('table_reports_pdf', False):
            return 'disabled'
        try:
            report = report_manifest(path)
        except (OSError, ValueError, TypeError):
            return 'invalid_manifest'
        try:
            digest = pdf_digest(path)
            if report is None:
                return 'legacy_unverified' if cfg.get('plots', {}).get('include') is None else 'unknown_identity'
            if report.get('status') != 'complete':
                return report.get('status', 'invalid_manifest')
            pages = report['pages']
            if (not pages or report.get('page_count') != len(pages) or
                    [page.get('page') for page in pages] != list(range(1, len(pages)+1)) or
                    any(not isinstance(page.get('figure_id'), str) for page in pages)):
                return 'invalid_manifest'
            if any(not table_page_enabled(cfg, page['figure_id']) for page in pages):
                return 'disabled'
            if digest != report.get('sha256'):
                return 'invalid_figure'
        except (OSError, ValueError, TypeError, KeyError, AttributeError):
            return 'invalid_figure'
    if path.suffix != '.pdf' and cfg.get("plots", {}).get("include") is not None:
        identifier = id_from_figure_path(path)
        if identifier is None or not figure_enabled(cfg, identifier):
            return 'disabled'
    parts = path.parts
    if family == 'interactions':
        kind = parts[parts.index('interactions')+1]
        if kind in cfg.get('plots', {}).get('interaction_exclude_kinds', []):
            return 'disabled'
    if records is None:
        provenance = root/'project_summary/figures/figure_provenance.json'
        records = json.loads(provenance.read_text()) if provenance.exists() else {}
    if path.suffix != '.pdf':
        try:
            with Image.open(path, formats=["PNG"]) as image:
                image.verify()
        except (OSError, SyntaxError, ValueError):
            return 'invalid_figure'
    record = records.get(_relative(path, root))
    if not record:
        return 'legacy_unverified'
    sources = record.get('source_files', [])
    if not isinstance(sources, list) or any(not isinstance(source, str) for source in sources):
        return 'invalid_manifest'
    if sources and any(not (root / source).is_file() for source in sources):
        return 'missing_source'
    return 'valid' if sources else 'legacy_unverified'


def statistics_status(cfg, root, family):
    """Explain eligibility with explicit legacy and publication policies."""
    from ..statistics.manifest import read_statistics_manifest
    stats = cfg.get('statistics', {})
    option = stats.get('families', {}).get(family, True)
    if (not stats.get('enabled', True) or option is False or
            (isinstance(option, dict) and not option.get('enabled', True))):
        return 'disabled'
    facts = read_statistics_manifest(Path(root)/'project_summary/statistics/manifest.json')
    if facts.state == 'missing':
        return 'legacy_unverified'
    if facts.state == 'invalid':
        return 'invalid_manifest'
    if facts.state == 'incomplete':
        return 'incomplete'
    if facts.family_disabled(family):
        return 'unavailable'
    if facts.state == 'legacy' or not facts.methods_present:
        return 'legacy_unverified'
    return 'valid'


def current_statistics(cfg, root, family):
    """Do not ingest disabled, incomplete, or malformed statistical products."""
    return statistics_status(cfg, root, family) in {'valid', 'legacy_unverified'}
