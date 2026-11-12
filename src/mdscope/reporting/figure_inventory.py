"""Inventory and provenance for manually generated scientific figures."""
from __future__ import annotations

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


def figure_inventory(root, system_repouts, cfg=None, include_rejected=False):
    """List every generated PNG with its origin and preferred presentation."""
    from filecmp import cmp
    root=Path(root)
    paths = _figure_paths(root, system_repouts)
    import json
    provenance = Path(root)/'project_summary/report/figure_provenance.json'
    records = json.loads(provenance.read_text()) if cfg is not None and provenance.exists() else {}
    rows=[]
    seen_png = {}
    for path in sorted(set(paths)):
        name=path.name
        kind=("primary" if name.startswith("lipid_contact_matrix_") or name.endswith(PRIMARY_SUFFIXES) else
              "supplement" if name.endswith(SUPPLEMENT_SUFFIXES) else "other")
        status = figure_status(path,cfg,root,system_repouts,records=records) if cfg is not None else 'unverified'
        if status in {'valid', 'legacy_unverified'}:
            peers = seen_png.setdefault(path.stat().st_size, [])
            if any(cmp(path, peer, shallow=False) for peer in peers):
                status = 'duplicate_figure'
            else:
                peers.append(path)
        if cfg is not None and status not in {'valid', 'legacy_unverified'} and not include_rejected:
            continue
        record = records.get(str(path.relative_to(root)),{})
        rows.append({"figure_id":str(path.relative_to(root)),
                     "display_label":record.get('display_label', path.stem.replace('_',' ')),
                     "source_metric":record.get('source_metric'),"units":record.get('units'),
                     "generated_at":record.get('generated_at'),
                     "source_files":json.dumps(record.get('source_files',[])),
                     "sampling":json.dumps(record.get('sampling',[])),
                     "status":status,"family":figure_family(path),"figure":str(path.relative_to(root)),
                     "section":str(path.parent.relative_to(root)),
                     "kind":kind,"bytes":path.stat().st_size})
    return pd.DataFrame(rows,columns=["figure","figure_id","display_label","section","kind","bytes","family","status","source_metric","units","generated_at","source_files","sampling"])


def figure_family(path):
    parts = Path(path).parts
    for token, family in (('shared_states','shared_states'),('sdc4_mechanisms','sdc4'),
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
    return sorted(set(paths))


def figure_snapshot(root, system_repouts):
    return {str(p.relative_to(root)): (p.stat().st_mtime_ns,p.stat().st_size)
            for p in _figure_paths(root, system_repouts)}


def _source_files(root, system_repouts, family):
    root = Path(root)
    folders = {'time_series':'statistics/time_series','profiles':'statistics/profiles',
               'interactions':'statistics/interactions','solvation':'statistics/solvation',
               'shared_states':'shared_states','diffusion':'diffusion'}
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
    from datetime import datetime, timezone
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
        name = str(path.relative_to(root))
        if before.get(name) == (path.stat().st_mtime_ns,path.stat().st_size):
            continue
        family = figure_family(path)
        if family not in sources:
            sources[family] = _source_files(root,system_repouts,family)
        files = sources[family]
        from .labels import _metric_meta
        metric = path.stem.split('__')[0]
        label, units = _metric_meta(metric)
        if units == metric:
            units = 'see figure axes/source schema'
        records[name] = dict(source_metric=metric,units=units,display_label=label,
            family=family, source_files=[str(p.relative_to(root)) for p in files],
            generated_at=datetime.now(timezone.utc).isoformat(), sampling=sampling)
    atomic_json(dest,records)


def figure_status(path, cfg, root, system_repouts, records=None):
    import json
    from PIL import Image
    from ..config import plot_enabled
    root, path = Path(root), Path(path)
    family = figure_family(path)
    if not plot_enabled(cfg, family):
        return 'disabled'
    parts = path.parts
    if family == 'interactions':
        kind = parts[parts.index('interactions')+1]
        if kind in cfg.get('plots', {}).get('interaction_exclude_kinds', []):
            return 'disabled'
    if records is None:
        provenance = root/'project_summary/figures/figure_provenance.json'
        records = json.loads(provenance.read_text()) if provenance.exists() else {}
    try:
        with Image.open(path, formats=["PNG"]) as image:
            image.verify()
    except (OSError, SyntaxError, ValueError):
        return 'invalid_figure'
    record = records.get(str(path.relative_to(root)))
    if not record:
        return 'legacy_unverified'
    sources = record.get('source_files', [])
    if not isinstance(sources, list) or any(not isinstance(source, str) for source in sources):
        return 'invalid_manifest'
    if sources and any(not (root / source).is_file() for source in sources):
        return 'missing_source'
    return 'valid' if sources else 'legacy_unverified'


def statistics_status(cfg, root, family):
    """Explain eligibility, preserving explicitly identified legacy exports."""
    import json
    stats = cfg.get('statistics', {})
    option = stats.get('families', {}).get(family, True)
    if not stats.get('enabled', True) or option is False or (isinstance(option, dict) and not option.get('enabled', True)):
        return 'disabled'
    manifest = Path(root)/'project_summary/statistics/manifest.json'
    if not manifest.exists():
        return 'legacy_unverified'
    try:
        record = json.loads(manifest.read_text())
        if record.get('status', 'complete') != 'complete':
            return 'incomplete'
        if record.get('families', {}).get(family) == 'disabled':
            return 'unavailable'
        if 'methods' not in record:
            return 'legacy_unverified'
    except (ValueError, TypeError, KeyError, AttributeError, OSError):
        return 'invalid_manifest'
    return 'valid'


def current_statistics(cfg, root, family):
    """Do not ingest disabled, incomplete, or malformed statistical products."""
    return statistics_status(cfg, root, family) in {'valid', 'legacy_unverified'}
