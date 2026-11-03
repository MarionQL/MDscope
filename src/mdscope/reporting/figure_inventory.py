"""Inventory generated scientific figures and audit coverage of reportable analyses.

Only persisted artifacts are inspected; never reopen trajectories. Primary
figures are selected deterministically, with all figures preserved in the ZIP.
"""
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
    from ..config import config_hash
    root=Path(root)
    paths = _figure_paths(root, system_repouts)
    fingerprints = {}
    signature = _presentation_hash(cfg) if cfg is not None else None
    import json
    provenance = Path(root)/'project_summary/report/figure_provenance.json'
    records = json.loads(provenance.read_text()) if cfg is not None and provenance.exists() else {}
    rows=[]
    seen_png = set()
    for path in sorted(set(paths)):
        name=path.name
        kind=("primary" if name.startswith("lipid_contact_matrix_") or name.endswith(PRIMARY_SUFFIXES) else
              "supplement" if name.endswith(SUPPLEMENT_SUFFIXES) else "other")
        status = figure_status(path,cfg,root,system_repouts,records=records,fingerprints=fingerprints,signature=signature) if cfg is not None else 'unverified'
        if status == 'valid':
            digest = records[str(path.relative_to(root))].get('png_sha256')
            if digest in seen_png:
                status = 'duplicate_figure'
            seen_png.add(digest)
        if cfg is not None and status != 'valid' and not include_rejected:
            continue
        record = records.get(str(path.relative_to(root)),{})
        rows.append({"figure_id":config_hash(str(path.relative_to(root))),
                     "display_label":record.get('display_label', path.stem.replace('_',' ')),
                     "source_metric":record.get('source_metric'),"units":record.get('units'),
                     "generation_config_hash":record.get('generation_config_hash'),
                     "generated_at":record.get('generated_at'),
                     "source_files":json.dumps(record.get('source_files',[])),
                     "sampling":json.dumps(record.get('sampling',[])),
                     "status":status,"family":figure_family(path),"figure":str(path.relative_to(root)),
                     "section":str(path.parent.relative_to(root)),
                     "kind":kind,"bytes":path.stat().st_size})
    return pd.DataFrame(rows,columns=["figure","figure_id","display_label","section","kind","bytes","family","status","source_metric","units","generation_config_hash","generated_at","source_files","sampling"])


def report_coverage(root, system_repouts, cfg=None):
    """Audit which statistical analysis families have corresponding figures.

    This checks family-level presence, not scientific correctness or the
    completeness of every plotted metric.
    """
    root=Path(root)
    fig=figure_inventory(root,system_repouts,cfg)
    paths=fig.figure.tolist()
    rows=[]
    specs=(
        ("Time series","statistics/time_series/system_statistics.parquet",
         "project_summary/figures/time_series/"),
        ("Profiles","statistics/profiles/system_profile_statistics.parquet",
         "project_summary/figures/profiles/"),
        ("Interactions","statistics/interactions/",
         "project_summary/figures/interactions/"),
    )
    for family,source,visual in specs:
        path=root/"project_summary"/source
        has_data=(path.is_dir() and any(path.rglob("*.parquet"))) if source.endswith("/") else path.exists()
        has_figure=any(p.startswith(visual) for p in paths)
        rows.append({"analysis":family,"data_available":bool(has_data),
                     "figure_available":bool(has_figure),
                     "status":"covered" if has_data and has_figure else
                              "missing_figure" if has_data else "not_generated"})
    for family,filename in (
        ("Solvation","species_timeseries.parquet"),
        ("SDC4 mechanisms","timeseries.parquet"),
    ):
        folder="solvation" if family=="Solvation" else "sdc4_mechanisms"
        found=any((Path(rep)/"extensions"/folder/filename).exists()
                  for reps in system_repouts.values() for rep in reps)
        has_figure=any("solvation" in p for p in paths) if family=="Solvation" else any(
            "sdc4" in p or "lipid_" in p or "shared_states" in p for p in paths)
        rows.append({"analysis":family,"data_available":found,
                     "figure_available":has_figure,
                     "status":"covered" if found and has_figure else
                              "missing_figure" if found else "not_generated"})
    result = pd.DataFrame(rows)
    if cfg is not None:
        from ..config import presentation_enabled
        names = {'Time series':'time_series','Profiles':'profiles','Interactions':'interactions',
                 'Solvation':'solvation','SDC4 mechanisms':'sdc4'}
        rejected = figure_inventory(root,system_repouts,cfg,include_rejected=True)
        for index, row in result.iterrows():
            family = names[row.analysis]
            enabled = presentation_enabled(cfg,'plots',family) and presentation_enabled(cfg,'report',family)
            result.loc[index,'expected'] = bool(enabled and (row.data_available or cfg.get('analyses',{}).get({'sdc4':'sdc4_mechanisms'}.get(family,family),{}).get('enabled',False)))
            if not enabled:
                result.loc[index,'status'] = 'disabled'
            elif family in {'time_series','profiles','interactions'} and row.data_available and not current_statistics(cfg,root,family):
                result.loc[index,'status'] = 'invalid_data'
            elif row.data_available and not row.figure_available and not rejected.loc[rejected.family.eq(family)].empty:
                result.loc[index,'status'] = 'invalid_figure'
            elif result.loc[index,'expected'] and not row.data_available:
                result.loc[index,'status'] = 'missing_data'
    return result


def metric_coverage(root, cfg=None, system_repouts=None):
    """Check every scalar statistic has a generated standard time-series plot.

    This is exact per (group, metric), not just a family-level presence test.
    It intentionally reports missing plots instead of silently dropping data.
    """
    root=Path(root)
    source=root/"project_summary"/"statistics"/"time_series"/"system_statistics.parquet"
    cols=["group","metric","figure","status"]
    if not source.exists():
        return pd.DataFrame(columns=cols)
    table=pd.read_parquet(source,columns=["group","metric"]).drop_duplicates()
    dest=root/"project_summary"/"figures"/"time_series"
    records=[]
    eligible = set(figure_inventory(root,system_repouts or {},cfg).figure) if cfg is not None else None
    for item in table.itertuples(index=False):
        group="" if pd.isna(item.group) else str(item.group)
        metric=str(item.metric)
        stem=metric+(f"__group_{group}" if group else "")
        stem="".join(c if c.isalnum() or c in "-_." else "_" for c in stem).strip("_")
        relative=f"project_summary/figures/time_series/{stem}__timeseries.png"
        records.append({"group":group,"metric":metric,"figure":relative,
                        "status":("covered" if relative in eligible else "missing_figure") if eligible is not None
                        else "covered" if (dest/(stem+"__timeseries.png")).exists() else "missing_figure"})
    result = pd.DataFrame(records,columns=cols)
    if cfg is not None and not current_statistics(cfg,root,'time_series'):
        result['status'] = 'excluded_statistics'
    return result



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


def _presentation_hash(cfg):
    from ..config import config_hash
    from ..core.provenance import package_code_hash
    return config_hash({'config': {k: cfg.get(k) for k in ('plots','statistics','analyses','trajectory','time','selections','systems')},
                        'code_hash': package_code_hash()})


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


def _source_fingerprint(paths, root):
    from ..config import config_hash
    return config_hash([(str(p.relative_to(root)),p.stat().st_size,p.stat().st_mtime_ns) for p in paths]) if paths else None


def record_generated_figures(cfg, root, system_repouts, before):
    """Stamp only newly written PNGs; never bless an unchanged stale figure.

    Source fingerprints use path/size/mtime, not a claimed content digest.
    PNG content hashes detect edited/replaced images. No trajectory is read.
    """
    import json, hashlib
    from datetime import datetime, timezone
    from ..config import atomic_json
    from .report_inventory import _summary_inventory
    root = Path(root)
    dest = root/'project_summary/report/figure_provenance.json'
    records = json.loads(dest.read_text()) if dest.exists() else {}
    sampling = _summary_inventory(root,cfg,system_repouts).to_dict('records')
    sources = {}
    signature = _presentation_hash(cfg)
    for path in _figure_paths(root,system_repouts):
        name = str(path.relative_to(root))
        if before.get(name) == (path.stat().st_mtime_ns,path.stat().st_size):
            continue
        family = figure_family(path)
        if family not in sources:
            files = _source_files(root,system_repouts,family)
            sources[family] = (files,_source_fingerprint(files,root))
        files, fingerprint = sources[family]
        from .labels import _metric_meta
        metric = path.stem.split('__')[0]
        label, units = _metric_meta(metric)
        if units == metric:
            units = 'see figure axes/source schema'
        records[name] = dict(source_metric=metric,units=units,display_label=label,
            family=family, generation_config_hash=signature,
            png_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            source_fingerprint=fingerprint, source_files=[str(p.relative_to(root)) for p in files],
            generated_at=datetime.now(timezone.utc).isoformat(), sampling=sampling)
    atomic_json(dest,records)


def figure_status(path, cfg, root, system_repouts, records=None, fingerprints=None, signature=None):
    import json, hashlib
    from ..config import presentation_enabled
    root, path = Path(root), Path(path)
    family = figure_family(path)
    if not presentation_enabled(cfg,'plots',family) or not presentation_enabled(cfg,'report',family):
        return 'disabled'
    parts = path.parts
    if family == 'interactions':
        kind = parts[parts.index('interactions')+1]
        if any(kind in cfg.get(domain,{}).get('interaction_exclude_kinds',[]) for domain in ('plots','report')):
            return 'disabled'
    if records is None:
        provenance = root/'project_summary/report/figure_provenance.json'
        records = json.loads(provenance.read_text()) if provenance.exists() else {}
    record = records.get(str(path.relative_to(root)))
    if record is None:
        return 'unverified'
    if not all(key in record for key in ('generation_config_hash','source_fingerprint','png_sha256')):
        return 'invalid_manifest'
    if record['generation_config_hash'] != (signature or _presentation_hash(cfg)):
        return 'stale_config'
    if fingerprints is None:
        fingerprints = {}
    if family not in fingerprints:
        fingerprints[family] = _source_fingerprint(_source_files(root,system_repouts,family),root)
    if fingerprints[family] is None:
        return 'missing_source'
    if record['source_fingerprint'] != fingerprints[family]:
        return 'stale_source'
    if record['png_sha256'] != hashlib.sha256(path.read_bytes()).hexdigest():
        return 'changed_figure'
    return 'valid'


def statistics_status(cfg, root, family):
    """Explain eligibility, preserving explicitly identified legacy exports."""
    import json
    from ..config import config_hash
    from ..statistics import analysis_config_hash, input_fingerprint
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
        if record.get('statistics_hash') != config_hash(stats):
            return 'stale_config'
        if 'analysis_config_hash' in record and record['analysis_config_hash'] != analysis_config_hash(cfg):
            return 'stale_config'
        if 'input_fingerprint' in record:
            if record['input_fingerprint'] != input_fingerprint(root, record['replicate_paths']):
                return 'stale_source'
        else:
            return 'legacy_unverified'
    except (ValueError, TypeError, KeyError, AttributeError, OSError):
        return 'invalid_manifest'
    return 'valid'


def current_statistics(cfg, root, family):
    """Do not ingest disabled, incomplete, or mismatched statistical products."""
    return statistics_status(cfg, root, family) in {'valid', 'legacy_unverified'}
