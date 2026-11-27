"""Paginated analysis overviews from saved series and saved contrast intervals.

Drawing never estimates uncertainty. Frame distributions stay replica-specific;
contrast intervals are copied verbatim from the statistics producer.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import figure_enabled, merge, replicate_label, system_label
from ..config.defaults import DEFAULTS
from ..registry import scalar_contracts, statistical_variable_owner
from .figure_inventory import current_statistics
from .plot_ids import figure_dpi, figure_path, save_figure


def _save_page(cfg, root, identifier, fig, page):
    path = figure_path(cfg, Path(root) / 'project_summary/figures/overview.png', identifier=identifier)
    path = path.with_stem(path.stem + f'@page_{page:03d}')
    fig.tight_layout(rect=(0, .04, 1, .94))
    save_figure(fig, cfg, path, dpi=figure_dpi(cfg, identifier, cfg.get('plots', {}).get('dpi', 180)),
                bbox_inches='tight')
    plt.close(fig)


def draw_replica_distributions(ax, series, *, unit):
    """One violin per independently preserved replica; diamonds are means."""
    labels = []
    for position, (label, values) in enumerate(series, 1):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        labels.append(label + f' (n={len(values)})')
        if not len(values):
            ax.text(position, .5, 'Missing', rotation=90, ha='center',
                    transform=ax.get_xaxis_transform(), color='.4')
        else:
            if len(values) > 1 and np.ptp(values) > 0:
                ax.violinplot([values], positions=[position], showextrema=False)
            ax.scatter([position], [values.mean()], marker='D', color='black', s=18)
    ax.set_xticks(range(1, len(labels) + 1), labels, rotation=35, ha='right', fontsize=8)
    ax.set_ylabel(unit)


def draw_effect_intervals(ax, rows, *, unit):
    """Saved signed differences; invalid or unavailable bounds stay unavailable."""
    labels = []
    for position, row in enumerate(rows):
        labels.append(row['label'])
        value, low, high = (row[key] for key in ('effect', 'low', 'high'))
        valid = np.isfinite([low, high]).all() and low <= high
        excludes = valid and (low > 0 or high < 0)
        color = '#9c3d10' if excludes else '#245a81'
        if valid:
            ax.plot([low, high], [position, position], color=color, linewidth=2)
        if np.isfinite(value):
            ax.scatter(value, position, color=color, marker='D' if excludes else 'o', s=26)
        if not valid:
            ax.text(.99, position, 'CI unavailable', transform=ax.get_yaxis_transform(),
                    ha='right', fontsize=8, color='.4')
    ax.axvline(0, color='.5', linewidth=.8)
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    limits = ax.get_xlim()
    extent = max(abs(limits[0]), abs(limits[1]))
    ax.set_xlim(-extent, extent)
    ax.set_xlabel('System 2 − system 1 (' + unit + ')')


def _distributions(cfg, root, system_repouts):
    values = {}
    requested = [item for item in scalar_contracts(cfg, include_figures=False)
                 if figure_enabled(cfg, f"overview:{item['analysis']}:distributions")]
    for system, reps in system_repouts.items():
        for rep in reps:
            tables = {}
            for item in requested:
                path = Path(rep) / item['path_within_replica']
                if not path.is_file():
                    continue
                if path not in tables:
                    tables[path] = pd.read_parquet(path)
                data = tables[path]
                group, metric = item.get('group'), item['metric']
                if group is not None:
                    if 'group' not in data:
                        continue
                    data = data.loc[data.group == group]
                if metric not in data:
                    continue
                values.setdefault((item['analysis'], item['units'][metric], metric, str(group or 'global')), {})[
                    (system, Path(rep).name)] = pd.to_numeric(data[metric], errors='coerce').to_numpy(float)
    # Saved dynamic variables use the same ownership resolver as table reports.
    # Only requested extension sources are opened; raw identities remain exact.
    from ..statistics.time_series import collect_time_series
    from ..contracts.outputs import FUSED_OUTPUTS
    from ..registry import STRUCTURAL_SCALAR_METRICS
    registered_identities = {(item.get("group"), item["metric"]) for item in requested}
    for system, reps in system_repouts.items():
        for rep in reps:
            sources = {'summary'}
            for source in ('solvation', 'hydration', 'lipid_contacts'):
                if figure_enabled(cfg, f'overview:{source}:distributions'):
                    sources.add(source)
            if 'solvation' in sources:
                sources.add('solvation_bridges')
            series = collect_time_series(rep, sources=sources)
            for owner in ('sdc4_mechanisms', 'local_environment'):
                if not figure_enabled(cfg, f'overview:{owner}:distributions'):
                    continue
                path = Path(rep) / FUSED_OUTPUTS[owner]
                if path.is_file():
                    data = pd.read_parquet(path)
                    for metric in data.columns:
                        if metric != 'time_ps' and pd.api.types.is_numeric_dtype(data[metric]):
                            series[(None, metric)] = pd.DataFrame({'value': data[metric]})
            for (group, metric), data in series.items():
                if metric in STRUCTURAL_SCALAR_METRICS.values() and (group, metric) not in registered_identities:
                    continue
                owner, unit = statistical_variable_owner(cfg, 'time_series', {'metric': metric, 'group': group})
                if owner == 'unmapped' or not figure_enabled(cfg, f'overview:{owner}:distributions'):
                    continue
                if 'see metric' in unit or 'producer-defined' in unit:
                    from .labels import _metric_meta
                    _, label_unit = _metric_meta(metric)
                    unit = label_unit if label_unit and label_unit != metric else 'source units unspecified'
                values.setdefault((owner, unit, metric, str(group or 'global')), {})[(system, Path(rep).name)] = data.value.to_numpy(float)
    by_owner = {}
    for (owner, unit, metric, group), replicas in sorted(values.items()):
        series = [(f'{system_label(cfg, system)}/{replicate_label(cfg, system, Path(rep).name)}', replicas.get((system, Path(rep).name), []))
                  for system, reps in system_repouts.items() for rep in reps]
        by_owner.setdefault(owner, []).append((unit, metric, group, series))
    for owner, panels in by_owner.items():
        bounds = {}
        for unit, metric, group, series in panels:
            finite = np.concatenate([np.asarray(values, float)[np.isfinite(values)] for _, values in series])
            if len(finite):
                low, high = float(finite.min()), float(finite.max())
                previous = bounds.get((unit, metric), (low, high))
                bounds[(unit, metric)] = (min(low, previous[0]), max(high, previous[1]))
        # Four metric/group facets per page; no mixing of raw scales or units.
        for offset in range(0, len(panels), 4):
            chunk = panels[offset:offset+4]
            fig, axes = plt.subplots(len(chunk), 1, figsize=(max(9, .8 * max(len(p[3]) for p in chunk)), 4.4 * len(chunk)), squeeze=False)
            for ax, (unit, metric, group, series) in zip(axes.flat, chunk):
                draw_replica_distributions(ax, series, unit=unit)
                ax.set_title(f'{metric} · {group}', fontsize=11)
                if (unit, metric) in bounds:
                    low, high = bounds[(unit, metric)]
                    padding = .05 * (high-low) if high > low else max(.05, abs(low)*.05)
                    ax.set_ylim(low-padding, high+padding)
            fig.suptitle(owner.replace('_', ' ') + ' · replica distributions')
            fig.text(.5, .01, 'Each violin weights saved finite frames equally within one replica; diamonds = replica means. No inferential CI.', ha='center', fontsize=8)
            _save_page(cfg, root, f'overview:{owner}:distributions', fig, offset//4+1)


def _effects(cfg, root):
    source = Path(root) / 'project_summary/statistics/time_series/system_comparisons.parquet'
    if not source.is_file() or not current_statistics(cfg, root, 'time_series'):
        return
    data = pd.read_parquet(source)
    columns = ['mean_difference_system2_minus_system1', 'mean_difference_ci_low', 'mean_difference_ci_high']
    if not {'metric', 'system_1', 'system_2', *columns} <= set(data):
        return
    facets = {}
    for row in data.to_dict('records'):
        if pd.isna(row.get('group')):
            row['group'] = None
        owner, unit = statistical_variable_owner(cfg, 'time_series', row)
        if owner == 'unmapped' or not figure_enabled(cfg, f'overview:{owner}:effects'):
            continue
        # A namespace-level unit note does not prove physical or normalization
        # compatibility. Isolate each such estimand, even if display suffixes
        # happen to resolve to the same dimension.
        metric_facet = ''
        if unit == 'unspecified' or 'see metric' in unit or 'producer-defined' in unit:
            metric_facet = str(row['metric'])
            from .labels import _metric_meta
            _, display_unit = _metric_meta(metric_facet)
            unit = display_unit if display_unit and display_unit != metric_facet else 'source units unspecified'
        level = row.get('ci_level')
        level_label = f'{100*float(level):g}%' if level is not None and pd.notna(level) else 'unspecified level'
        key = (owner, unit, row['system_1'], row['system_2'], level_label, metric_facet)
        group = row.get('group')
        group_label = str(group) if pd.notna(group) else 'global'
        facets.setdefault(key, []).append(dict(label=str(row['metric']) + ' · ' + group_label,
                                                **dict(zip(('effect', 'low', 'high'), [float(row[c]) for c in columns]))))
    pages = {}
    for (owner, unit, first, second, level, metric_facet), rows in sorted(facets.items()):
        rows.sort(key=lambda row: row['label'])
        for offset in range(0, len(rows), 24):
            chunk = rows[offset:offset+24]
            fig, ax = plt.subplots(figsize=(12, max(4, 2 + .38 * len(chunk))))
            draw_effect_intervals(ax, chunk, unit=unit)
            ax.set_title(f'{owner}: {second} − {first}\nSaved {level} bootstrap intervals; diamonds: CI excludes zero (not an FDR decision)', fontsize=11)
            pages[owner] = pages.get(owner, 0) + 1
            _save_page(cfg, root, f'overview:{owner}:effects', fig, pages[owner])


def plot_overviews(cfg, root, system_repouts):
    """Opt-in overview IDs; legacy omitted selectors do not gain new figures."""
    options = cfg.get('plots', {})
    if options.get('include') is None and options.get('presets') is None:
        return
    config = merge(DEFAULTS, cfg)
    selected = options.get('include')
    if selected is None or any(identifier.startswith('overview:') and identifier.endswith(':distributions') for identifier in selected):
        _distributions(config, root, system_repouts)
    if selected is None or any(identifier.startswith('overview:') and identifier.endswith(':effects') for identifier in selected):
        _effects(config, root)
