"""Analysis-owned, saved-statistics tables; inference remains in statistics/.

Selection is applied to source identities before consolidation. CSV retains every
source column; PDFs paginate rows and columns without ranking or truncation.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import plot_enabled
from ..contracts.statistics import STATISTICS_CANDIDATES
from ..registry import statistical_variable_owner
from .figure_inventory import current_statistics, statistics_status
from ..contracts.plot_validation import safe_plot_stem
from .table_reports import (
    REPORT_NAME,
    add_report_data,
    add_table_page,
    file_digest,
    table_page_enabled,
    table_selection,
)


def _source_id(family, row, variant):
    if family == 'profiles':
        stem = safe_plot_stem(row['profile'])
    elif family == 'time_series':
        group = row.get('group')
        stem = safe_plot_stem(str(row['metric']) + (f'__group_{group}' if pd.notna(group) else ''))
    else:
        stem = f"{safe_plot_stem(row['interaction_type'])}:{safe_plot_stem(row['metric'])}"
    return f'{family}:{stem}:{variant}'


def _condition_join(comparisons, conditions, family):
    if conditions.empty or comparisons.empty:
        return comparisons
    keys = {'profiles': ['profile', 'position', 'position_label'],
            'time_series': ['metric', 'group'],
            'interactions': ['interaction_type', 'interaction', 'metric']}[family]
    if not set(keys) <= set(comparisons) or not {*keys, 'system'} <= set(conditions):
        raise ValueError(f'Cannot join {family} statistics: exact identity columns required: {keys}')
    result = comparisons
    for side in (1, 2):
        columns = {column: f'condition_{side}_{column}' for column in conditions if column not in keys}
        right = conditions.rename(columns=columns)
        result = result.merge(right, left_on=[*keys, f'system_{side}'],
                              right_on=[*keys, f'condition_{side}_system'],
                              how='left', validate='many_to_one', sort=False)
    return result



def _lipid_identity(profile):
    parts = profile.split('__')[1:]
    leaflet = parts.pop() if parts[-1] in {'upper', 'lower'} else 'combined'
    tail = parts.pop()
    return '__'.join(parts) or 'all lipids', tail, leaflet


def _profile_identity_columns(data):
    if 'profile' not in data:
        return data
    lipid = data.profile.str.startswith('lipid_order__')
    if lipid.any():
        identities = [_lipid_identity(profile) for profile in data.loc[lipid, 'profile']]
        for i, column in enumerate(('lipid', 'tail', 'leaflet')):
            data.loc[lipid, column] = [identity[i] for identity in identities]
    if 'position_label' in data:
        residues = data.position_label.astype('string').str.extract(
            r'^(?P<chain>[^:]+):(?P<residue_name>[A-Z]+)(?P<residue_number>-?[0-9]+)$')
        for column in residues:
            if residues[column].notna().any():
                data[column] = residues[column]
    return data

def collect_statistical_reports(cfg, root):
    """Return selected full-precision rows, with exact-key condition joins."""
    root = Path(root)
    frames = []
    for family in ('time_series', 'profiles', 'interactions'):
        if not plot_enabled(cfg, family):
            continue
        selected = cfg.get('plots', {}).get('include')
        if selected is not None and not any(identifier.startswith(f'{family}:') and
                identifier.endswith((':bootstrap_table', ':comparison_table')) and
                table_page_enabled(cfg, identifier) for identifier in selected):
            continue
        directory = root / 'project_summary/statistics' / family
        # Full lag curves/episodes are separate products, not bootstrap summary
        # rows. Their integrated diagnostics already live in replica statistics.
        auxiliary = {'autocorrelation.parquet', 'hs_contact_episodes.parquet', 'residue_regions.parquet'}
        paths = [directory / name for name in STATISTICS_CANDIDATES[family]
                 if name not in auxiliary and (directory / name).is_file()]
        if not current_statistics(cfg, root, family):
            if paths and statistics_status(cfg, root, family) not in {'disabled', 'unavailable'}:
                raise ValueError(f'Cannot publish tables from {statistics_status(cfg, root, family)} statistics: {directory}')
            continue
        wants_comparisons = selected is None or any(identifier.startswith(f'{family}:') and
                            identifier.endswith(':comparison_table') for identifier in selected)
        wants_bootstrap = cfg.get('plots', {}).get('bootstrap_tables', True) and (selected is None or any(
                            identifier.startswith(f'{family}:') and identifier.endswith(':bootstrap_table') for identifier in selected))
        tables = {path.name: pd.read_parquet(path) for path in paths
                  if (wants_comparisons if 'comparisons' in path.name else
                      wants_bootstrap or (wants_comparisons and path.name in
                          {'system_profiles.parquet', 'system_statistics.parquet', 'system_interaction_metrics.parquet'}))}
        condition_name = {'profiles': 'system_profiles.parquet',
                          'time_series': 'system_statistics.parquet',
                          'interactions': 'system_interaction_metrics.parquet'}[family]
        conditions = tables.get(condition_name, pd.DataFrame())
        required = {'profiles': {'profile'}, 'time_series': {'metric'},
                    'interactions': {'interaction_type', 'metric', 'interaction'}}[family]
        for filename, data in tables.items():
            if data.empty:
                continue
            if not required <= set(data):
                raise ValueError(f'Malformed saved statistics {directory / filename}: missing identity columns {sorted(required - set(data))}')
            comparison = 'comparisons' in filename
            variant = 'comparison_table' if comparison else 'bootstrap_table'
            data = data.copy()
            data['source_id'] = [_source_id(family, row, variant) for row in data.to_dict('records')]
            data = data.loc[data.source_id.map(lambda identifier: table_page_enabled(cfg, identifier))].copy()
            if family == 'interactions':
                data = data.loc[~data.interaction_type.isin(cfg.get('plots', {}).get('interaction_exclude_kinds', []))]
            if data.empty:
                continue
            if comparison:
                data = _condition_join(data, conditions, family)
                if not conditions.empty:
                    data['condition_source_table'] = (directory / condition_name).relative_to(root).as_posix()
            identity_columns = [column for column in ('profile', 'metric', 'group', 'interaction_type') if column in data]
            ownership = {tuple(None if pd.isna(value) else value for value in row.values()): statistical_variable_owner(cfg, family, row)
                         for row in data[identity_columns].drop_duplicates().to_dict('records')}
            owners = [ownership[tuple(None if pd.isna(value) else value for value in row.values())] for row in data[identity_columns].to_dict('records')]
            data = _profile_identity_columns(data)
            data['analysis'] = [item[0] for item in owners]
            data['units'] = [item[1] for item in owners]
            data['family'] = family
            data['section'] = ('comparison' if comparison else 'diagnostics' if 'autocorrelation' in filename or
                               'profile_statistics' in filename else 'replica' if filename.startswith('replicate') else 'condition')
            data['source_table'] = (directory / filename).relative_to(root).as_posix()
            if comparison:
                alpha = cfg.get('statistics', {}).get('hypothesis_tests', {}).get('alpha', .05)
                for prefix in ('difference', 'mean_difference', 'median_difference'):
                    low, high = f'{prefix}_ci_low', f'{prefix}_ci_high'
                    if low in data and high in data:
                        lo, hi = data[low], data[high]
                        data[f'{prefix}_ci_excludes_zero'] = (np.isfinite(lo) & np.isfinite(hi) &
                            (lo <= hi) & ((lo > 0) | (hi < 0)))
                for q in ('q_value', 'mean_q_value', 'median_q_value'):
                    if q in data:
                        data[f'{q}_below_alpha'] = np.isfinite(data[q]) & (data[q] < alpha)
                flag = 'mean_difference_ci_excludes_zero' if family == 'time_series' else 'difference_ci_excludes_zero'
                data['ci_excludes_zero'] = data.get(flag, False)
                q = 'mean_q_value_below_alpha' if family == 'time_series' else 'q_value_below_alpha'
                data['q_below_alpha'] = data.get(q, False)
            frames.append(data)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def _display(value):
    if pd.isna(value):
        return 'N/A'
    if isinstance(value, (float, np.floating)):
        return f'{value:.4g}'
    return str(value)


def _interval(low, high):
    return f'[{_display(low)}, {_display(high)}]'


def _pages(frame, title, path, identifiers, *, alpha=.05, columns_per_page=7):
    """Fixed paper, wrapped cells, repeated identities and column continuation."""
    frame = frame.copy()
    frame.insert(0, 'row', np.arange(1, len(frame) + 1))
    columns = [column for column in frame.columns[1:] if not column.endswith(('_excludes_zero', '_below_alpha'))]
    leading = [column for column in columns if column in {'profile', 'metric', 'group', 'interaction_type', 'interaction', 'position', 'position_label', 'system_1', 'system_2', 'system', 'replicate'}]
    summary = [column for column in columns if column not in leading and not column.startswith('condition_')]
    columns = [*leading, *summary, *[column for column in columns if column.startswith('condition_')]]
    for band in range(0, len(columns), columns_per_page):
        chosen = ['row', *columns[band:band + columns_per_page]]
        shown = frame[chosen]
        wrapped = [[textwrap.fill(_display(value), width=24, break_long_words=True)
                    for value in row] for row in shown.itertuples(index=False, name=None)]
        headers = [textwrap.fill(column.replace('_', ' ') if column[:1].islower() else column, width=22) for column in chosen]
        header_height = max(h.count('\n') + 1 for h in headers)
        line_budget = 52 - header_height
        start = 0
        while start < len(wrapped):
            end, lines = start, 0
            while end < len(wrapped):
                height = max(cell.count('\n') + 1 for cell in wrapped[end])
                if end > start and lines + height > line_budget:
                    break
                lines += height
                end += 1
            fig, ax = plt.subplots(figsize=(16.54, 11.69))
            try:
                fig.subplots_adjust(left=.035, right=.965, top=.88, bottom=.09)
                ax.axis('off')
                fig.text(.035, .955, title.replace('_', ' '), fontsize=15, weight='bold')
                fig.text(.035, .922, f'Rows {start + 1}-{end} / {len(frame)} | Columns {band + 1}-{min(band + columns_per_page, len(columns))} / {len(columns)}', fontsize=10)
                heights = [max(cell.count('\n') + 1 for cell in row) for row in wrapped[start:end]]
                total = sum(heights) + header_height
                table = ax.table(cellText=wrapped[start:end], colLabels=headers,
                                 cellLoc='left', bbox=[0, 1 - total / 54, 1, total / 54],
                                 colWidths=[.04] + [.96 / (len(chosen) - 1)] * (len(chosen) - 1))
                table.auto_set_font_size(False)
                table.set_fontsize(9)
                for (r, c), cell in table.get_celld().items():
                    cell.set_height((header_height if r == 0 else heights[r - 1]) / total)
                    cell.set_edgecolor('#d9e0e6')
                    cell.set_linewidth(.4)
                    if r == 0:
                        cell.set_facecolor('#243746')
                        cell.get_text().set_color('white')
                        cell.get_text().set_weight('bold')
                    else:
                        cell.set_facecolor('#f0f4f7' if r % 2 else 'white')
                        text = cell.get_text()
                        value = frame.iloc[start + r - 1]
                        if text.get_text() == 'N/A' or 'N/A' in text.get_text():
                            text.set_color('#777777')
                        elif (bool(value.get(f'{chosen[c]}_excludes_zero', False)) if
                              pd.notna(value.get(f'{chosen[c]}_excludes_zero', False)) else False) or frame.attrs.get('red_cells', {}).get((start + r - 1, chosen[c]), False):
                            text.set_color('#c62828')
                        elif chosen[c].endswith('q_value') and bool(value.get(f'{chosen[c]}_below_alpha', False)):
                            text.set_weight('bold')
                fig.text(.035, .035, f'Red: finite ordered contrast CI strictly excludes zero. Bold q: q < {alpha:g}. These are separate criteria. N/A: unavailable.\nDelta = condition 2 - condition 1. Intervals use the saved ci_level; no inference is recomputed. Full precision and source identities: companion CSV.', fontsize=9, color='#465565')
                add_table_page(fig, path, identifier=identifiers[0], source_ids=identifiers)
            finally:
                plt.close(fig)
            start = end


def _lipid_matrix(data):
    comp = data.loc[(data.family == 'profiles') & (data.section == 'comparison') & data.profile.fillna('').str.startswith('lipid_order__')].copy()
    if comp.empty:
        return pd.DataFrame()
    rows, red_cells = {}, {}
    for source in comp.to_dict('records'):
        species, tail, leaf = _lipid_identity(source['profile'])
        key = (species, source.get('position'), source.get('position_label'), source['system_1'], source['system_2'])
        row = rows.setdefault(key, dict(lipid=species, position=source.get('position'),
                                       position_label=source.get('position_label'),
                                       system_1=source['system_1'], system_2=source['system_2']))
        column = f'{leaf.capitalize()} {tail}'
        red_cells[(key, column)] = bool(source['ci_excludes_zero'])
        row[column] = _interval(source.get('difference_ci_low'), source.get('difference_ci_high'))
    matrix = pd.DataFrame(rows.values())
    columns = ['lipid', 'position', 'position_label', 'system_1', 'system_2']
    columns += [f'{leaf} {tail}' for tail in ('sn1', 'sn2', 'n_acyl', 'sphingosine')
                for leaf in ('Lower', 'Upper', 'Combined') if tail in {'sn1', 'sn2'} or
                any(any(column.endswith(f' {tail}') for column in row) for row in rows.values())]
    matrix = matrix.reindex(columns=columns)
    matrix.attrs['red_cells'] = {(i, column): flag for (key, column), flag in red_cells.items()
                                for i, current in enumerate(rows) if current == key}
    return matrix


def render_statistical_reports(cfg, root):
    """Publish one PDF/CSV per registered analysis from selected saved products."""
    data = collect_statistical_reports(cfg, root)
    if data.empty:
        return
    root = Path(root)
    alpha = cfg.get('statistics', {}).get('hypothesis_tests', {}).get('alpha', .05)
    for analysis, rows in data.groupby('analysis', sort=True):
        order = [column for column in ('section', 'family', 'units', 'profile', 'metric', 'group',
                 'interaction_type', 'interaction', 'system_1', 'system_2', 'system', 'replicate',
                 'position', 'position_label', 'source_table') if column in rows]
        rows = rows.sort_values(order, kind='stable', na_position='last').dropna(axis=1, how='all').reset_index(drop=True)
        identifiers = sorted(rows.source_id.unique())
        path = root / 'project_summary/figures/statistical_reports' / safe_plot_stem(analysis) / 'consolidated.png'
        matrix = _lipid_matrix(rows) if 'profile' in rows else pd.DataFrame()
        if not matrix.empty:
            for tails in (('sn1', 'sn2'), ('n_acyl', 'sphingosine')):
                cells = [f'{leaf} {tail}' for tail in tails for leaf in ('Lower', 'Upper', 'Combined')]
                cells = [column for column in cells if column in matrix]
                if not cells:
                    continue
                indices = matrix.index[matrix[cells].notna().any(axis=1)]
                compact = matrix.loc[indices, ['lipid', 'position_label', 'system_1', 'system_2', *cells]].copy()
                compact['contrast'] = compact.system_2 + ' - ' + compact.system_1
                compact = compact[['lipid', 'position_label', 'contrast', *cells]]
                compact.attrs['red_cells'] = {(i, column): matrix.attrs['red_cells'].get((original, column), False)
                                               for i, original in enumerate(indices) for column in cells}
                _pages(compact, f'{analysis} | pointwise lipid order contrast CI matrix (-S_CH)', path,
                       sorted(rows.loc[rows.section == 'comparison', 'source_id'].unique()), alpha=alpha,
                       columns_per_page=9)
        for (section, family, units), subset in rows.groupby(['section', 'family', 'units'], sort=True, dropna=False):
            ids = sorted(subset.source_id.unique())
            # Source metadata lives in CSV/manifest; every scientific column is
            # retained in the PDF, in column continuations with repeated row IDs.
            table = subset.drop(columns=['analysis', 'section', 'family', 'units', 'source_id', 'source_table', 'condition_source_table'], errors='ignore').dropna(axis=1, how='all')
            for low in list(table):
                high = low.removesuffix('_low') + '_high'
                if low.endswith('_ci_low') and high in table:
                    table[low.removesuffix('_low')] = [_interval(a, b) for a, b in zip(table[low], table[high])]
                    table = table.drop(columns=[low, high])
            title = f'{analysis} | {section} | {family} | units: {units}'
            diagnostics = [column for column in table if any(token in column for token in
                           ('tau_int', 'effective_n', 'block_', 'blocks_', 'sampling_', 'bootstrap_'))]
            if section != 'diagnostics' and diagnostics:
                identity = [column for column in table if column in {'profile', 'metric', 'group',
                            'interaction_type', 'interaction', 'position', 'position_label', 'system',
                            'system_1', 'system_2', 'replicate'}]
                _pages(table.drop(columns=diagnostics), title, path, ids, alpha=alpha)
                if cfg.get('plots', {}).get('bootstrap_tables', True):
                    _pages(table[[*identity, *diagnostics]], title + ' | diagnostics', path, ids, alpha=alpha)
            else:
                _pages(table, title, path, ids, alpha=alpha)
        sources = set(rows.source_table.unique())
        if 'condition_source_table' in rows:
            sources.update(rows.condition_source_table.dropna().unique())
        sources = sorted(sources)
        manifest = root / 'project_summary/statistics/manifest.json'
        if manifest.is_file():
            sources.append(manifest.relative_to(root).as_posix())
        add_report_data(path.with_name(REPORT_NAME), rows, {
            'analysis': analysis, 'source_ids': identifiers, 'selection': table_selection(cfg),
            'sources': [{'path': source, 'sha256': file_digest(root / source)} for source in sources],
            'alpha': alpha, 'row_count': len(rows), 'selection_policy': 'exact source IDs; legacy IDs select rows, never duplicate reports'})
