"""Reusable, unit-explicit paginated heatmaps; all values are supplied by callers."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def heatmap_pages(values, row_labels, column_labels, *, title, units,
                  vmin, vmax, cmap='viridis', rows_per_page=36, columns_per_page=24):
    """Yield every rectangular tile with shared color limits and exact slices.

    Tiles preserve input identity order, including zero and missing values.
    Callers can annotate using the returned slices without estimating anything.
    """
    values = np.asarray(values, dtype=float)
    if values.shape != (len(row_labels), len(column_labels)):
        raise ValueError('Heatmap values and labels must have identical shape')
    colors = plt.get_cmap(cmap).copy()
    colors.set_bad('#e0e0e0')
    page = 0
    for row in range(0, len(row_labels), rows_per_page):
        for column in range(0, len(column_labels), columns_per_page):
            page += 1
            rs = slice(row, row + rows_per_page)
            cs = slice(column, column + columns_per_page)
            labels_y, labels_x = list(row_labels[rs]), list(column_labels[cs])
            longest = max(map(lambda label: len(str(label)), labels_y), default=0)
            fig, ax = plt.subplots(figsize=(max(8, .45*len(labels_x)+3) + min(6, .07*longest),
                                            max(4, .30*len(labels_y)+2)))
            image = ax.imshow(np.ma.masked_invalid(values[rs, cs]), aspect='auto',
                              interpolation='nearest', vmin=vmin, vmax=vmax, cmap=colors)
            ax.set_xticks(range(len(labels_x)), labels_x, rotation=55, ha='right', fontsize=8)
            ax.set_yticks(range(len(labels_y)), labels_y, fontsize=8)
            ax.set_title(title, fontsize=11)
            fig.colorbar(image, ax=ax, label=units)
            yield fig, ax, rs, cs, page
