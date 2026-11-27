"""Condition-level SDC4 figures and shared conformational-state diagnostics."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import figure_enabled, plot_enabled
from .plot_ids import figure_dpi, figure_path, save_figure
from .sdc4_condition_figures import plot_sdc4_condition_figures


def _plot_shared_state_figures(cfg, root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    state_dir = Path(root) / "project_summary/shared_states"
    if figure_enabled(cfg, "shared_states:state_populations"):
        pop_path = state_dir / "replicate_populations.parquet"
        if pop_path.is_file():
            pop = pd.read_parquet(pop_path)
            if not pop.empty:
                systems = list(dict.fromkeys(pop["system"]))
                states = sorted(pop["state"].unique())
                fig, ax = plt.subplots(figsize=(max(7, len(states)*2.2), 4.5))
                width = .8 / max(len(systems), 1)
                for j, system in enumerate(systems):
                    sub = pop[pop.system == system]
                    means = [sub.loc[sub.state == state, "fraction"].mean() for state in states]
                    xpos = np.arange(len(states)) + (j - (len(systems)-1)/2)*width
                    ax.bar(xpos, means, width=width*.9, alpha=.6, label=system)
                    for i, state in enumerate(states):
                        values = sub.loc[sub.state == state, "fraction"].to_numpy(float)
                        offsets = np.linspace(-width*.22, width*.22, len(values))
                        ax.scatter(np.full(len(values), xpos[i]) + offsets, values,
                                   color="black", s=13, zorder=3)
                ax.set_xticks(np.arange(len(states)), [f"State {s}" for s in states])
                ax.set_ylabel("Fraction of analyzed frames")
                ax.set_ylim(0, 1)
                ax.set_title("Shared-state populations\n(points = independent replicates)")
                ax.legend()
                fig.tight_layout()
                save_figure(fig, cfg, figure_path(cfg, state_dir / "state_populations.png"),
                            dpi=figure_dpi(cfg, "shared_states:state_populations", 200))
                plt.close(fig)

    if figure_enabled(cfg, "shared_states:state_sampling_diagnostics"):
        diag_path = state_dir / "sampling_diagnostics.parquet"
        if diag_path.is_file():
            diag = pd.read_parquet(diag_path)
            if not diag.empty:
                fig, ax = plt.subplots(figsize=(max(6, len(diag)*1.1), 4))
                labels = [f"{row.system}/{row.replicate}" for row in diag.itertuples()]
                ax.bar(np.arange(len(diag)), diag.first_second_half_js_divergence)
                ax.set_xticks(np.arange(len(diag)), labels, rotation=45, ha="right")
                ax.set_ylabel("Jensen–Shannon divergence (bits)")
                ax.set_title("First- vs second-half state populations")
                fig.tight_layout()
                save_figure(fig, cfg, figure_path(cfg, state_dir / "state_sampling_diagnostics.png"),
                            dpi=figure_dpi(cfg, "shared_states:state_sampling_diagnostics", 200))
                plt.close(fig)


def plot_sdc4_mechanisms(cfg, root, system_repouts):
    if not cfg.get('plots',{}).get('enabled',True):
        return
    if plot_enabled(cfg, 'shared_states'):
        _plot_shared_state_figures(cfg,root)
    if (cfg.get('analyses',{}).get('sdc4_mechanisms',{}).get('enabled',False)
            and plot_enabled(cfg, 'sdc4')):
        plot_sdc4_condition_figures(cfg,root,system_repouts)
