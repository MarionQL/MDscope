"""Syndecan lipid-region saved-data renderer readability regressions."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from mdscope.reporting.sdc4_condition_figures import (
    _region_figures, _region_label,
)

REGIONS = [
    "TM_A", "TM_B",
    "motif_NxIP_PROA_1", "motif_NxIP_PROB_1",
    "motif_acidic_patch_PROA_1", "motif_acidic_patch_PROB_1",
    "motif_proline_hydrophobic_PROA_1", "motif_proline_hydrophobic_PROB_1",
]
LIPIDS = ["CER160", "CHL1", "DOPC", "DPPC", "POPC", "POPE", "POPI", "POPS", "PSM"]


def _saved_reps(tmp_path, *, all_zero=False, many=True):
    reps = []
    for rep_id in range(1, 4):
        rep = tmp_path / "control" / f"rep{rep_id}"
        folder = rep / "extensions/sdc4_mechanisms"
        folder.mkdir(parents=True)
        rows = []
        for region in REGIONS:
            for species in (LIPIDS if many else LIPIDS[:2]):
                for leaflet in ("upper", "lower"):
                    value = 0. if all_zero else (
                        .2 + rep_id * .07 if species == "POPC" else
                        .05 + rep_id * .01
                    )
                    rows.append(dict(region=region, species=species,
                                     leaflet=leaflet, moiety="any",
                                     contact_fraction=value,
                                     local_composition_enrichment=value,
                                     matched_local_contact_fraction=value,
                                     matched_local_composition_enrichment=value))
        pd.DataFrame(rows).to_parquet(folder / "lipid_region_specificity.parquet")
        reps.append(rep)
    return reps


def _config():
    return {"plots": {"enabled": True,
                      "include": ["sdc4:control:lipid_region_mean_contact_fraction"]},
            "systems": {"control": {"label": "Control"}}}


def test_region_names_shorten_without_changing_data_identity():
    assert _region_label("TM_A") == "TM A"
    assert _region_label("motif_NxIP_PROB_1") == "NxIP B · 1"
    assert _region_label("motif_acidic_patch_PROA_1") == "Acidic A · 1"
    assert _region_label("motif_proline_hydrophobic_PROB_1") == "Proline B · 1"


def test_three_replicas_and_many_lipids_use_one_external_legend(tmp_path, monkeypatch):
    from mdscope.reporting import sdc4_condition_figures as renderer

    reps = _saved_reps(tmp_path)
    captured = []
    def inspect(fig, cfg, out, system, stem, **kwargs):
        assert fig.get_figwidth() >= 16
        assert fig.get_figheight() >= 10
        assert len(fig.axes) == 3
        assert len(fig.legends) == 1
        assert len(fig.legends[0].get_texts()) == len(LIPIDS) * 2
        assert all(ax.get_legend() is None for ax in fig.axes)
        assert all(ax.get_ylim() == pytest.approx(fig.axes[0].get_ylim())
                   for ax in fig.axes)
        assert all(ax.get_position().height > .15 for ax in fig.axes)
        labels = [tick.get_text() for tick in fig.axes[-1].get_xticklabels()]
        assert labels == [_region_label(name) for name in sorted(REGIONS)]
        captured.append(stem)
        plt.close(fig)
    monkeypatch.setattr(renderer, "_save", inspect)

    _region_figures(_config(), "control", reps,
                    tmp_path / "project_summary/figures/sdc4/control")
    assert captured == ["lipid_region_mean_contact_fraction"]
    assert all((rep / "extensions/sdc4_mechanisms/"
                "lipid_region_replicate_summary.csv").exists() for rep in reps)


def test_all_zero_region_plot_is_not_rendered_and_old_png_is_archived(tmp_path):
    reps = _saved_reps(tmp_path, all_zero=True)
    out = tmp_path / "project_summary/figures/sdc4/control"
    out.mkdir(parents=True)
    stale = out / "lipid_region_mean_contact_fraction.png"
    stale.write_bytes(b"previous flat-line png")
    _region_figures(_config(), "control", reps, out)
    assert not stale.exists()
    archive = tmp_path / ".archives/uninformative_sdc4_region_figures"
    moved = list(archive.glob("control_lipid_region_mean_contact_fraction*.png"))
    assert len(moved) == 1
    assert moved[0].read_bytes() == b"previous flat-line png"
    assert (reps[0] / "extensions/sdc4_mechanisms/"
            "lipid_region_replicate_summary.csv").is_file()
