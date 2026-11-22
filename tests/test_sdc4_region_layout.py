"""Syndecan lipid-region saved-data renderer readability regressions."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    assert all((tmp_path / "project_summary/reporting_exports/sdc4/control" /
                rep.name / "lipid_region_replicate_summary.csv").exists() for rep in reps)
    assert all(not (rep / "extensions/sdc4_mechanisms/"
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
    assert (tmp_path / "project_summary/reporting_exports/sdc4/control" /
            reps[0].name / "lipid_region_replicate_summary.csv").is_file()


def test_plot_only_rerender_never_rewrites_legacy_replica_csvs(tmp_path):
    reps = _saved_reps(tmp_path, many=False)
    originals = {}
    for rep in reps:
        csv = rep / "extensions/sdc4_mechanisms/lipid_region_replicate_summary.csv"
        csv.write_bytes(b"previous-scientific-output")
        originals[csv] = (csv.read_bytes(), csv.stat().st_mtime_ns)
    numerical = {path: (path.read_bytes(), path.stat().st_mtime_ns)
                 for rep in reps for path in rep.rglob('*') if path.is_file()}
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms
    cfg = _config()
    cfg['analyses'] = {'sdc4_mechanisms': {'enabled': True}}
    for _ in range(2):
        plot_sdc4_mechanisms(cfg, tmp_path, {'control': reps})
    assert {path: (path.read_bytes(), path.stat().st_mtime_ns)
            for rep in reps for path in rep.rglob('*') if path.is_file()} == numerical
    for csv, state in originals.items():
        assert (csv.read_bytes(), csv.stat().st_mtime_ns) == state
    export = tmp_path / "project_summary/reporting_exports/sdc4/control/rep1/lipid_region_replicate_summary.csv"
    frame = pd.read_csv(export)
    assert frame.columns.tolist() == [
        "region", "species", "leaflet", "mean_contact_fraction",
        "mean_local_enrichment", "mean_matched_local_contact_fraction",
        "mean_matched_local_enrichment", "n_matched_local_frames",
        "n_matched_enrichment_frames", "n_frames_with_species",
    ]
    value = frame.loc[(frame.region == "TM_A") &
                      (frame.species == "CER160") &
                      (frame.leaflet == "upper"), "mean_contact_fraction"]
    assert value.iloc[0] == pytest.approx(.06)


def test_missing_or_empty_saved_regions_do_not_write_replica_csvs(tmp_path):
    from mdscope.reporting.sdc4_exports import export_region_replica_summaries
    rep = tmp_path / "control/rep1"
    rep.mkdir(parents=True)
    export_region_replica_summaries(tmp_path, "control", [rep], {})
    export_region_replica_summaries(tmp_path, "control", [rep], {
        rep.name: pd.DataFrame({"moiety": ["glycan"], "region": ["TM_A"]})
    })
    assert not list((tmp_path / "project_summary").rglob("*.csv"))


def test_reporting_summary_preserves_missing_matched_fields(tmp_path):
    from mdscope.reporting.sdc4_exports import export_region_replica_summaries
    rep = tmp_path / "control/rep1"
    rep.mkdir(parents=True)
    tab = pd.DataFrame({
        "moiety": ["any", "any"], "region": ["TM_A", "TM_A"],
        "species": ["POPC", "POPC"], "leaflet": ["upper", "upper"],
        "contact_fraction": [.2, .4],
        "local_composition_enrichment": [1.0, 2.0],
    })
    export_region_replica_summaries(tmp_path, "control", [rep], {rep.name: tab})
    export = tmp_path / "project_summary/reporting_exports/sdc4/control/rep1/lipid_region_replicate_summary.csv"
    row = pd.read_csv(export).iloc[0]
    assert row.mean_contact_fraction == pytest.approx(.3)
    assert row.mean_local_enrichment == pytest.approx(1.5)
    assert row.n_frames_with_species == 2
    assert row.n_matched_local_frames == 0
    assert pd.isna(row.mean_matched_local_contact_fraction)


@pytest.mark.parametrize('missing', [True, False])
def test_empty_reporting_export_archives_previous_presentation_csv(tmp_path, missing):
    from mdscope.reporting.sdc4_exports import export_region_replica_summaries
    rep=tmp_path/'control/rep1'
    output=tmp_path/'project_summary/reporting_exports/sdc4/control/rep1/lipid_region_replicate_summary.csv'
    output.parent.mkdir(parents=True); output.write_bytes(b'previous presentation summary')
    tables={} if missing else {rep.name: pd.DataFrame({'moiety':['glycan']})}
    export_region_replica_summaries(tmp_path,'control',[rep],tables)
    assert not output.exists()
    archived=list((tmp_path/'.archives').rglob('*.csv'))
    assert len(archived)==1 and archived[0].read_bytes()==b'previous presentation summary'


def test_all_missing_replicas_archive_stale_export_through_production_renderer(tmp_path):
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms
    rep=tmp_path/'control/rep1'
    export=tmp_path/'project_summary/reporting_exports/sdc4/control/rep1/lipid_region_replicate_summary.csv'
    export.parent.mkdir(parents=True); export.write_bytes(b'old presentation summary')
    cfg=_config(); cfg['analyses']={'sdc4_mechanisms':{'enabled':True}}
    plot_sdc4_mechanisms(cfg,tmp_path,{'control':[rep]})
    assert not export.exists()
    assert next((tmp_path/'.archives').rglob('*.csv')).read_bytes()==b'old presentation summary'
