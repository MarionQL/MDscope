"""Visual-review fixtures exercise the actual plotting functions, not fake drawings."""

import json
from pathlib import Path
from runpy import run_path

from PIL import Image
import pytest


def _scripts():
    root = Path(__file__).resolve().parents[1] / "scripts"
    preview = run_path(str(root / "plot_preview.py"))
    gallery = run_path(str(root / "figure_gallery.py"))
    return preview, gallery


def test_synthetic_plot_preview_is_read_only_with_respect_to_trajectory(tmp_path, monkeypatch):
    import MDAnalysis as mda

    def forbidden(*args, **kwargs):
        raise AssertionError("Preview must not load topology or trajectory")
    monkeypatch.setattr(mda, "Universe", forbidden)

    preview, _ = _scripts()
    root = tmp_path / "preview"
    report = preview["render_preview"](root)
    assert "SYNTHETIC" in report["notice"]
    assert report["plots"]
    expected = {
        "rmsd_A__timeseries.png",
        "lipid_order__POPC__sn1__profile.png",
        "lipid_order__POPC__sn1__system_difference.png",
        "occupancy_heatmap.png", "lifetime_distribution.png",
        "shell_rdf.png", "contact_survival.png", "water_and_ions.png",
        "lipid_enrichment_heatmap.png", "diffusion_coefficients.png",
        "state_populations.png", "nxip_height_occlusion_joint.png",
    }
    assert expected <= {Path(name).name for name in report["plots"]}
    assert len(report["plots"]) >= 30
    for name in report["plots"]:
        with Image.open(root / name) as im:
            im.verify()
        with Image.open(root / name) as im:
            assert im.width >= 300 and im.height >= 150
            # A separate header band must not overwrite the Matplotlib figure.
            assert im.getpixel((0, 0)) == (255, 255, 255)
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        preview["render_preview"](root)


def test_preview_plot_axes_and_legend_use_scientific_units(tmp_path, monkeypatch):
    from mdscope.reporting import plotting

    preview, _ = _scripts()
    original = plotting._save
    observed = []

    def check_before_save(fig, path, cfg):
        if Path(path).name == "rmsd_A__timeseries.png":
            observed.append(True)
            axes = [ax for ax in fig.axes if ax.has_data()]
            assert axes and all(ax.get_xlabel() == "Time (ns)" for ax in axes)
            assert all(ax.get_ylabel() == "A" for ax in axes)
            assert all(ax.get_legend() is not None for ax in axes)
        original(fig, path, cfg)

    monkeypatch.setattr(plotting, "_save", check_before_save)
    preview["render_preview"](tmp_path / "figure")
    assert observed, "No RMSD time-series renderer was exercised"


def test_html_gallery_displays_before_after_and_no_scientific_claims(tmp_path):
    preview, gallery = _scripts()
    base, head = tmp_path / "base", tmp_path / "head"
    preview["render_preview"](base)
    preview["render_preview"](head)
    output = tmp_path / "index.html"
    n = gallery["build_gallery"](base, head, output)
    text = output.read_text()
    assert n >= 4
    assert "SYNTHETIC EXAMPLES" in text
    assert "Previous code" in text and "Updated code" in text
    assert "rmsd_A__timeseries.png" in text
    assert "Unchanged pixels" in text
    assert json.loads((base / "preview_manifest.json").read_text())["plots"] == (
        json.loads((head / "preview_manifest.json").read_text())["plots"]
    )


def test_csv_samples_include_every_column_but_only_bounded_rows(tmp_path):
    import csv
    from runpy import run_path

    preview, _ = _scripts()
    root = tmp_path / "preview"
    report = preview["render_preview"](root)
    assert report["table_previews"]
    catalog = root / report["table_inventory_csv"]
    with catalog.open(newline="") as handle:
        entries = list(csv.DictReader(handle))
    assert len(entries) >= 15
    named = next(item for item in entries
                 if item["source_parquet"] == "condition_A/rep1/summary/timeseries.parquet")
    assert int(named["n_rows"]) == 61
    assert int(named["sample_rows"]) == 5
    with (root / named["columns_csv"]).open(newline="") as handle:
        cols = list(csv.DictReader(handle))
    assert "time_ps" in {row["column"] for row in cols}
    assert "rmsd_A" in {row["column"] for row in cols}
    with (root / named["sample_csv"]).open(newline="") as handle:
        sample = list(csv.DictReader(handle))
    assert len(sample) == 5
    assert set(sample[0]) == {row["column"] for row in cols}
    assert all("SYNTHETIC" not in row["time_ps"] for row in sample)
    assert "SYNTHETIC" in (root / "table_previews/README.txt").read_text()


def test_preview_gallery_groups_condition_figures_without_confusing_revisions(tmp_path):
    from runpy import run_path

    preview, _ = _scripts()
    root = tmp_path / "current"
    preview["render_preview"](root)
    builder = run_path(str(Path(__file__).resolve().parents[1]
                           / "scripts/preview_gallery.py"))["build_preview_gallery"]
    main = tmp_path / "index.html"
    n = builder(root, main)
    html = main.read_text()
    assert n >= 30
    assert "condition-centered figure preview" in html
    assert 'id="condition_A"' in html
    assert 'id="condition_B"' in html
    assert 'id="project"' in html
    assert "CSV table previews" in html
    assert "table_inventory.csv" in html
    assert "before/after comparison" in html
    # Both replicas are clearly tagged within their condition.
    assert "rep1" in html and "rep2" in html
    # No automatic cross-labeling of one replica as a different condition.
    a = html.split('<section id="condition_A">')[1].split("</section>")[0]
    assert "condition_B/rep2/extensions" not in a

