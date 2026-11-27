"""Visual-review fixtures exercise the actual plotting functions, not fake drawings."""

import json
from pathlib import Path
from runpy import run_path

import pytest
from PIL import Image


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
    # Stable IDs survive the documented directory/filename migration.
    expected_ids = {'time_series:rmsd_A:timeseries', 'profiles:lipid_order__POPC__sn1:profile', 'profiles:lipid_order__POPC__sn1:system_difference', 'interactions:salt_bridges:occupancy_heatmap', 'solvation:shell_rdf', 'solvation:contact_survival', 'shared_states:state_populations', 'diffusion:diffusion_coefficients', 'sdc4:condition_A:nxip_height_occlusion_joint'}
    assert expected_ids <= set(report['plot_ids'].values())
    assert len(report['plot_ids']) == len(report['plots'])
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
    assert "time_series@rmsd_A@timeseries.png" in text
    assert "Unchanged pixels" in text
    assert json.loads((base / "preview_manifest.json").read_text())["plots"] == (
        json.loads((head / "preview_manifest.json").read_text())["plots"]
    )
    assert "Compared:" not in text
    gallery["build_gallery"](base, head, output, revisions="base <abc1234> vs head def5678")
    assert "Compared: base &lt;abc1234&gt; vs head def5678." in output.read_text()


def test_csv_samples_include_every_column_but_only_bounded_rows(tmp_path):
    import csv

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


def test_single_plot_developer_review_summarizes_extras_without_repeating_them(tmp_path):
    preview, gallery = _scripts()
    base, head = tmp_path / "before", tmp_path / "after"
    before = preview["render_preview"](base)
    after = preview["render_preview"](head, include=["time_series:rmsd_A:timeseries"])
    output = tmp_path / "review/selection/index.html"
    count = gallery["build_gallery"](
        base, head, output, only_name="time_series@rmsd_A@timeseries.png")
    html = output.read_text()
    assert count == 1
    assert f"Previous code generated {len(before['plots']) - 1} other figure(s)" in html
    assert f"updated code generated {len(after['plots']) - 1} other figure(s)" in html
    assert html.count("<section>") == 1
    assert "time_series@rmsd_A@timeseries.png" in html
    assert "<h2>project_summary/figures/solvation/" not in html



def test_ci_artifact_cleanup_keeps_selected_ids_after_path_migration(tmp_path, monkeypatch):
    import yaml
    workflow = yaml.safe_load((Path(__file__).parents[1] / '.github/workflows/figure-review.yml').read_text())
    run = next(step['run'] for step in workflow['jobs']['before-after-plots']['steps']
               if step.get('name') == 'Render identical synthetic inputs with base and head renderers')
    cleanup = run.split("python - <<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
    root = tmp_path/'mdscope-plot-review'
    old = 'project_summary/figures/time_series/rmsd_A__timeseries.png'
    new = 'project_summary/figures/rmsd/secondary/time_series@rmsd_A@timeseries.png'
    for relative, selected in [('review/previous', old), ('review/previous-selected', old), ('review/updated-selected', new)]:
        folder = root/relative
        extra = 'extra.png'
        for name in (selected, extra):
            path = folder/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'png fixture')
        (folder/'preview_manifest.json').write_text(json.dumps({'plots': [selected, extra],
            'plot_ids': {selected: 'time_series:rmsd_A:timeseries', extra: 'other'}}))
    monkeypatch.setenv('RUNNER_TEMP', str(tmp_path))
    monkeypatch.setenv('REVISIONS', 'synthetic test')
    exec(compile(cleanup, 'figure-review-cleanup', 'exec'), {})
    folder = root/'review/updated-selected'
    record = json.loads((folder/'preview_manifest.json').read_text())
    assert record['plots'] == [new]
    assert record['unselected_generated_plots'] == 1
    assert (folder/new).is_file()
    assert not (folder/'extra.png').exists()
