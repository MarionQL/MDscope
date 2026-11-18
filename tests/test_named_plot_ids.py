"""Individual-figure allowlists must not change metric semantics or load excluded families."""
from copy import deepcopy
from pathlib import Path
from runpy import run_path

from PIL import Image
import pytest
import yaml

from mdscope.config import DEFAULTS, figure_enabled, load_config, plot_enabled
from mdscope.reporting.plot_ids import id_from_figure_path, validate_plot_ids
from mdscope.reporting.plotting import plot_all, plot_interactions
from mdscope.reporting.figure_inventory import figure_status


def _fake_data(tmp_path):
    root = Path(__file__).resolve().parents[1]
    functions = run_path(str(root / "scripts" / "plot_preview.py"))
    return functions["create_fixture"](tmp_path)


def test_exact_figure_allowlist_preserves_legacy_defaults():
    cfg = deepcopy(DEFAULTS)
    cfg["plots"]["enabled"] = True
    assert plot_enabled(cfg, "time_series")
    assert figure_enabled(cfg, "time_series:rmsd_A:timeseries")
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    assert plot_enabled(cfg, "time_series")
    assert not plot_enabled(cfg, "profiles")
    assert figure_enabled(cfg, "time_series:rmsd_A:timeseries")
    assert not figure_enabled(cfg, "time_series:rmsd_A:distribution")
    cfg["plots"]["families"] = {"time_series": False}
    assert not figure_enabled(cfg, "time_series:rmsd_A:timeseries")


@pytest.mark.parametrize("item", [
    "time_series:rmsd_A:bogus", "profiles:foo:bad",
    "interactions:kind:unregistered", "sdc4:anything",
    "solvation:no_such_plot", "diffusion:nope",
    "unknown:foo:timeseries", "time_series::timeseries",
])
def test_unsupported_plot_ids_are_rejected(item):
    with pytest.raises(ValueError, match="unsupported plot ID"):
        validate_plot_ids([item])


def test_named_plot_allowlist_is_validated_before_trajectory_work(tmp_path, project):
    path, root = project("selected_fig")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"] = {"enabled": True, "include": [
        "time_series:rmsd_A:timeseries", "solvation:shell_rdf",
        "profiles:residue_membrane_distance__membrane:tiny:distance_heatmap",
    ]}
    path.write_text(yaml.safe_dump(cfg))
    result = load_config(path)
    assert result["plots"]["include"] == cfg["plots"]["include"]
    assert not root.exists()


def test_only_selected_time_series_is_saved(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, tmp_path, repouts)
    produced = sorted((tmp_path / "project_summary/figures").rglob("*.png"))
    assert [path.name for path in produced] == ["rmsd_A__timeseries.png"]
    assert id_from_figure_path(produced[0]) == "time_series:rmsd_A:timeseries"
    with Image.open(produced[0]) as figure:
        figure.verify()


def test_selected_distribution_does_not_write_time_series(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:distribution"]
    cfg["plots"]["distributions"] = True
    plot_all(cfg, tmp_path, repouts)
    produced = sorted((tmp_path / "project_summary/figures").rglob("*.png"))
    assert [path.name for path in produced] == ["rmsd_A__distribution.png"]


def test_interaction_filter_is_pushed_before_loading_tables(tmp_path, monkeypatch):
    import mdscope.reporting.plotting as module

    seen = []
    monkeypatch.setattr(module, "interaction_tables", lambda *args, **kwargs:
                        (seen.append(kwargs["include_kinds"]) or {}))
    cfg = {"plots": {"enabled": True, "include": [
        "interactions:hydrogen_bonds:occupancy_heatmap",
    ]}, "statistics": {"enabled": False}, "systems": {
        "condition": {"replicates": [{"label": "rep1"}]},
    }}
    plot_interactions(cfg, tmp_path, {"condition": [tmp_path / "rep1"]})
    assert seen == [{"hydrogen_bonds"}]


def test_old_generated_figure_excluded_when_no_longer_selected(tmp_path):
    rep = tmp_path / "condition" / "rep1"
    rep.mkdir(parents=True)
    dest = tmp_path / "project_summary/figures/time_series/rmsd_A__timeseries.png"
    dest.parent.mkdir(parents=True)
    Image.new("RGB", (80, 50), "white").save(dest)
    cfg = {"plots": {"enabled": True, "include": [
        "time_series:radius_of_gyration_A:timeseries",
    ]}, "statistics": {"enabled": False}}
    assert figure_status(dest, cfg, tmp_path, {"condition": [rep]}) == "disabled"


def test_solvation_not_selected_never_opens_rdf(tmp_path, monkeypatch):
    import mdscope.reporting.solvation_figures as module

    def forbidden(*args, **kwargs):
        raise AssertionError("solvation RDF should not be loaded")
    monkeypatch.setattr(module, "solvation_rdf_figure", forbidden)
    monkeypatch.setattr(module, "solvation_timeseries_figure", lambda *args: None)
    cfg = {"plots": {"enabled": True, "include": ["solvation:water_and_ions"]}}
    module.plot_solvation(cfg, tmp_path, {"condition": []})


def test_sdc4_replica_selection_does_not_open_unselected_maps(tmp_path):
    import numpy as np
    import pandas as pd
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms

    rep = tmp_path / "condition" / "rep1"
    folder = rep / "extensions/sdc4_mechanisms"
    folder.mkdir(parents=True)
    # A corrupt non-selected map must never be opened.
    (folder / "lipid_density_maps.npz").write_bytes(b"not npz")
    (folder / "lipid_region_specificity.parquet").write_bytes(b"not parquet")
    t = np.linspace(0, 1, 20)
    pd.DataFrame({
        "time_ps": t * 1000,
        "sdc4_ectodomain_A_height_A": 17 + 2*t,
        "sdc4_ectodomain_B_height_A": 18 + t,
        "sdc4_motif_NxIP_PROA_1_steric_occlusion_fraction": t,
        "sdc4_motif_NxIP_PROB_1_steric_occlusion_fraction": 1-t,
    }).to_parquet(folder / "timeseries.parquet")
    cfg = {
        "plots": {"enabled": True,
                  "include": ["sdc4:condition:rep1:nxip_height_occlusion_joint"]},
        "analyses": {"sdc4_mechanisms": {"enabled": True}},
    }
    plot_sdc4_mechanisms(cfg, tmp_path, {"condition": [rep]})
    assert (folder / "nxip_height_occlusion_joint.png").is_file()
    assert not (folder / "vvhplvpl_proline_ramachandran.png").exists()
    assert id_from_figure_path(folder / "nxip_height_occlusion_joint.png") == (
        "sdc4:condition:rep1:nxip_height_occlusion_joint")


def test_shared_state_sampling_only_does_not_require_populations(tmp_path):
    import pandas as pd
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms

    out = tmp_path / "project_summary/shared_states"
    out.mkdir(parents=True)
    pd.DataFrame({
        "system": ["condition"], "replicate": ["rep1"],
        "first_second_half_js_divergence": [0.015],
    }).to_parquet(out / "sampling_diagnostics.parquet")
    cfg = {
        "plots": {"enabled": True,
                  "include": ["shared_states:state_sampling_diagnostics"]},
        "analyses": {"sdc4_mechanisms": {"enabled": True}},
    }
    plot_sdc4_mechanisms(cfg, tmp_path, {"condition": []})
    assert (out / "state_sampling_diagnostics.png").is_file()
    assert not (out / "state_populations.png").exists()


def test_lipid_contact_include_one_exact_heatmap(tmp_path):
    import pandas as pd
    from mdscope.reporting.membrane_figures import plot_lipid_contacts

    rep = tmp_path / "condition" / "rep1"
    matrix = rep / "extensions/lipid_contacts/matrix.parquet"
    matrix.parent.mkdir(parents=True)
    pd.DataFrame({
        "time_ps": [0.0, 1000.0], "leaflet": ["upper"]*2,
        "species_a": ["POPC"]*2, "species_b": ["DOPC"]*2,
        "neighbors_per_molecule": [1.3, 1.6], "enrichment": [1.15, 1.22],
        "pair_count_undirected": [10, 12],
        "molecules_a": [20, 20], "molecules_b": [20, 20],
    }).to_parquet(matrix)
    cfg = {
        "plots": {"enabled": True,
                  "include": ["lipid_contacts:lipid_contact_matrix_condition_upper_enrichment"]},
        "analyses": {"lipid_contacts": {"enabled": True}},
    }
    plot_lipid_contacts(cfg, tmp_path, {"condition": [rep]})
    dest = tmp_path / "project_summary/figures/lipid_contacts"
    names = sorted(p.name for p in dest.iterdir())
    assert names == ["lipid_contact_matrix_condition_upper_enrichment.png"]
    assert id_from_figure_path(dest / names[0]) == cfg["plots"]["include"][0]


def test_specialized_figure_ids_are_validated():
    for item in [
        "sdc4:condition:rep1:nxip_height_occlusion_joint",
        "sdc4:condition:rep1:lipid_density_maps",
        "shared_states:state_sampling_diagnostics",
        "lipid_contacts:lipid_contact_matrix_condition_upper_enrichment",
        "lipid_contacts:lipid_contact_difference_first__second_upper_neighbors_per_molecule",
    ]:
        validate_plot_ids([item])
    for item in [
        "sdc4:condition:rep1:unknown",
        "lipid_contacts:other.png", "shared_states:some_chart",
    ]:
        with pytest.raises(ValueError, match="unsupported plot ID"):
            validate_plot_ids([item])


def test_single_metric_ignores_unselected_solvation_and_lipid_tables(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    for paths in repouts.values():
        for rep in paths:
            for folder, name in (
                ("lipid_contacts", "matrix.parquet"),
                ("solvation", "species_timeseries.parquet"),
                ("solvation", "bridge_pair_timeseries.parquet"),
                ("hydration", "aggregate_counts.parquet"),
            ):
                path = rep / "extensions" / folder / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"invalid high-cardinality parquet that must not be read")
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, tmp_path, repouts)
    assert sorted(path.name for path in (
        tmp_path / "project_summary/figures").rglob("*.png")) == [
        "rmsd_A__timeseries.png"
    ]


def test_selected_lipid_contact_time_series_still_reads_needed_matrix(tmp_path):
    import pandas as pd

    cfg, repouts = _fake_data(tmp_path)
    rep = repouts["condition_A"][0]
    file = rep / "extensions/lipid_contacts/matrix.parquet"
    file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "time_ps": [0., 1000.],
        "leaflet": ["upper", "upper"],
        "species_a": ["POPC", "POPC"], "species_b": ["DOPC", "DOPC"],
        "neighbors_per_molecule": [1.2, 1.3], "enrichment": [1.1, 1.2],
    }).to_parquet(file)
    cfg["plots"]["include"] = [
        "time_series:enrichment__group_lipid_contacts_upper_POPC_DOPC:timeseries"
    ]
    plot_all(cfg, tmp_path, repouts)
    produced = [path.name for path in (
        tmp_path / "project_summary/figures").rglob("*.png")]
    assert produced == ["enrichment__group_lipid_contacts_upper_POPC_DOPC__timeseries.png"]


def test_individual_dpi_changes_only_selected_visual_resolution(tmp_path):
    root_base, root_hires = tmp_path / "base", tmp_path / "hires"
    cfg, repouts = _fake_data(root_base)
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, root_base, repouts)
    image_path = "project_summary/figures/time_series/rmsd_A__timeseries.png"
    with Image.open(root_base / image_path) as image:
        base_pixels = image.size

    hires, reps_hi = _fake_data(root_hires)
    hires["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    hires["plots"]["options"] = {"time_series:rmsd_A:timeseries": {"dpi": 250}}
    plot_all(hires, root_hires, reps_hi)
    with Image.open(root_hires / image_path) as image:
        hires_pixels = image.size
    assert hires_pixels[0] > base_pixels[0] * 1.5
    assert hires_pixels[1] > base_pixels[1] * 1.5


@pytest.mark.parametrize("options,diagnostic", [
    ({"time_series:rmsd_A:timeseries": {"dpi": 30}}, "dpi"),
    ({"time_series:rmsd_A:timeseries": {"dpi": True}}, "dpi"),
    ({"time_series:rmsd_A:timeseries": {"legend": True}}, "supports only dpi"),
    ({"unknown:plot": {"dpi": 150}}, "unsupported plot ID"),
])
def test_per_plot_options_reject_bad_keys_and_values(options, diagnostic, tmp_path, project):
    path, _ = project("per_plot_bad")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"]["options"] = options
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match=diagnostic):
        load_config(path)

