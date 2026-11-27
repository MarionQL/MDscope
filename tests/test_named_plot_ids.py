"""Individual-figure allowlists must not change metric semantics or load excluded families."""
from copy import deepcopy
from pathlib import Path
from runpy import run_path

import pytest
import yaml
from PIL import Image

from mdscope.config import DEFAULTS, figure_enabled, load_config, plot_enabled
from mdscope.contracts.plot_validation import validate_plot_ids
from mdscope.reporting.figure_inventory import figure_status
from mdscope.reporting.plot_ids import id_from_figure_path
from mdscope.reporting.plotting import plot_all, plot_interactions


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
    assert [path.name for path in produced] == ["time_series@rmsd_A@timeseries.png"]
    assert id_from_figure_path(produced[0]) == "time_series:rmsd_A:timeseries"
    with Image.open(produced[0]) as figure:
        figure.verify()


def test_selected_distribution_does_not_write_time_series(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:distribution"]
    cfg["plots"]["distributions"] = True
    plot_all(cfg, tmp_path, repouts)
    produced = sorted((tmp_path / "project_summary/figures").rglob("*.png"))
    assert [path.name for path in produced] == ["time_series@rmsd_A@distribution.png"]


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
    dest = tmp_path / "project_summary/figures/rmsd/secondary/time_series@rmsd_A@timeseries.png"
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
    monkeypatch.setattr(module, "solvation_timeseries_figure", lambda *args, **kwargs: None)
    cfg = {"plots": {"enabled": True, "include": ["solvation:water_and_ions"]}}
    module.plot_solvation(cfg, tmp_path, {"condition": []})


def test_sdc4_condition_selection_does_not_open_unselected_maps(tmp_path):
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
                  "include": ["sdc4:condition:nxip_height_occlusion_joint"]},
        "analyses": {"sdc4_mechanisms": {"enabled": True}},
    }
    plot_sdc4_mechanisms(cfg, tmp_path, {"condition": [rep]})
    output = tmp_path / "project_summary/figures/sdc4_mechanisms/secondary"
    assert (output / "sdc4@condition@nxip_height_occlusion_joint.png").is_file()
    assert not (output / "vvhplvpl_proline_ramachandran.png").exists()
    assert not list(folder.glob("*.png"))
    assert id_from_figure_path(output / "sdc4@condition@nxip_height_occlusion_joint.png") == (
        "sdc4:condition:nxip_height_occlusion_joint")


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
    assert (tmp_path / "project_summary/figures/shared_states/tertiary/shared_states@state_sampling_diagnostics.png").is_file()
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
    names = sorted(p.name for p in dest.rglob("*.png"))
    assert names == ["lipid_contacts@lipid_contact_matrix_condition_upper_enrichment.png"]
    assert id_from_figure_path(next(dest.rglob("*.png"))) == cfg["plots"]["include"][0]


def test_specialized_figure_ids_are_validated():
    for item in [
        "sdc4:condition:nxip_height_occlusion_joint",
        "sdc4:condition:lipid_density_maps",
        "shared_states:state_sampling_diagnostics",
        "lipid_contacts:lipid_contact_matrix_condition_upper_enrichment",
        "lipid_contacts:lipid_contact_difference_first__second_upper_neighbors_per_molecule",
    ]:
        validate_plot_ids([item])
    for item in [
        "sdc4:condition:unknown",
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
        "time_series@rmsd_A@timeseries.png"
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
    assert produced == ["time_series@enrichment__group_lipid_contacts_upper_POPC_DOPC@timeseries.png"]


def test_individual_dpi_changes_only_selected_visual_resolution(tmp_path):
    root_base, root_hires = tmp_path / "base", tmp_path / "hires"
    cfg, repouts = _fake_data(root_base)
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, root_base, repouts)
    image_path = "project_summary/figures/rmsd/secondary/time_series@rmsd_A@timeseries.png"
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
    ({"time_series:rmsd_A:timeseries": {"legend": True}}, "legend"),
    ({"unknown:plot": {"dpi": 150}}, "unsupported plot ID"),
])
def test_per_plot_options_reject_bad_keys_and_values(options, diagnostic, tmp_path, project):
    path, _ = project("per_plot_bad")
    cfg = yaml.safe_load(path.read_text())
    cfg["plots"]["options"] = options
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match=diagnostic):
        load_config(path)


def test_describe_plot_uses_declared_renderers_and_saved_input_contracts():
    from mdscope.contracts.plot_validation import describe_plot

    for identifier, expected in (
        ("time_series:rmsd_A:timeseries", "summary/timeseries.parquet"),
        ("profiles:lipid_order__POPC__sn1:profile", "replicate_profiles.parquet"),
        ("interactions:salt_bridges:occupancy_heatmap", "occupancy.parquet"),
        ("solvation:shell_rdf", "system_rdf.parquet"),
        ("shared_states:state_populations", "replicate_populations.parquet"),
        ("sdc4:condition_A:nxip_height_occlusion_joint", "timeseries.parquet"),
        ("lipid_contacts:lipid_contact_matrix_condition_A_upper_enrichment", "matrix.parquet"),
        ("diffusion:diffusion_coefficients", "replicate_coefficients.parquet"),
    ):
        record = describe_plot(identifier)
        assert record["id"] == identifier
        assert record["renderer"].startswith("reporting.")
        assert record["destination"].endswith(".png")
        assert record["status"] == "conditional_on_saved_data"
        assert "dpi" in record["available_overrides"]
        assert any(expected in source for source in record["source_candidates"])
        if identifier.startswith("sdc4:"):
            assert record["destination"].startswith("project_summary/figures/sdc4_mechanisms/")

    with pytest.raises(ValueError, match="unsupported plot ID"):
        describe_plot("sdc4:condition_A:unsupported")


def test_figure_preview_manifest_can_describe_exact_output_ids(tmp_path):
    from mdscope.reporting.plot_ids import id_from_figure_path

    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, tmp_path, repouts)
    figure = next((tmp_path / "project_summary/figures").rglob("*.png"))
    identifier = id_from_figure_path(figure)
    from mdscope.registry import figure_descriptor
    record = figure_descriptor(cfg, identifier)
    assert record["destination"] == "project_summary/figures/rmsd/secondary/time_series@rmsd_A@timeseries.png"


def test_single_time_series_does_not_read_unselected_statistics(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    stat = tmp_path / "project_summary/statistics/time_series/replicate_statistics.parquet"
    stat.write_bytes(b"invalid statistics file that must remain unread")
    cfg["plots"]["include"] = ["time_series:rmsd_A:timeseries"]
    plot_all(cfg, tmp_path, repouts)
    assert (tmp_path / "project_summary/figures/rmsd/secondary/time_series@rmsd_A@timeseries.png").exists()


def test_single_profile_does_not_read_unselected_contrast(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    contrast = tmp_path / "project_summary/statistics/profiles/system_comparisons.parquet"
    contrast.write_bytes(b"invalid contrast that must remain unread")
    cfg["plots"]["include"] = ["profiles:lipid_order__POPC__sn1:profile"]
    plot_all(cfg, tmp_path, repouts)
    assert (tmp_path / "project_summary/figures/membrane_order/primary/profiles@lipid_order__POPC__sn1@profile.png").exists()


def test_selected_occupancy_skips_event_presence_and_lifetimes(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    for reps in repouts.values():
        for rep in reps:
            source = rep / "interactions/salt_bridges"
            (source / "presence.parquet").write_bytes(b"invalid high-cardinality presence")
            (source / "lifetimes.parquet").write_bytes(b"invalid lifetimes")
    cfg["plots"]["include"] = ["interactions:salt_bridges:occupancy_heatmap"]
    plot_all(cfg, tmp_path, repouts)
    assert (tmp_path /
            "project_summary/figures/salt_bridges/primary/interactions@salt_bridges@occupancy_heatmap.png").is_file()


def test_selected_interaction_kind_honors_explicit_exclusion(tmp_path):
    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["interactions:salt_bridges:occupancy_heatmap"]
    cfg["plots"]["interaction_exclude_kinds"] = ["salt_bridges"]
    plot_all(cfg, tmp_path, repouts)
    assert not list((tmp_path / "project_summary/figures").rglob("*.png"))


def test_sdc4_condition_figures_keep_replica_traces_separate(tmp_path, monkeypatch):
    from mdscope.reporting import sdc4_condition_figures as module
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms

    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = [
        "sdc4:condition_A:nxip_height_occlusion_joint",
        "sdc4:condition_A:lipid_density_maps",
        "sdc4:condition_A:lipid_region_mean_contact_fraction",
    ]
    observed = []
    original_save = module._save

    def check_legend(fig, config, out, system, stem, **kwargs):
        if stem == "nxip_height_occlusion_joint":
            for ax in fig.axes:
                legend = ax.get_legend()
                assert legend is not None
                labels = [text.get_text() for text in legend.get_texts()]
                assert set(labels) == {"rep1", "rep2"}
                # Separate scatter collections: no pooled regression.
                assert len(ax.collections) == 2
            observed.append(True)
        return original_save(fig, config, out, system, stem, **kwargs)

    monkeypatch.setattr(module, "_save", check_legend)
    plot_sdc4_mechanisms(cfg, tmp_path, repouts)

    root = tmp_path / "project_summary/figures/sdc4_mechanisms"
    assert sorted(id_from_figure_path(path) for path in root.rglob("*.png")) == sorted(cfg['plots']['include'])
    assert observed
    assert not (tmp_path / "project_summary/figures/sdc4/condition_B").exists()
    for paths in repouts.values():
        for repout in paths:
            assert not list((repout / "extensions/sdc4_mechanisms").glob("*.png"))
            assert (repout / "extensions/sdc4_mechanisms/timeseries.parquet").is_file()


def test_condition_sdc4_missing_replica_is_visible_not_silently_averaged(tmp_path, monkeypatch):
    from mdscope.reporting import sdc4_condition_figures as module
    from mdscope.reporting.sdc4_figures import plot_sdc4_mechanisms

    cfg, repouts = _fake_data(tmp_path)
    cfg["plots"]["include"] = ["sdc4:condition_A:lipid_region_mean_contact_fraction"]
    path = repouts["condition_A"][1] / "extensions/sdc4_mechanisms/lipid_region_specificity.parquet"
    path.unlink()
    inspected = []
    original = module._save

    def check_missing(fig, config, out, system, stem, **kwargs):
        assert "1/2" in fig._suptitle.get_text()
        assert any("Not measured" in text.get_text()
                   for ax in fig.axes for text in ax.texts)
        inspected.append(stem)
        return original(fig, config, out, system, stem, **kwargs)

    monkeypatch.setattr(module, "_save", check_missing)
    plot_sdc4_mechanisms(cfg, tmp_path, repouts)
    assert inspected == ["lipid_region_mean_contact_fraction"]


def test_stale_per_replica_sdc4_png_is_preserved_but_not_active(tmp_path):
    rep = tmp_path / "condition" / "rep1"
    legacy = rep / "extensions/sdc4_mechanisms/nxip_height_occlusion_joint.png"
    legacy.parent.mkdir(parents=True)
    Image.new("RGB", (100, 80), "white").save(legacy)
    cfg = {"plots": {"enabled": True}}
    assert figure_status(legacy, cfg, tmp_path, {"condition": [rep]}) == "disabled"
    assert legacy.is_file(), "Never remove existing analysis results as a plotting side-effect"



def test_descriptive_plotting_does_not_require_completed_statistics(tmp_path, project):
    """A selected saved-data RMSD plot does not imply running bootstrap."""
    from mdscope.pipeline.runner import finalize_project, run_one

    path, root = project("descriptive_only", statistics=False)
    options = yaml.safe_load(path.read_text())
    options["plots"] = {
        "enabled": True,
        "include": ["time_series:rmsd_A:timeseries"],
    }
    path.write_text(yaml.safe_dump(options))
    run_one(path, 0)
    assert not (root / "project_summary/statistics/manifest.json").exists()

    finalize_project(path, phase="plotting")
    figures = sorted((root / "project_summary/figures").rglob("*.png"))
    assert [figure.name for figure in figures] == ["time_series@rmsd_A@timeseries.png"]
    assert not (root / "project_summary/statistics/manifest.json").exists()


def test_statistical_plot_still_requires_completed_statistics(tmp_path, project):
    """A bootstrap figure cannot quietly render without inferential inputs."""
    from mdscope.pipeline.runner import finalize_project, run_one

    path, root = project("statistical_plot_without_stats", statistics=False)
    options = yaml.safe_load(path.read_text())
    options["plots"] = {
        "enabled": True,
        "include": ["time_series:rmsd_A:bootstrap_table"],
    }
    path.write_text(yaml.safe_dump(options))
    run_one(path, 0)
    with pytest.raises(RuntimeError, match="completed statistics"):
        finalize_project(path, phase="plotting")
    assert not list(root.rglob("*.png"))


def test_exact_descriptive_lipid_matrix_can_render_without_statistics(tmp_path, project):
    from mdscope.pipeline.runner import finalize_project, run_one
    from mdscope.planning import build_plan

    path, root = project("lipid_description", statistics=False)
    options = yaml.safe_load(path.read_text())
    import MDAnalysis as mda
    from mdscope_test_support import membrane_universe

    universe = membrane_universe()
    universe.add_TopologyAttr("bonds", [(8, 9)])
    topology, trajectory = tmp_path / "membrane.pdb", tmp_path / "membrane.dcd"
    universe.atoms.write(topology)
    with mda.Writer(str(trajectory), n_atoms=len(universe.atoms)) as writer:
        for _ in universe.trajectory:
            writer.write(universe.atoms)
    for system_options in options["systems"].values():
        system_options["topology"] = str(topology)
        for replica in system_options["replicates"]:
            replica["trajectory"] = str(trajectory)
    for name in ("alignment", "rmsd", "rmsf"):
        options["analyses"][name] = {"enabled": False}
    options["trajectory"] = {"start": 0, "stop": 5, "stride": 1, "chunk_frames": 2}
    options["analyses"]["lipid_contacts"] = {
        "enabled": True, "selection": "resname POPC",
        "species_selections": {"POPC": "resname POPC and name P"},
    }
    system = next(iter(options["systems"]))
    identifier = f"lipid_contacts:lipid_contact_matrix_{system}_upper_neighbors_per_molecule"
    options["plots"] = {"enabled": True, "include": [identifier]}
    path.write_text(yaml.safe_dump(options))
    run_one(path, 0)
    finalize_project(path, phase="plotting")
    from mdscope.registry import figure_descriptor
    assert (root / figure_descriptor(options, identifier)["destination"]).is_file()
    assert not (root / "project_summary/statistics/manifest.json").exists()
    plan = build_plan(load_config(path), path)
    figure = next(f for family in plan["figure_families"] for f in family["specific_figures"] if f["id"] == identifier)
    assert not figure["requires_statistics"]
