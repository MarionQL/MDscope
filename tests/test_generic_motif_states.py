"""YAML motifs feed saved-product state fitting without an SDC4 analysis."""
import json
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from test_solvation_motifs import _protein

from mdscope.config import resolve_config
from mdscope.mechanisms.shared_states import fit_shared_states
from mdscope.observables.solvation import _targets


def _spec():
    return {"enabled": True, "preset": None, "features": ["water", "size"],
            "feature_sources": {
                "water": {"source": "solvation", "group": "solvation/water/motif/patch",
                          "metric": "mean_n_solvent", "units": "molecules"},
                "size": {"source": "summary", "group": None,
                         "metric": "rmsd_A", "units": "angstrom"}},
            "n_states": 2, "bootstrap_draws": 20}


def _saved(root, mismatch=False):
    rep = root / "sample" / "rep1"
    (rep / "summary").mkdir(parents=True)
    (rep / "extensions/solvation").mkdir(parents=True)
    times = np.arange(20, dtype=float)*10
    x = np.r_[np.linspace(1, 2, 10), np.linspace(8, 9, 10)]
    pd.DataFrame({"time_ps": times, "rmsd_A": x}).to_parquet(rep / "summary/timeseries.parquet")
    pd.DataFrame({"time_ps": times + int(mismatch), "solvent": "water",
                  "target_class": "motif", "site_class": "patch",
                  "mean_n_solvent": x*2}).to_parquet(rep / "extensions/solvation/species_timeseries.parquet")
    return {"sample": [rep]}


def test_generic_saved_features_fit_without_sdc4(tmp_path):
    paths = _saved(tmp_path)
    out = fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    assignments = pd.read_parquet(out / "frame_assignments.parquet")
    assert assignments.state.tolist() == [0]*10 + [1]*10
    model = json.loads((out / "model.json").read_text())
    assert model["feature_sources"] == _spec()["feature_sources"]
    assert model["alignment"] == "exact"
    assert not (out / "state_conditioned_nxip.parquet").exists()


def test_mismatched_feature_times_preserve_existing_outputs(tmp_path):
    paths = _saved(tmp_path, mismatch=True)
    out = tmp_path / "project_summary/shared_states"
    out.mkdir(parents=True)
    (out / "model.json").write_text("previous valid model")
    with pytest.raises(ValueError, match="time grid"):
        fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    assert (out / "model.json").read_text() == "previous valid model"
    assert list(out.iterdir()) == [out / "model.json"]


def test_shared_motif_references_resolve_chain_scoped_whole_residues():
    definitions = {"acidic": {"sequence": "DE", "selection": "protein",
                              "chains": ["PROB"], "overlap": False, "missing": "error"}}
    targets = _targets(_protein(), {"targets": {"motif": {"patch": {"motif": "acidic"}}}},
                       motifs=definitions)
    assert len(targets) == 1
    assert targets[0][3].resindices.tolist() == [2, 2, 3, 3]


def test_shared_config_is_validated_and_roundtrips():
    import yaml

    from mdscope.studio.model import import_document, new_document
    doc = new_document()
    doc["selections"]["motifs"] = {"acidic": {"sequence": "DE"}}
    doc["analyses"]["shared_states"] = _spec()
    resolved = resolve_config(doc, check_paths=False)
    assert resolved["analyses"]["shared_states"]["preset"] is None
    restored = import_document(yaml.safe_dump(doc))
    assert restored["selections"]["motifs"] == doc["selections"]["motifs"]
    bad = deepcopy(doc)
    bad["analyses"]["shared_states"]["features"] = []
    with pytest.raises(ValueError, match="nonempty"):
        resolve_config(bad, check_paths=False)


@pytest.mark.parametrize("change", [
    {"preset": "unknown"}, {"alignment": "nearest"}, {"features": []},
    {"features": ["water", "water"]}, {"features": ["time_ps"]},
    {"min_complete_fraction": 0}, {"min_complete_fraction": float("nan")},
    {"block_frames": 0}, {"n_states": True},
])
def test_invalid_state_choices_fail_at_yaml_boundary(change):
    from mdscope.studio.model import new_document
    doc = new_document()
    doc["analyses"]["shared_states"] = {**_spec(), **change}
    with pytest.raises(ValueError, match="shared_states"):
        resolve_config(doc, check_paths=False)


@pytest.mark.parametrize("change", [{"source": "arbitrary_file"}, {"units": ""},
                                   {"group": "wrong_identity"}, {"metric": ""}])
def test_invalid_feature_identity_is_rejected(change):
    from mdscope.state_feature_spec import validate_state_features
    spec = _spec()
    spec["feature_sources"]["water"].update(change)
    with pytest.raises(ValueError, match="feature water"):
        validate_state_features(spec)


def test_feature_gaps_remain_unavailable_and_disable_temporal_ci(tmp_path):
    paths = _saved(tmp_path)
    path = paths["sample"][0] / "summary/timeseries.parquet"
    table = pd.read_parquet(path)
    table.loc[9, "rmsd_A"] = np.nan
    table.to_parquet(path)
    out = fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    diagnostics = pd.read_parquet(out / "sampling_diagnostics.parquet").iloc[0]
    assert diagnostics.missing_feature_frames == 1
    assert diagnostics.complete_feature_fraction == .95
    assert diagnostics.observed_transitions == 0
    assert diagnostics.transitions_across_missing_gaps == 1
    assert pd.read_parquet(out / "replicate_populations.parquet").within_replicate_block_ci_low.isna().all()


def test_generic_rerun_archives_old_sdc4_outcomes(tmp_path):
    paths = _saved(tmp_path)
    old = tmp_path / "project_summary/shared_states"
    old.mkdir(parents=True)
    (old / "state_conditioned_nxip.parquet").write_bytes(b"previous outcome")
    out = fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    assert not (out / "state_conditioned_nxip.parquet").exists()
    assert (out.parent / "archive/shared_states/state_conditioned_nxip.parquet").read_bytes() == b"previous outcome"


@pytest.mark.parametrize("problem", ["absent", "constant", "all_missing", "duplicate_time"])
def test_invalid_saved_feature_does_not_replace_a_valid_model(tmp_path, problem):
    paths = _saved(tmp_path)
    path = paths["sample"][0] / "summary/timeseries.parquet"
    table = pd.read_parquet(path)
    if problem == "absent":
        table = table.drop(columns="rmsd_A")
    elif problem == "constant":
        table["rmsd_A"] = 1.
    elif problem == "all_missing":
        table["rmsd_A"] = np.nan
    else:
        table.loc[1, "time_ps"] = 0.
    table.to_parquet(path)
    out = tmp_path / "project_summary/shared_states"
    out.mkdir(parents=True)
    (out / "model.json").write_text("existing model")
    with pytest.raises(ValueError):
        fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    assert (out / "model.json").read_text() == "existing model"


def test_shared_motif_policies_and_topology_identity():
    from mdscope.core.motifs import resolve_motifs
    u = _protein()
    u.residues.resids = [1, 2, 1, 2]
    matches = resolve_motifs(u, {"patch": {"sequence": "DE"}})["patch"]
    assert list(matches) == ["patch_PROA_1", "patch_PROB_1"]
    assert [row["resindices"] for row in matches.values()] == [[0, 1], [2, 3]]
    assert all(row["resids"] == [1, 2] for row in matches.values())
    assert resolve_motifs(u, {"patch": {"sequence": "DE", "chains": [], "missing": "allow"}}) == {"patch": {}}
    with pytest.raises(ValueError, match="no match"):
        resolve_motifs(u, {"patch": {"sequence": "AA"}})
    u.residues.resnames = ["ALA"]*4
    u.atoms.chainIDs = ["A"]*8
    nonoverlap = resolve_motifs(u, {"repeat": {"pattern": "AA", "overlap": False}})
    overlap = resolve_motifs(u, {"repeat": {"pattern": "AA", "overlap": True}})
    assert len(nonoverlap["repeat"]) == 2
    assert len(overlap["repeat"]) == 3
    with pytest.raises(ValueError, match="overlaps"):
        _targets(u, {"targets": {"motif": {"repeat": {"motif": "repeat"}}}},
                 motifs={"repeat": {"pattern": "AA", "overlap": True}})
    # A selection gap or unrecognized residue must never join distant residues.
    assert resolve_motifs(u, {"gap": {"sequence": "AAA", "selection": "resindex 0 2 3", "missing": "allow"}}) == {"gap": {}}
    u.residues[1].resname = "UNK"
    assert resolve_motifs(u, {"unknown": {"pattern": "A.A", "missing": "allow"}}) == {"unknown": {}}


@pytest.mark.parametrize("definition", [{}, {"sequence": "DE", "pattern": "D.E"},
                                       {"sequence": "DX"}, {"pattern": "["},
                                       {"pattern": "A*"}, {"sequence": "DE", "overlap": "yes"},
                                       {"sequence": "DE", "missing": "zero"}])
def test_invalid_motif_definitions_are_rejected(definition):
    from mdscope.core.motifs import validate_motifs
    with pytest.raises(ValueError, match="motif"):
        validate_motifs({"patch": definition})


def test_motifs_preview_through_production_topology_preflight(tmp_path):
    from mdscope.planning import authoring_catalog, build_plan
    from mdscope.studio.model import new_document, topology_preflight
    u = _protein()
    # PDB uses one-character chain labels; the topology indices remain authoritative.
    u.atoms.chainIDs = ["A"]*4 + ["B"]*4
    topology = tmp_path / "protein.pdb"
    u.atoms.write(topology)
    doc = new_document()
    doc["systems"]["condition1"]["topology"] = str(topology)
    doc["selections"]["motifs"] = {"patch": {"sequence": "DE"}}
    doc["analyses"]["shared_states"] = _spec()
    cfg = resolve_config(doc, check_paths=False)
    counts = topology_preflight(doc, tmp_path)
    assert list(counts["condition1"]["motifs"]["patch"]) == ["patch_A_1", "patch_B_1"]
    plan = build_plan(cfg, tmp_path / "analysis.yaml")
    state = next(item for item in plan["analyses"] if item["analysis"] == "shared_states")
    assert state["state_model_inputs"]["feature_sources"] == _spec()["feature_sources"]
    assert authoring_catalog(cfg)["motifs"] == doc["selections"]["motifs"]
    assert not (tmp_path / "mdscope_results").exists()


def test_maintained_sdc4_configs_roundtrip_without_migration():
    from pathlib import Path

    import yaml

    from mdscope.studio.model import import_document
    for name in ("synd4_baseline.yaml", "synd4_full_mechanisms.yaml", "synd4_full_mechanisms_finalize_fast.yaml"):
        text = (Path(__file__).resolve().parents[1] / name).read_text()
        original = yaml.safe_load(text)
        restored = import_document(yaml.safe_dump(import_document(text)))
        assert resolve_config(restored, check_paths=False) == resolve_config(original, check_paths=False)


def test_yaml_motif_measurement_to_saved_state_finalize_and_figures(tmp_path):
    import MDAnalysis as mda
    import yaml

    from mdscope.pipeline.runner import finalize_project, run_all
    from mdscope.studio.model import new_document
    u = mda.Universe.empty(6, n_residues=4, atom_resindex=[0, 0, 1, 1, 2, 3], trajectory=True)
    u.add_TopologyAttr("names", ["CA", "CB", "CA", "CB", "OW", "OW"])
    u.add_TopologyAttr("resnames", ["ASP", "GLU", "SOL", "SOL"])
    u.add_TopologyAttr("resids", [1, 2, 1, 2])
    u.add_TopologyAttr("chainIDs", ["A"]*4 + ["W"]*2)
    u.atoms.positions = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [20, 20, 20]]
    u.dimensions = [60, 60, 60, 90, 90, 90]
    top, traj = tmp_path / "protein.pdb", tmp_path / "protein.xtc"
    u.atoms.write(top)
    with mda.Writer(str(traj), n_atoms=6, dt=10) as writer:
        for frame in range(20):
            u.atoms[5].position = [20, 20, 20] if frame < 10 else [1, 1, 1]
            u.trajectory.ts.frame = frame
            u.trajectory.ts.time = frame*10
            writer.write(u)
    doc = new_document()
    doc["project"]["output"] = str(tmp_path / "results")
    doc["systems"] = {"protein": {"topology": str(top), "replicates": [{"label": "rep1", "trajectory": str(traj)}]}}
    doc["trajectory"] = {"stride": 1, "chunk_frames": 7}
    doc["selections"]["motifs"] = {"patch": {"sequence": "DE", "chains": ["A"]}}
    doc["analyses"]["solvation"].update(enabled=True, water_selection="resname SOL and name OW",
        targets={"motif": {"patch": {"motif": "patch"}}})
    state = _spec()
    state["features"] = ["water"]
    state["feature_sources"].pop("size")
    doc["analyses"]["shared_states"] = state
    doc["plots"] = {"enabled": True, "include": ["shared_states:state_populations", "shared_states:state_sampling_diagnostics"]}
    path = tmp_path / "analysis.yaml"
    path.write_text(yaml.safe_dump(doc))
    root = run_all(path)
    counts = pd.read_parquet(root / "protein/rep1/extensions/solvation/counts.parquet")
    assert counts.n_solvent.tolist() == [1]*10 + [2]*10
    assigned = pd.read_parquet(root / "project_summary/shared_states/frame_assignments.parquet")
    assert assigned.state.tolist() == [0]*10 + [1]*10
    metadata = json.loads((root / "protein/rep1/extensions/solvation/metadata.json").read_text())
    assert metadata["motifs"]["patch"]["patch_A_1"]["resindices"] == [0, 1]
    finalize_project(path, phase="plotting")
    assert (root / "project_summary/shared_states/state_populations.png").is_file()
    assert (root / "project_summary/shared_states/state_sampling_diagnostics.png").is_file()


def test_replica_balanced_scaling_ignores_excess_length(tmp_path):
    paths = _saved(tmp_path)
    first = paths["sample"][0]
    second = first.parent / "rep2"
    import shutil
    shutil.copytree(first, second)
    for relative in ("summary/timeseries.parquet", "extensions/solvation/species_timeseries.parquet"):
        table = pd.read_parquet(second / relative)
        extra = table.copy()
        extra["time_ps"] += 200
        column = "rmsd_A" if "rmsd_A" in table else "mean_n_solvent"
        extra[column] += 10
        pd.concat([table, extra], ignore_index=True).to_parquet(second / relative)
    paths["sample"].append(second)
    out = fit_shared_states({"analyses": {"shared_states": _spec()}}, tmp_path, paths)
    model = json.loads((out / "model.json").read_text())
    assert model["training_frames_per_replicate"] == 20
    first_x = pd.read_parquet(first / "summary/timeseries.parquet").rmsd_A.to_numpy()
    second_x = pd.read_parquet(second / "summary/timeseries.parquet").rmsd_A.to_numpy()
    balanced = np.r_[first_x, second_x[np.linspace(0, 39, 20).astype(int)]]
    assert model["scaler_median"][1] == np.median(balanced)
    assert model["scaler_sd"][1] == pytest.approx(np.std(balanced))
    assert [row["complete_frames"] for row in model["replicate_time_windows"]] == [20, 40]


def test_generic_motif_overlaps_gaps_modifications_and_segments():
    import MDAnalysis as mda

    from mdscope.core.motifs import resolve_motifs
    u = mda.Universe.empty(6, n_residues=6, n_segments=2, atom_resindex=np.arange(6),
                           residue_segindex=[0, 0, 0, 1, 1, 1], trajectory=True)
    u.add_TopologyAttr("names", ["CA"]*6)
    u.add_TopologyAttr("resnames", ["ALA"]*6)
    u.add_TopologyAttr("resids", [1, 2, 3, 1, 2, 3])
    u.add_TopologyAttr("chainIDs", ["A"]*6)
    u.add_TopologyAttr("segids", ["first", "second"])
    spec = {"pattern": "AA", "overlap": True}
    matches = resolve_motifs(u, {"pair": spec})["pair"]
    assert list(matches) == ["pair_A_seg0_1", "pair_A_seg0_2", "pair_A_seg1_1", "pair_A_seg1_2"]
    assert [row["resindices"] for row in matches.values()] == [[0, 1], [1, 2], [3, 4], [4, 5]]
    assert len(resolve_motifs(u, {"pair": {**spec, "overlap": False}})["pair"]) == 2
    assert resolve_motifs(u, {"pair": {**spec, "selection": "index 0 2 3 5", "missing": "allow"}}) == {"pair": {}}
    u.residues[1].resname = "UNK"
    assert resolve_motifs(u, {"triple": {"pattern": "A.A", "selection": "segid first", "missing": "allow"}}) == {"triple": {}}
    with pytest.raises(ValueError, match="overlaps"):
        _targets(u, {"targets": {"motif": {"pair": {"motif": "pair"}}}}, motifs={"pair": spec})


def test_disabled_states_do_not_read_inputs_or_create_outputs(tmp_path):
    assert fit_shared_states({"analyses": {"shared_states": {"enabled": False}}}, tmp_path, {}) is None
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("definition", [None, [], {"x": {"sequence": "DE", "pattern": "D.E"}},
    {"x": {"pattern": "["}}, {"x": {"pattern": ".*"}}, {"x": {"sequence": "de"}},
    {"x": {"sequence": "DE", "chains": "A"}}, {"x": {"sequence": "DE", "missing": "zero"}},
    {"x": {"sequence": "DE", "overlap": "yes"}}, {"x": {"sequence": "DE", "selection": ""}}])
def test_invalid_motif_definitions_fail_before_topology(definition):
    from mdscope.core.motifs import validate_motifs
    with pytest.raises(ValueError):
        validate_motifs(definition)


@pytest.mark.parametrize("pattern", [r"(A)\1", "(?i)AA"])
def test_overlap_preserves_regex_captures_and_flags(pattern):
    from mdscope.core.motifs import resolve_motifs
    u = _protein()
    u.residues.resnames = ["ALA"]*4
    assert len(resolve_motifs(u, {"pair": {"pattern": pattern, "overlap": True}})["pair"]) == 2


def test_zero_width_regex_fails_during_topology_preflight():
    from mdscope.core.motifs import resolve_motifs
    u = _protein()
    with pytest.raises(ValueError, match="consume residues"):
        resolve_motifs(u, {"zero": {"pattern": "(?=D)", "missing": "allow", "overlap": True}})


def test_chain_scope_uses_nonblank_segment_fallback():
    from mdscope.core.motifs import resolve_motifs
    u = _protein()
    u.add_TopologyAttr("segids", ["PROA"])
    u.atoms.chainIDs = [" "]*len(u.atoms)
    result = resolve_motifs(u, {"patch": {"sequence": "DE", "chains": ["PROA"]}})["patch"]
    assert len(result) == 2
    assert all(item["chain"] == "PROA" for item in result.values())


def test_studio_group_rename_migrates_state_inputs_and_delete_is_atomic():
    from mdscope.studio.model import change_group, new_document
    doc = new_document()
    doc["selections"]["groups"] = {"old": {"selection": "protein"}}
    spec = _spec()
    spec["feature_sources"]["size"]["group"] = "old"
    doc["analyses"]["shared_states"] = spec
    change_group(doc, "old", "new")
    assert doc["analyses"]["shared_states"]["feature_sources"]["size"]["group"] == "new"
    before = deepcopy(doc)
    with pytest.raises(ValueError, match="shared_states features"):
        change_group(doc, "new")
    assert doc == before


def test_production_finalize_failure_preserves_model_and_saved_configuration(tmp_path):
    from mdscope.pipeline.runner import finalize_project
    from mdscope.studio.model import new_document
    _saved(tmp_path, mismatch=True)
    old = tmp_path / "project_summary/shared_states"
    old.mkdir(parents=True)
    (old / "model.json").write_text("previous valid model")
    (tmp_path / "resolved_config.yaml").write_text("previous valid configuration")
    doc = new_document()
    doc["systems"] = {"sample": {"topology": "not_opened.pdb", "replicates": [{"label": "rep1", "trajectory": "not_opened.xtc"}]}}
    doc["project"]["output"] = str(tmp_path)
    doc["analyses"]["shared_states"] = _spec()
    with pytest.raises(ValueError, match="time grid"):
        finalize_project(tmp_path / "analysis.yaml", cfg=resolve_config(doc, check_paths=False))
    assert (old / "model.json").read_text() == "previous valid model"
    assert (tmp_path / "resolved_config.yaml").read_text() == "previous valid configuration"
    assert not (tmp_path / ".archives").exists()
    assert not list(tmp_path.glob(".shared-state-preflight-*"))


def test_state_feature_catalog_exposes_sources_identities_and_prerequisites():
    from mdscope.output_paths import SAVED_SCALAR_SOURCE_FILES
    from mdscope.planning import authoring_catalog, build_plan
    from mdscope.studio.model import new_document
    doc = new_document()
    doc["selections"]["motifs"] = {"patch": {"sequence": "DE"}}
    doc["analyses"]["solvation"].update(enabled=True, water_selection="resname SOL",
        ion_selections={"potassium": "resname POT"},
        targets={"motif": {"patch": {"motif": "patch"}}})
    doc["analyses"]["shared_states"] = _spec()
    cfg = resolve_config(doc, check_paths=False)
    catalog = authoring_catalog(cfg)["state_features"]
    assert any(row["group"] == "solvation/water/motif/patch" and row["metric"] == "mean_n_solvent"
               for row in catalog["candidates"])
    assert any(row["group"] == "solvation/ions/potassium/motif/patch"
               for row in catalog["candidates"])
    assert set(catalog["sources"]) == set(SAVED_SCALAR_SOURCE_FILES)
    state = next(row for row in build_plan(cfg, "analysis.yaml")["analyses"] if row["analysis"] == "shared_states")
    prerequisite = state["state_model_inputs"]["source_prerequisites"]["water"]
    assert prerequisite["paths_within_replica"] == ["extensions/solvation/species_timeseries.parquet"]
    assert prerequisite["identity_status"] == "template_until_saved_product_validation"
