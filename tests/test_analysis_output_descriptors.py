"""Static output descriptors stay authoritative across registry and producers (#260)."""
import ast
import subprocess
import sys
from pathlib import Path

import pytest

from mdscope.contracts.outputs import (
    FUSED_OUTPUTS,
    OUTPUT_DESCRIPTORS,
    RUN_SCOPED_OUTPUTS,
)
from mdscope.registry import FusedOutput, artifact_schema, fused_output, required_artifacts


def test_each_persisted_output_has_one_descriptor():
    assert set(OUTPUT_DESCRIPTORS) == set(FUSED_OUTPUTS) | set(RUN_SCOPED_OUTPUTS)
    assert all(name == descriptor.name for name, descriptor in OUTPUT_DESCRIPTORS.items())
    assert all(descriptor.primary_path or descriptor.named_run
               for descriptor in OUTPUT_DESCRIPTORS.values())
    assert FUSED_OUTPUTS["lipid_enrichment"] == "extensions/lipid_footprint/timeseries.parquet"
    assert FUSED_OUTPUTS["sdc4_mechanisms"] == "extensions/sdc4_mechanisms/timeseries.parquet"
    assert RUN_SCOPED_OUTPUTS["sasa"] == ("profiles/sasa", "residue_sasa.parquet")
    assert RUN_SCOPED_OUTPUTS["residue_membrane_distance"] == (
        "profiles/residue_membrane_distance", "distance_matrix_A.npy"
    )


@pytest.mark.parametrize("name,companions", [
    ("membrane_order", ("site_mapping.parquet", "profile_index.parquet")),
    ("lipid_contacts", ("matrix.parquet", "selection_audit.json")),
    ("solvation", ("rdf.parquet", "target_residues.parquet", "metadata.json",
                   "species_timeseries.parquet", "target_summary.parquet",
                   "shell_candidates.parquet")),
    ("hydration", ("aggregate_counts.parquet", "target_metadata.parquet")),
    ("sdc4_mechanisms", ("selection_audit.json", "lipid_density_maps.npz",
                         "hs_conditioned_solvation.parquet")),
])
def test_registry_reads_descriptor_companions(tmp_path, name, companions):
    entry = fused_output(tmp_path, name)
    cfg = {"analyses": {name: {"enabled": True}}}
    expected = (entry.path, *(entry.path.parent / item for item in companions))
    assert OUTPUT_DESCRIPTORS[name].companions == companions
    assert required_artifacts(entry, cfg) == expected


@pytest.mark.parametrize("name,filename,columns", [
    ("solvation", "counts.parquet", {"frame", "time_ps", "target_id", "n_solvent"}),
    ("hydration", "counts.parquet", {"frame", "time_ps", "compartment", "target_id", "water_count"}),
    ("hydration", "glycan_attachments.parquet", {"chain", "protein_atom_index", "n_sugars"}),
    ("membrane_order", "site_mapping.parquet", {"resname", "tail", "carbon_position", "atom_name"}),
    ("lipid_contacts", "matrix.parquet", {"leaflet", "species_a", "species_b", "enrichment"}),
    ("sdc4_mechanisms", "lipid_density_maps.npz", None),
])
def test_descriptor_schema_is_registry_schema(tmp_path, name, filename, columns):
    entry = FusedOutput(name, tmp_path / "primary.parquet")
    expected = OUTPUT_DESCRIPTORS[name].schema(tmp_path / filename)
    assert artifact_schema(entry, tmp_path / filename) == expected
    if columns is None:
        assert {"counts", "density_per_nm2", "nframes"} <= set(expected["array_keys"])
    else:
        assert columns <= set(expected["columns"])


def test_producers_use_exact_descriptor_column_contracts():
    from mdscope.observables.hydration import HydrationObserver, ATTACHMENTS, OUTPUT
    from mdscope.mechanisms.residue_occlusion import OCCLUSION_FILE
    from mdscope.observables.solvation import COUNT_COLUMNS
    from mdscope.observables.solvation_common import EPISODE_COLUMNS
    from mdscope.contracts.outputs import (
        GLYCAN_ATTACHMENTS_FILE, HS_SUGAR_HYDRATION_PATH,
        HYDRATION_COUNT_COLUMNS, RESIDUE_OCCLUSION_FILE,
        SOLVATION_COUNT_COLUMNS, SOLVATION_EPISODE_COLUMNS,
    )
    assert HydrationObserver.COUNT_COLUMNS is HYDRATION_COUNT_COLUMNS
    assert COUNT_COLUMNS is SOLVATION_COUNT_COLUMNS
    assert EPISODE_COLUMNS is SOLVATION_EPISODE_COLUMNS
    assert ATTACHMENTS == GLYCAN_ATTACHMENTS_FILE
    assert OUTPUT.as_posix() == HS_SUGAR_HYDRATION_PATH
    assert OCCLUSION_FILE == RESIDUE_OCCLUSION_FILE


def test_contract_module_is_dependency_free():
    """Static output metadata must not import MDAnalysis or a higher mdscope layer."""
    path = Path(__file__).resolve().parents[1] / "src/mdscope/contracts/outputs.py"
    module = ast.parse(path.read_text(encoding="utf-8"))
    imports = [node for node in ast.walk(module) if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert not any(isinstance(node, ast.ImportFrom) and node.level for node in imports)
    assert not any(isinstance(node, ast.Import) and
                   any(alias.name.startswith(("mdscope", "MDAnalysis", "pyarrow"))
                       for alias in node.names) for node in imports)
    code = ("import sys; from mdscope.contracts.outputs import OUTPUT_DESCRIPTORS; "
            "assert 'mdscope.registry' not in sys.modules; "
            "assert 'mdscope.observables.hydration' not in sys.modules; "
            "assert 'mdscope.observables.solvation' not in sys.modules")
    result = subprocess.run([sys.executable, "-c", code], text=True,
                            capture_output=True, timeout=45, check=False)
    assert result.returncode == 0, result.stderr


def test_named_run_output_identity_still_matches_historical_paths(tmp_path):
    assert fused_output(tmp_path, "sasa", run_id="default").path == (
        tmp_path / "profiles/sasa/residue_sasa.parquet")
    assert fused_output(tmp_path, "sasa", run_id="named").path == (
        tmp_path / "profiles/sasa/named/residue_sasa.parquet")
    assert fused_output(tmp_path, "residue_membrane_distance", run_id="default").path == (
        tmp_path / "profiles/residue_membrane_distance/default/distance_matrix_A.npy")
