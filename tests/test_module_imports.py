"""Import smoke checks for the canonical MDscope modules."""
def test_canonical_modules_import():
    from mdscope.pipeline import streaming
    from mdscope import registry
    from mdscope.observables import membrane, density_1d
    from mdscope.observables import rmsd, rmsf
    from mdscope.observables import dssp, sasa, membrane_order
    from mdscope.observables import radius_of_gyration
    from mdscope.observables import residue_membrane_distance
    from mdscope.core import alignment, lipid_geometry, profiling
    from mdscope.observables import msd, bilayer_normal, lipid_enrichment
    from mdscope.pipeline.runner import run_all, run_one, finalize_project
    from mdscope.cli import main
    assert callable(streaming.run_streaming_pass)
    assert callable(msd.msd_metrics)
    assert callable(bilayer_normal.fit_local_normals)
    assert callable(lipid_enrichment.enrichment_frame)
    assert callable(alignment._fit_transform)
    assert callable(lipid_geometry.anchor_groups)
    assert callable(profiling.PhaseProfiler)
    assert callable(registry.specialized_tasks)
    assert all(callable(f) for f in (run_all, run_one, finalize_project, main))
