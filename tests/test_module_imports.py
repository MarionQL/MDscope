"""Smoke checks that all public MDscope entrypoints still import after module extraction."""

def test_module_exports_and_compatibility():
    from mdscope.analysis import timeseries, profiles, interactions
    from mdscope.analysis import (streaming,rmsd,rmsf,msd,membrane,density_1d,
        alignment,radius_of_gyration,dssp,sasa,residue_membrane_distance,
        salt_bridges,hydrogen_bonds,ligand_contacts,interaction_tables,
        membrane_order,bilayer_normal,local_environment,registry)
    from mdscope.analysis.runner import run_all,run_one,finalize_project
    from mdscope.analysis.cli import main
    from mdscope.analysis.statistics import analyze_project
    from mdscope.analysis.plotting import plot_all

    assert timeseries.run_streaming_pass is streaming.run_streaming_pass
    assert timeseries.Welford is rmsf.Welford
    assert timeseries.temporal_rmsf_profiles is rmsf.temporal_rmsf_profiles
    assert timeseries.temporal_density_profiles is density_1d.temporal_density_profiles
    assert profiles.run_residue_membrane_distance is residue_membrane_distance.run_residue_membrane_distance
    assert interactions.run_hydrogen_bonds is hydrogen_bonds.run_hydrogen_bonds
    from mdscope.analysis.observables import msd as nested_msd
    from mdscope.analysis.observables import bilayer_normal as nested_normals
    from mdscope.analysis.observables import lipid_enrichment as nested_enrichment
    from mdscope.analysis.reporting import coupling_screen as nested_coupling
    from mdscope.analysis import lipid_enrichment, coupling_screen
    assert msd.msd_metrics is nested_msd.msd_metrics
    assert bilayer_normal.fit_local_normals is nested_normals.fit_local_normals
    assert lipid_enrichment.enrichment_frame is nested_enrichment.enrichment_frame
    assert coupling_screen.screen_couplings is nested_coupling.screen_couplings
    assert callable(registry.specialized_tasks)
