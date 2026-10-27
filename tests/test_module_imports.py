"""Smoke checks that all public MDscope entrypoints still import after module extraction."""

def test_module_exports_and_compatibility():
    from mdscope.analysis import timeseries, profiles, interactions
    from mdscope.analysis import (streaming,rmsd,rmsf,msd,membrane,density_1d,
        alignment,radius_of_gyration,dssp,sasa,residue_membrane_distance,
        salt_bridges,hydrogen_bonds,ligand_contacts,interaction_tables,
        membrane_order,local_environment,registry)
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
    assert callable(registry.specialized_tasks)
