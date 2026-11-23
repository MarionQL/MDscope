import json
from pathlib import Path
from runpy import run_path

summarize = run_path(str(Path(__file__).resolve().parents[1] / "scripts" / "benchmark_workflow.py"))["summarize"]


def test_workflow_benchmark_aggregates_nonoverlapping_scopes(tmp_path):
    rep=tmp_path/"carb"/"rep1"
    rep.mkdir(parents=True)
    (rep/"streaming_phase_profile.json").write_text(json.dumps({
        "frames_completed":20,
        "timings":[{"stage":"trajectory_decode","seconds":2.5,
                    "seconds_per_frame":.125}]}))
    (rep/"stage_timings.json").write_text(json.dumps({
        "stages":[{"stage":"streaming","seconds":10.0}]}))
    out=tmp_path/"project_summary"
    out.mkdir()
    (out/"finalization_timings.json").write_text(json.dumps({
        "stages":[{"stage":"pdf_and_bundle","seconds":3.0}]}))
    rows=summarize(tmp_path)
    assert len(rows)==3
    assert rows[0]["stage"]=="streaming"
    assert {row["scope"] for row in rows}=={
        "streaming_phase","replicate_stage","finalization"}


def test_performance_benchmark_rejects_short_inputs_before_analysis(tmp_path):
    import os
    import subprocess
    import sys

    import MDAnalysis as mda

    source = tmp_path / 'input'
    source.mkdir()
    u = mda.Universe.empty(1, trajectory=True)
    u.atoms.positions = [[1., 2., 3.]]
    for condition in ('carb', 'nocarb'):
        with mda.Writer(str(source / f'{condition}-sdc4-rep1_test.xtc'), n_atoms=1) as writer:
            writer.write(u.atoms)
            writer.write(u.atoms)
    script = Path(__file__).resolve().parent / 'real_data/performance_pipeline_benchmark.py'
    result = subprocess.run(
        [sys.executable, str(script), '--frames', '3', '--draws', '2'],
        env={**os.environ, 'MDSCOPE_REAL_DATA': str(source),
             'MDSCOPE_REAL_OUTPUT': str(tmp_path / 'results')},
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 2
    assert 'requested 3 frames but only 2 are available' in result.stderr
    assert not list((tmp_path / 'results').glob('benchmark_frames_*'))
