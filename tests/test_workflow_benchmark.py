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
