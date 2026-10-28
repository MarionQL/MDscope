"""Profile accounting and SLURM launch contracts."""
import json
from pathlib import Path

from mdscope.analysis.profiling import PhaseProfiler


def test_profiler_is_opt_in_and_produces_correct_stream(tmp_path, monkeypatch):
    monkeypatch.delenv("MDSCOPE_PROFILE", raising=False)
    assert not PhaseProfiler.from_environment().enabled
    monkeypatch.setenv("MDSCOPE_PROFILE", "1")
    profiler = PhaseProfiler.from_environment()
    assert profiler.enabled
    assert list(profiler.frames(iter(range(4)))) == [0, 1, 2, 3]
    assert profiler.call("square", lambda value: value * value, 4) == 16
    assert profiler.calls["trajectory_decode"] == 4
    assert profiler.calls["square"] == 1
    assert profiler.seconds["square"] >= 0.
    output = tmp_path / "streaming_phase_profile.json"
    profiler.write(output, frames=4)
    data = json.loads(output.read_text())
    assert data["frames_completed"] == 4
    assert {row["stage"] for row in data["timings"]} == {
        "trajectory_decode", "square",
    }


def test_slurm_launcher_routes_progress_to_stdout_and_loads_checkout():
    script = (Path(__file__).resolve().parents[1] / "run_mdscope.sh").read_text()
    assert "exec 2>&1" in script
    assert "PYTHONUNBUFFERED=1" in script
    assert 'MDSCOPE_PROGRESS_EVERY_PERCENT' in script
    assert 'MDSCOPE_PROFILE' in script
    assert "synd4_full_mechanisms.yaml" in script
    assert "python -u -m MDscope.analysis.cli replicate" in script
    assert "MDscope commit:" in script
