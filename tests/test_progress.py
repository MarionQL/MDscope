"""Batch progress should report completed work, ETA, and bounded cadence."""
import io
from threading import Event

import pytest

from mdscope.pipeline import progress as mod


class _Stderr(io.StringIO):
    def isatty(self):
        return False


def test_slurm_reports_completed_frames_elapsed_and_eta(monkeypatch):
    stderr = _Stderr()
    monkeypatch.setattr(mod.sys, "stderr", stderr)
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_FRAMES", "25")
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_PERCENT", "1")
    clock = [0.]
    monkeypatch.setattr(mod, "perf_counter", lambda: clock[0])

    def work():
        for i in range(100):
            yield i

    for _ in mod.track(work(), total=100, desc="system/rep1"):
        clock[0] += 2.
    lines = stderr.getvalue().splitlines()
    assert len(lines) == 101  # percentage-only: start plus 1..100
    assert "0% (0/100)" in lines[0]
    assert "25.0% (25/100)" in lines[25]
    assert "elapsed 00:00:50" in lines[25]
    assert "0.50 frames/s" in lines[25]
    assert "ETA 00:02:30" in lines[25]
    assert "100.0% (100/100)" in lines[-1]
    assert "ETA 00:00:00" in lines[-1]


def test_slurm_reports_each_one_percent_even_with_large_frame_interval(monkeypatch):
    stderr = _Stderr()
    monkeypatch.setattr(mod.sys, "stderr", stderr)
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_FRAMES", "10000")
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_PERCENT", "1")
    list(mod.track(range(1900), total=1900, desc="rep"))
    lines = stderr.getvalue().splitlines()
    assert len(lines) == 101  # start and each 1% up to 100%
    assert "(19/1900)" in lines[1]
    assert "(190/1900)" in lines[10]
    assert "(1900/1900)" in lines[-1]


def test_slurm_unknown_total_and_disabled(monkeypatch):
    stderr = _Stderr()
    monkeypatch.setattr(mod.sys, "stderr", stderr)
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_FRAMES", "2")
    assert list(mod.track(range(5), desc="unknown")) == list(range(5))
    assert "unknown: 2 frames" in stderr.getvalue()
    assert "unknown: 4 frames" in stderr.getvalue()
    stderr.seek(0)
    stderr.truncate(0)
    assert list(mod.track(range(5), total=5, enabled=False)) == list(range(5))
    assert stderr.getvalue() == ""


def test_known_total_does_not_print_frame_interval_updates(monkeypatch):
    stderr = _Stderr()
    monkeypatch.setattr(mod.sys, "stderr", stderr)
    monkeypatch.setenv("MDSCOPE_PROGRESS_EVERY_FRAMES", "25")
    monkeypatch.delenv("MDSCOPE_PROGRESS_EVERY_PERCENT", raising=False)
    list(mod.track(range(10000), total=10000, desc="rep"))
    lines = stderr.getvalue().splitlines()
    assert len(lines) == 101
    assert "(100/10000)" in lines[1]
    assert "(250/10000)" not in stderr.getvalue()
    assert "(10000/10000)" in lines[-1]



def test_nondefault_percent_cadence_is_monotonic_bounded_and_finishes(monkeypatch):
    import re
    import math
    stderr = _Stderr()
    monkeypatch.setattr(mod.sys, 'stderr', stderr)
    monkeypatch.setenv('MDSCOPE_PROGRESS_EVERY_PERCENT', '7')
    list(mod.track(range(123), total=123, desc='rep'))
    completed = [int(m.group(1)) for m in re.finditer(r'\((\d+)/123\)', stderr.getvalue())]
    assert completed[0] == 0 and completed[-1] == 123
    assert completed == sorted(set(completed))
    assert len(completed) <= math.ceil(100/7)+2


def test_stage_heartbeat_reports_during_work_and_labels_failure(monkeypatch):
    heartbeat = Event()
    messages = []

    def log(message, *, flush):
        assert flush
        messages.append(message)
        if "still running" in message:
            heartbeat.set()

    monkeypatch.setattr(mod, "print", log, raising=False)
    with pytest.raises(ValueError, match="failed input"):
        with mod.stage_progress("Statistics profiles", interval=.01):
            assert heartbeat.wait(timeout=5)
            raise ValueError("failed input")
    assert messages[0] == "Statistics profiles: started"
    assert any(msg.startswith("Statistics profiles: still running") for msg in messages)
    assert messages[-1].startswith("Statistics profiles: failed")
    with mod.stage_progress("disabled", enabled=False):
        pass
    assert not any("disabled" in msg for msg in messages)
