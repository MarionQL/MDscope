"""Lightweight opt-in timing of I/O and analysis phases in long SLURM jobs.

Enable by setting MDSCOPE_PROFILE=1. Timings are best-effort wall durations
including underlying Python/NumPy calls and can be compared across replicates.
"""
from __future__ import annotations

import os
from collections import defaultdict
from time import perf_counter

from ..config import atomic_json


class PhaseProfiler:
    def __init__(self, enabled=False):
        self.enabled = bool(enabled)
        self.seconds = defaultdict(float)
        self.calls = defaultdict(int)

    @classmethod
    def from_environment(cls):
        value = os.environ.get("MDSCOPE_PROFILE", "0").strip().lower()
        return cls(value not in {"", "0", "false", "off", "no"})

    def add(self, name, seconds):
        self.seconds[name] += float(seconds)
        self.calls[name] += 1

    def call(self, name, function, *args, **kwargs):
        if not self.enabled:
            return function(*args, **kwargs)
        started = perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            self.add(name, perf_counter() - started)

    def frames(self, iterable):
        """Measure trajectory decoding/iterator next() without preloading frames."""
        if not self.enabled:
            yield from iterable
            return
        iterator = iter(iterable)
        while True:
            started = perf_counter()
            try:
                frame = next(iterator)
            except StopIteration:
                return
            self.add("trajectory_decode", perf_counter() - started)
            yield frame

    def write(self, path, *, frames):
        if not self.enabled:
            return
        total = sum(self.seconds.values())
        rows = [
            {
                "stage": name,
                "seconds": float(seconds),
                "calls": self.calls[name],
                "seconds_per_frame": float(seconds / frames) if frames else None,
                "fraction_of_measured_time": float(seconds / total) if total else None,
            }
            for name, seconds in sorted(
                self.seconds.items(), key=lambda item: -item[1]
            )
        ]
        atomic_json(path, {
            "frames_completed": int(frames),
            "measured_seconds": float(total),
            "timings": rows,
            "note": "Instrumentation measures selected phases, not entire wall time; "
                    "unmeasured processing and startup/finalization are excluded.",
        })
        for row in rows:
            print(
                f"PROFILE {row['stage']}: {row['seconds']:.3f}s "
                f"({row['seconds_per_frame']:.6f}s/frame, "
                f"{row['calls']} calls)",
                flush=True,
            )
