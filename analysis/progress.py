"""Progress reporting for MD trajectories (interactive tqdm or batch milestones).

On Slurm, emit newline-delimited, flushed stderr progress with elapsed time,
processing speed and estimated time remaining. No additional trajectory reads.
Override the percentage cadence with MDSCOPE_PROGRESS_EVERY_PERCENT or the
unknown-length frame cadence with MDSCOPE_PROGRESS_EVERY_FRAMES (positive integers).
"""
from __future__ import annotations

import os
import sys
from contextlib import contextmanager
from time import perf_counter


def _cadence(name: str, default: int) -> int:
    value = int(os.environ.get(name, str(default)))
    if value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


@contextmanager
def progress(iterable, *, total=None, desc="", enabled=True):
    """Yield a wrapped iterable without changing trajectory access or memory use.

    Interactive terminals use tqdm with 1% minimum visual updates. Batch jobs
    with known totals report every 1% (default), while unknown-length jobs
    report every 25 completed frames. ETA uses measured mean throughput.
    """
    if not enabled:
        yield iterable
        return
    try:
        from tqdm.auto import tqdm
    except ImportError:
        tqdm = None
    if tqdm is not None and sys.stderr.isatty():
        with tqdm(iterable, total=total, desc=desc, unit="frame", mininterval=0.0,
                  miniters=max(1, (total + 99) // 100) if total else 1) as bar:
            yield bar
        return

    every_frames = _cadence("MDSCOPE_PROGRESS_EVERY_FRAMES", 25) if not total else None
    every_percent = _cadence("MDSCOPE_PROGRESS_EVERY_PERCENT", 1) if total else None

    class Milestones:
        def __iter__(self):
            started = perf_counter()
            last_percent_bucket = 0
            if total is not None and total > 0:
                print(f"{desc}: 0% (0/{total}) | elapsed 00:00:00 | ETA calculating",
                      file=sys.stderr, flush=True)
            else:
                print(f"{desc}: started | elapsed 00:00:00",
                      file=sys.stderr, flush=True)
            for count, frame in enumerate(iterable, 1):
                # Yield before measuring/reporting so the frame's analysis
                # computation is included in elapsed time and ETA.
                yield frame
                bucket = ((100 * count // total) // every_percent
                          if total is not None and total > 0 else 0)
                report = (
                    (total is None or total <= 0) and count % every_frames == 0
                    or (total is not None and total > 0
                        and (bucket > last_percent_bucket or count == total))
                )
                if not report:
                    continue
                last_percent_bucket = bucket
                elapsed = max(0., perf_counter() - started)
                rate = count / elapsed if elapsed > 0 else 0.
                speed = f"{rate:.2f} frames/s" if rate > 0 else "calculating"
                if total is not None and total > 0:
                    remaining = max(0, total - count)
                    eta = _duration(remaining / rate) if rate > 0 else "calculating"
                    print(
                        f"{desc}: {100 * count / total:.1f}% ({count}/{total})"
                        f" | elapsed {_duration(elapsed)}"
                        f" | {speed} | ETA {eta}",
                        file=sys.stderr, flush=True,
                    )
                else:
                    print(
                        f"{desc}: {count} frames | elapsed {_duration(elapsed)}"
                        f" | {speed}",
                        file=sys.stderr, flush=True,
                    )

    yield Milestones()


def track(iterable, *, total=None, desc="", enabled=True):
    """Generator wrapper for existing loops without changing body indentation."""
    with progress(iterable, total=total, desc=desc, enabled=enabled) as wrapped:
        yield from wrapped
