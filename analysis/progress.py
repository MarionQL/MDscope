"""Progress reporting for long MD trajectories (optional tqdm, safe batch fallback)."""
from __future__ import annotations
import sys
from contextlib import contextmanager

@contextmanager
def progress(iterable, *, total=None, desc="", enabled=True):
    """Yield the wrapped iterable; no additional trajectory seeks or memory use.

    tqdm is optional. On HPC batch nodes print milestones rather than carriage returns.
    """
    if not enabled:
        yield iterable
        return
    try:
        from tqdm.auto import tqdm
    except ImportError:
        tqdm = None
    if tqdm is not None and sys.stderr.isatty():
        with tqdm(iterable, total=total, desc=desc, unit="frame") as bar:
            yield bar
        return
    class Milestones:
        def __iter__(self):
            last = -1
            for count, frame in enumerate(iterable, 1):
                if total is not None and total>0:
                    percent = int(count*100/total)//10
                    if percent != last and percent<=10:
                        print(f"{desc}: {min(count*100//total,100)}% ({count}/{total})",
                              file=sys.stderr, flush=True)
                        last = percent
                elif count % 1000 == 0:
                    print(f"{desc}: {count} frames",file=sys.stderr,flush=True)
                yield frame
    yield Milestones()

def track(iterable, *, total=None, desc="", enabled=True):
    """Generator wrapper for long pre-existing loops without changing body indentation."""
    with progress(iterable,total=total,desc=desc,enabled=enabled) as wrapped:
        yield from wrapped
