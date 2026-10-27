"""Progress reporting for long MD trajectories (optional tqdm, safe batch fallback)."""
from __future__ import annotations
import sys
import time
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
    if tqdm is not None:
        with tqdm(iterable, total=total, desc=desc, unit="frame",
                  disable=not sys.stderr.isatty()) as bar:
            yield bar
        if not sys.stderr.isatty():
            print(f"{desc}: completed", file=sys.stderr, flush=True)
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
