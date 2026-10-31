"""Fused-only scheduler must decode each frame once for all observers."""
from types import SimpleNamespace
from mdscope.analysis.pipeline import fused_all


class _Frames:
    def __getitem__(self, selection):
        return range(7)


class _Universe:
    def __init__(self, *args):
        self.trajectory = _Frames()


def test_fused_only_updates_two_families_in_one_pass(tmp_path, monkeypatch):
    events = []
    opened = []

    def universe(*args):
        opened.append(args)
        return _Universe(*args)

    class Observer:
        def __init__(self, *args):
            self.name = type(self).__name__
        def prune_from(self, part):
            events.append((self.name, "prune", part))
        def update(self, ts, frame):
            events.append((self.name, "frame", frame))
        def flush(self, part):
            events.append((self.name, "flush", part))
        def finalize(self):
            return {}

    class Membrane(Observer):
        pass

    class Extensions(Observer):
        pass

    monkeypatch.setattr(fused_all.mda, "Universe", universe)
    monkeypatch.setattr(fused_all, "membrane_enabled", lambda cfg: ("membrane_order",))
    monkeypatch.setattr(fused_all, "extension_enabled", lambda cfg: ("hydrogen_bonds",))
    monkeypatch.setattr(fused_all, "distance_enabled", lambda cfg: ())
    monkeypatch.setattr(fused_all, "sdc4_enabled", lambda cfg: False)
    monkeypatch.setattr(fused_all, "FusedMembrane", Membrane)
    monkeypatch.setattr(fused_all, "FusedExtensions", Extensions)

    rep = SimpleNamespace(topology="top.psf", trajectory="frames.xtc", system="test", label="rep1")
    cfg = {"trajectory": {"start": None, "stop": None, "stride": 1, "chunk_frames": 3},
           "progress": {"enabled": False}}
    fused_all.run_fused_only(rep, cfg, tmp_path)
    assert len(opened) == 1
    for cls in ("Membrane", "Extensions"):
        assert [event[2] for event in events if event[:2] == (cls, "frame")] == list(range(7))
        assert [event[2] for event in events if event[:2] == (cls, "flush")] == [0, 1, 2]
