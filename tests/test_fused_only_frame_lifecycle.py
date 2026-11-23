"""Fused-only scheduler must decode each frame once for all observers."""
from types import SimpleNamespace

from mdscope.pipeline import fused_all


class _Frames:
    closed = False

    def __getitem__(self, selection):
        return range(7)

    def close(self):
        self.closed = True


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


def test_pending_instances_use_one_pass_without_initializing_completed_runs(tmp_path, monkeypatch):
    """Partial recovery retains per-analysis and named-instance selection."""
    from mdscope.config import DEFAULTS
    from mdscope.registry import fused_output

    opened, selected, updates = [], {}, []
    def universe(*args):
        opened.append(args)
        return _Universe(*args)

    class Observer:
        def __init__(self, u, cfg, root, **options):
            selected[type(self).__name__] = options
        def prune_from(self, part):
            assert part == 0
        def update(self, ts, frame):
            updates.append((type(self).__name__, frame))
        def flush(self, part):
            pass
        def finalize(self, *args):
            return {}

    class Membrane(Observer):
        pass
    class Extensions(Observer):
        pass
    class Distances(Observer):
        pass
    class Mechanism(Observer):
        pass

    cfg = {
        "trajectory": {"start": None, "stop": None, "stride": 1, "chunk_frames": 3},
        "progress": {"enabled": False},
        "analyses": {
            **DEFAULTS["analyses"],
            "sasa": {"enabled": True, "runs": {
                "baseline": {"selection": "protein"},
                "environment": {"selection": "protein", "reference_run": "baseline"},
            }},
            "residue_membrane_distance": {"enabled": True, "runs": {
                "near": {"membrane_selection": "name P"},
                "far": {"membrane_selection": "name P"},
            }},
        },
    }
    monkeypatch.setattr(fused_all.mda, "Universe", universe)
    monkeypatch.setattr(fused_all, "membrane_enabled",
                        lambda cfg: ("bilayer_normal", "lipid_enrichment"))
    monkeypatch.setattr(fused_all, "extension_enabled",
                        lambda cfg: ("sasa", "hydrogen_bonds"))
    monkeypatch.setattr(fused_all, "distance_enabled", lambda cfg: ("near", "far"))
    monkeypatch.setattr(fused_all, "sdc4_enabled", lambda cfg: True)
    monkeypatch.setattr(fused_all, "FusedMembrane", Membrane)
    monkeypatch.setattr(fused_all, "FusedExtensions", Extensions)
    monkeypatch.setattr(fused_all, "FusedDistances", Distances)
    monkeypatch.setattr(fused_all, "SDC4Mechanisms", Mechanism)

    rep = SimpleNamespace(topology="top.psf", trajectory="frames.xtc",
                          system="test", label="rep1")
    only = (
        fused_output(tmp_path, "bilayer_normal"),
        fused_output(tmp_path, "sasa", run_id="environment"),
        fused_output(tmp_path, "residue_membrane_distance", run_id="near"),
    )
    fused_all.run_fused_only(rep, cfg, tmp_path, only=only)
    assert len(opened) == 1
    assert set(selected) == {"Membrane", "Extensions", "Distances"}
    assert selected["Membrane"]["names"] == ("bilayer_normal",)
    assert selected["Extensions"]["names"] == ("sasa",)
    assert selected["Extensions"]["sasa_run_ids"] == ("environment",)
    assert selected["Distances"]["run_ids"] == ("near",)
    assert all(len([x for x in updates if x[0] == cls]) == 7
               for cls in selected)
    # Scratch chunk directories are always cleaned after finalization.
    assert not list(tmp_path.glob(".fused-only-*/"))


def test_only_lipid_footprint_is_a_valid_observer_only_recovery(tmp_path, monkeypatch):
    """A sole missing footprint must not be skipped and must not run normals."""
    from mdscope.registry import fused_output

    selected, opened = [], []
    class Membrane:
        def __init__(self, u, cfg, root, **options):
            selected.append(options["names"])
        def prune_from(self, part): pass
        def update(self, ts, frame): pass
        def flush(self, part): pass
        def finalize(self): return {}

    made = []
    monkeypatch.setattr(fused_all.mda, "Universe",
                        lambda *args: (opened.append(args) or made.append(_Universe(*args)) or made[-1]))
    monkeypatch.setattr(fused_all, "membrane_enabled",
                        lambda cfg: ("bilayer_normal", "lipid_enrichment"))
    monkeypatch.setattr(fused_all, "extension_enabled", lambda cfg: ())
    monkeypatch.setattr(fused_all, "distance_enabled", lambda cfg: ())
    monkeypatch.setattr(fused_all, "sdc4_enabled", lambda cfg: False)
    monkeypatch.setattr(fused_all, "FusedMembrane", Membrane)
    rep = SimpleNamespace(topology="t", trajectory="x", system="s", label="r")
    cfg = {"trajectory": {"start": None, "stop": None, "stride": 1,
                          "chunk_frames": 2}, "progress": {"enabled": False}}
    fused_all.run_fused_only(
        rep, cfg, tmp_path, only=(fused_output(tmp_path, "lipid_enrichment"),))
    assert selected == [("lipid_enrichment",)]
    assert len(opened) == 1
    assert made[0].trajectory.closed  # released, not left for the GC
