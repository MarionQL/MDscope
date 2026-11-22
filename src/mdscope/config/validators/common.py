"""Shared selection and interval validators."""
import math
import re

def _validate_time_window(window, label):
    if window is not None and (
            not isinstance(window, (list, tuple)) or len(window) != 2
            or any(type(value) not in (int, float) or not math.isfinite(value) for value in window)
            or window[0] >= window[1]):
        raise ValueError(f"{label} requires two increasing finite times")


def _validate_tail_spec(tail_spec, label):
    if not isinstance(tail_spec, dict) or not (
        tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names")
    ):
        raise ValueError(f"{label} requires carbon_pattern or carbon_names")
    if tail_spec.get("carbon_pattern"):
        try:
            re.compile(str(tail_spec["carbon_pattern"]))
        except re.error as exc:
            raise ValueError(f"{label}.carbon_pattern is not a valid regular expression: {exc}") from exc
    names, positions = tail_spec.get("carbon_names"), tail_spec.get("positions")
    if names and positions is not None and len(names) != len(positions):
        raise ValueError(f"{label}.positions must match carbon_names")


def _validate_tail_naming(order):
    """Tail naming shared by membrane_order and SDC4 local order."""
    tails = order.get("tails", {}) or {}
    if not tails:
        raise ValueError("membrane_order.tails must define at least one tail")
    for tail, tail_spec in tails.items():
        _validate_tail_spec(tail_spec, f"membrane_order.tails.{tail}")
    for resname, overrides in (order.get("lipid_overrides") or {}).items():
        if not isinstance(overrides, dict):
            raise ValueError(f"membrane_order.lipid_overrides.{resname} must map tail names")
        for tail, tail_spec in overrides.items():
            label = f"membrane_order.lipid_overrides.{resname}.{tail}"
            if tail not in tails:
                raise ValueError(f"{label}: tail {tail!r} is not defined in membrane_order.tails")
            if not isinstance(tail_spec, dict):
                raise ValueError(f"{label} must be a mapping")
            if tail_spec.get("carbon_pattern") or tail_spec.get("carbon_names"):
                _validate_tail_spec(tail_spec, label)
