"""Stable plot identifiers and reporting-specific figure path presentation."""
from __future__ import annotations

from pathlib import Path

from ..contracts.plot_validation import (
    DYNAMIC_VARIANTS,
    SDC4_CONDITION_STEMS,
    STATIC_IDS,
    describe_plot,
)
from ..contracts.plot_validation import (
    validate_figure_options as validate_figure_options,
)


def figure_dpi(cfg, plot_id, default):
    if plot_id is None:
        return default
    return cfg.get("plots", {}).get("options", {}).get(plot_id, {}).get("dpi", default)



def requires_finalization(identifier: str) -> bool:
    """Derived descriptive figures need publication without requiring inference."""
    descriptor = describe_plot(identifier)
    return descriptor["requires_statistics"] or descriptor["family"] in {"diffusion", "shared_states"}


def figure_family(path):
    """Figure family of a saved figure path (one of FIGURE_FAMILIES)."""
    parts = Path(path).parts
    if Path(path).parent.name in {"primary", "secondary", "tertiary"}:
        identifier = id_from_figure_path(Path(path))
        if identifier:
            return identifier.split(":", 1)[0]
    for token, family in (('shared_states','shared_states'),('sdc4','sdc4'),('sdc4_mechanisms','sdc4'),
                          ('solvation','solvation'),('diffusion','diffusion'),
                          ('lipid_contacts','lipid_contacts'),('membrane','lipid_enrichment'),
                          ('interactions','interactions'),('profiles','profiles'),
                          ('time_series','time_series')):
        if token in parts:
            return family
    return 'other'


def id_from_figure_path(path: Path) -> str | None:
    """Recover stable IDs for the supported persisted output filenames."""
    path = Path(path)
    if path.suffix != ".png":
        return None
    parts = path.parts
    name = path.stem
    if len(parts) >= 3 and path.parent.name in {"primary", "secondary", "tertiary"}:
        from urllib.parse import unquote

        from ..contracts.plot_validation import validate_plot_ids
        identifier = ":".join(unquote(part) for part in name.split("@"))
        # Page identities are presentation companions of the same selected ID.
        if identifier.rsplit(":", 1)[-1].startswith("page_"):
            identifier = identifier.rsplit(":", 1)[0]
        try:
            validate_plot_ids([identifier])
        except ValueError:
            return None
        return identifier
    if (path.parent.name.startswith("residue_membrane_distance__")
            and name.endswith("__distance_heatmap")):
        system = name.removesuffix("__distance_heatmap")
        return f"profiles:{path.parent.name}:{system}:distance_heatmap"
    if "time_series" in parts or "profiles" in parts:
        family = "time_series" if "time_series" in parts else "profiles"
        if "__" in name:
            stem, variant = name.rsplit("__", 1)
            if variant in DYNAMIC_VARIANTS[family]:
                return f"{family}:{stem}:{variant}"
    if "interactions" in parts and len(parts) >= 2:
        kind = path.parent.name
        if name in {"occupancy_heatmap", "lifetime_distribution"}:
            return f"interactions:{kind}:{name}"
        if "__" in name:
            metric, variant = name.rsplit("__", 1)
            if variant in {"bootstrap_table", "comparison_table"}:
                return f"interactions:{kind}:{metric}:{variant}"
    if "shared_states" in parts:
        identifier = f"shared_states:{name}"
        return identifier if identifier in STATIC_IDS else None
    if "lipid_contacts" in parts:
        identifier = f"lipid_contacts:{name}"
        if name.startswith(("lipid_contact_matrix_", "lipid_contact_difference_",
                            "lipid_contact_class_matrix_", "lipid_contact_class_difference_")):
            return identifier
    if (len(parts) >= 5 and parts[-3] == "sdc4" and
            parts[-4] == "figures" and name in SDC4_CONDITION_STEMS):
        return f"sdc4:{parts[-2]}:{name}"
    for family in ("solvation", "diffusion", "lipid_enrichment"):
        identifier = f"{family}:{name}"
        if identifier in STATIC_IDS:
            return identifier
    return None



def figure_path(cfg, path, *, identifier=None):
    """Translate a renderer's legacy destination through the shared catalog."""
    from ..registry import figure_descriptor
    path = Path(path)
    identifier = identifier or id_from_figure_path(path)
    if identifier is None:
        raise ValueError(f"Unmapped figure destination: {path}; declare a stable plot ID")
    if 'project_summary' not in path.parts:
        raise ValueError(f"Figure destination must be below project_summary: {path}")
    index = path.parts.index('project_summary')
    root = Path(*path.parts[:index])
    destination = root / figure_descriptor(cfg, identifier)['destination']
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination


def apply_figure_options(fig, options):
    """User-authored display text, never replacements for units or estimands."""
    if 'title' in options:
        fig.suptitle(options['title'])
    annotations = options.get('annotations', [])
    if annotations:
        import textwrap
        text = '\n'.join(textwrap.fill(item, width=110) for item in annotations)
        fig.text(.01, .005, text, va='bottom', ha='left', fontsize=8)
        fig.subplots_adjust(bottom=max(.16, .03 + .025 * len(text.splitlines())) )


def save_figure(fig, cfg, path, *, identifier=None, **kwargs):
    """One destination and presentation owner for all disk PNG writers."""
    path = Path(path)
    identifier = identifier or id_from_figure_path(path)
    if path.parent.name not in {'primary', 'secondary', 'tertiary'}:
        path = figure_path(cfg, path, identifier=identifier)
    options = cfg.get('plots', {}).get('options', {}).get(identifier, {})
    apply_figure_options(fig, options)
    fig.savefig(path, **kwargs)
