#!/usr/bin/env python3
"""Browse one synthetic MDscope output by condition; not a code diff."""
from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from html import escape
from pathlib import Path

CSS = """
body{font:15px/1.5 system-ui,sans-serif;max-width:1380px;margin:2rem auto;padding:0 1.4rem}
h1{font-size:1.6rem}h2{margin-top:2.3rem;font-size:1.3rem}
h3{font-size:1rem;margin:1.3rem 0 .5rem}
nav{display:flex;flex-wrap:wrap;gap:.9rem}
nav a{padding:.45rem .8rem;border:1px solid #ccc;border-radius:7px;text-decoration:none}
.notice{border:2px solid #555;padding:.8rem;font-weight:700}
figure{margin:0;max-width:100%}figure img{max-width:100%;height:auto}
.panels{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(340px,100%),1fr));gap:1.2rem}
.panel{border:1px solid #ddd;border-radius:8px;padding:.7rem;min-width:0}
figcaption,.muted{color:#555;font-size:.84rem;overflow-wrap:anywhere}
table{border-collapse:collapse;width:100%;font-size:.84rem}
td,th{padding:.45rem;border-bottom:1px solid #ddd;text-align:left}
td{overflow-wrap:anywhere}a{color:#1b527e} .small{font-size:.86rem}
"""


def _figure_group(name: str, identifier=None) -> tuple[str, str, str]:
    """Group by *condition*, not a replica folder's apparent scope."""
    parts = Path(name).parts
    if identifier and identifier.startswith('sdc4:'):
        return identifier.split(':')[1], parts[-1], 'all replicates'
    if len(parts) > 2 and parts[0].startswith("condition_"):
        return parts[0], parts[-1], parts[1]
    if len(parts) == 5 and parts[:3] == ("project_summary", "figures", "sdc4"):
        return parts[3], parts[-1], "all replicates"
    return "project", parts[-1], "all conditions"


def build_preview_gallery(preview: Path, output: Path) -> int:
    preview, output = Path(preview).resolve(), Path(output).resolve()
    info = json.loads((preview / "preview_manifest.json").read_text())
    groups: dict[str, dict[str, list[tuple[str,str]]]] = defaultdict(lambda: defaultdict(list))
    for name in info["plots"]:
        group, plot, replica = _figure_group(name, info.get("plot_ids", {}).get(name))
        groups[group][plot].append((name, replica))
    if not groups:
        raise ValueError("No figures in preview manifest")
    output.parent.mkdir(parents=True, exist_ok=True)
    def link(path: Path) -> str:
        return escape(Path(os.path.relpath(path, output.parent)).as_posix(), quote=True)
    sections = []
    order = sorted(key for key in groups if key != "project")
    if "project" in groups:
        order.append("project")
    for group in order:
        label = "All conditions / comparisons" if group == "project" else group.replace("_", " ").title()
        cards = []
        for filename, variants in sorted(groups[group].items()):
            panels = []
            for name, replica in sorted(variants):
                src = link(preview / name)
                cap = (f'<figcaption>{escape(replica)} · {escape(name)}</figcaption>')
                panels.append(f'<figure class="panel"><a href="{src}"><img loading="lazy" '
                              f'src="{src}" alt="{escape(name, quote=True)}"></a>{cap}</figure>')
            cards.append(f'<h3>{escape(filename.replace("_", " ").removesuffix(".png"))}</h3>'
                         f'<div class="panels">{"".join(panels)}</div>')
        sections.append(f'<section id="{escape(group, quote=True)}"><h2>{escape(label)}</h2>'
                        f'{"".join(cards)}</section>')
    inventory = preview / "table_previews/table_inventory.csv"
    tables = []
    if inventory.is_file():
        with inventory.open(newline="") as handle:
            for row in csv.DictReader(handle):
                schema = link(preview / row["columns_csv"])
                sample = link(preview / row["sample_csv"])
                tables.append(
                    f'<tr><td>{escape(row["source_parquet"])}</td>'
                    f'<td>{escape(row["n_rows"])}</td>'
                    f'<td>{escape(row["n_columns"])}</td>'
                    f'<td><a href="{schema}">Columns and types</a></td>'
                    f'<td><a href="{sample}">First {escape(row["sample_rows"])} rows</a></td></tr>'
                )
    table_html = (
        '<section id="data"><h2>Readable saved-data examples</h2>'
        '<p class="small">Samples are not full exports. Keep original Parquet as the '
        'authoritative numerical data. CSVs show column names/types and the first five rows.</p>'
        f'<p><a href="{link(inventory)}">Full Parquet inventory (CSV)</a></p>'
        '<table><thead><tr><th>Saved Parquet</th><th>Rows</th><th>Columns</th>'
        '<th>Schema CSV</th><th>Sample CSV</th></tr></thead><tbody>'
        + "".join(tables) + '</tbody></table></section>'
        if inventory.exists() else ""
    )
    nav = " ".join(
        f'<a href="#{escape(key, quote=True)}">{escape("All conditions" if key == "project" else key)}</a>'
        for key in order
    )
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport"
content="width=device-width,initial-scale=1"><title>MDscope preview — conditions</title>
<style>{CSS}</style></head><body>
<h1>MDscope — condition-centered figure preview</h1>
<p class="notice">SYNTHETIC EXAMPLE — NOT MOLECULAR DYNAMICS RESULTS</p>
<p>This is a single preview of the current implementation, <strong>not a before/after comparison</strong>.
Replicates belong to their condition. Project-wide figures may include multiple conditions;
those are grouped separately below. Saved numeric data are kept per replica.</p>
<nav>{nav} <a href="#data">CSV table previews</a></nav>
<p class="small">For code-review comparisons, open the separate <code>review/index.html</code>
in this artifact. “Before” and “after” refer to source revisions, not simulated conditions.</p>
{"".join(sections)}
{table_html}
</body></html>
"""
    output.write_text(html)
    return sum(len(items) for group in groups.values() for items in group.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    opts = parser.parse_args()
    print(f"Indexed {build_preview_gallery(opts.preview, opts.output)} figures")


if __name__ == "__main__":
    main()
