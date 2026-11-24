#!/usr/bin/env python3
"""Build a local HTML gallery comparing preview PNGs from two code revisions."""
from __future__ import annotations

import argparse
import json
import os
from html import escape
from pathlib import Path

from PIL import Image, ImageChops


def build_gallery(base: Path, head: Path, output: Path, *, only_name: str | None = None,
                  revisions: str | None = None):
    base, head, output = Path(base).resolve(), Path(head).resolve(), Path(output).resolve()
    for label, directory in (("base", base), ("head", head)):
        if not (directory / "preview_manifest.json").is_file():
            raise FileNotFoundError(f"{label} missing preview_manifest.json: {directory}")
    base_manifest = json.loads((base / "preview_manifest.json").read_text())
    head_manifest = json.loads((head / "preview_manifest.json").read_text())
    base_paths, head_paths = base_manifest["plots"], head_manifest["plots"]
    all_paths = sorted(set(base_paths) | set(head_paths))
    extra_note = ""
    if only_name is not None:
        requested = [name for name in all_paths if Path(name).name == only_name]
        if len(requested) != 1:
            raise ValueError(f"Selection example {only_name!r} must match one figure")
        all_paths = requested
        previous_extra = len(base_paths) - sum(name in base_paths for name in requested)
        updated_extra = len(head_paths) - sum(name in head_paths for name in requested)
        extra_note = (f'<p>One-figure selection test: requested <strong>{escape(only_name)}</strong>. '
                      f'Previous code generated {previous_extra} other figure(s); '
                      f'updated code generated {updated_extra} other figure(s). '
                      'The unrequested pictures are not repeated in this gallery.</p>')
    if not all_paths:
        raise ValueError("No preview plots to compare")
    output.parent.mkdir(parents=True, exist_ok=True)
    image_dir = output.parent / "diff"
    image_dir.mkdir(parents=True, exist_ok=True)

    def relative(path):
        return escape(Path(os.path.relpath(path, output.parent)).as_posix(), quote=True)

    rows = []
    for index, name in enumerate(all_paths):
        left, right = base / name, head / name
        panels = []
        for label, path in (("Previous code", left), ("Updated code", right)):
            if path.is_file():
                url = relative(path)
                panels.append(f'<div class="panel"><strong>{label}</strong>'
                              f'<a href="{url}"><img loading="lazy" src="{url}" '
                              f'alt="{escape(label + " " + name, quote=True)}"></a></div>')
            else:
                panels.append(f'<div class="panel"><strong>{label}</strong>'
                              '<p>Not generated</p></div>')
        diff = ""
        status = ("Added in head" if right.is_file() and not left.is_file() else
                  "Removed in head" if left.is_file() and not right.is_file() else
                  "Both revisions")
        plot_id = (head_manifest.get("plot_ids", {}).get(name) or
                   base_manifest.get("plot_ids", {}).get(name))
        if left.is_file() and right.is_file():
            with Image.open(left) as a, Image.open(right) as b:
                if a.size == b.size:
                    delta = ImageChops.difference(a.convert("RGB"), b.convert("RGB"))
                    if delta.getbbox():
                        path = image_dir / f"figure-{index:03}.png"
                        delta.save(path)
                        diff = (f'<div class="panel"><strong>Absolute pixel difference</strong>'
                                f'<a href="{relative(path)}"><img loading="lazy" '
                                f'src="{relative(path)}" alt="Pixel differences"></a></div>')
                    else:
                        diff = '<div class="panel"><strong>Unchanged pixels</strong></div>'
                else:
                    diff = '<div class="panel"><strong>Different image dimensions</strong></div>'
        metadata = f'<p>{escape(status)}' + (f' · Plot ID: {escape(plot_id)}' if plot_id else '') + '</p>'
        rows.append(f'<section><h2>{escape(name)}</h2>{metadata}<div class="row">'
                    f'{"".join(panels)}{diff}</div></section>')

    if only_name is None:
        reports = []
        for label, directory, manifest in (("Previous code", base, base_manifest), ("Updated code", head, head_manifest)):
            for name in manifest.get('reports', []):
                path = directory / name
                if path.is_file():
                    reports.append(f'<li>{label}: <a href="{relative(path)}">{escape(name)}</a></li>')
        if reports:
            rows.insert(0, '<section><h2>Saved statistical reports (PDF)</h2><ul>' + ''.join(reports) + '</ul></section>')

    revision_note = (f'<strong>Compared: {escape(revisions)}.</strong>' if revisions else '')
    html = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>MDscope figure review</title>
<style>
body{font:15px system-ui,sans-serif;margin:2rem auto;padding:0 1rem;max-width:1700px}
h1,h2{font-weight:650}h2{font-size:1rem;overflow-wrap:anywhere}
.notice{border:2px solid #555;padding:12px;font-weight:700}
section{border-top:1px solid #bbb;padding:18px 0}.row{display:flex;gap:15px;align-items:start}
.panel{flex:1;min-width:0}.panel img{width:100%;height:auto;border:1px solid #ddd}
.panel strong{display:block;margin-bottom:6px}
@media(max-width:750px){.row{flex-direction:column}.panel{width:100%}}
</style></head><body><h1>MDscope developer review — previous vs updated code</h1>
<p class="notice">SYNTHETIC EXAMPLES — NOT MOLECULAR DYNAMICS RESULTS.</p>
<p><strong>This compares source-code versions, not condition A against condition B or replica 1 against replica 2.</strong>
"Previous code" is the comparison base revision; "Updated code" is the proposed revision.
''' + revision_note + '''
Both receive the same deterministic synthetic Parquet inputs. When no rendering code
changed, the images will intentionally look identical. Pixel differences are diagnostic,
not scientific validation. Inspect axes, labels, units and legends.</p>
<p>For the nonduplicated condition-centered preview, open <code>../index.html</code>.</p>
''' + extra_note + "\n".join(rows) + "\n</body></html>\n"
    output.write_text(html)
    return len(all_paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--head", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--only-name", help="For selection reviews, show only the requested figure")
    parser.add_argument("--revisions", help="Text naming the compared base and head commits")
    args = parser.parse_args()
    number = build_gallery(args.base, args.head, args.output, only_name=args.only_name,
                           revisions=args.revisions)
    print(f"Wrote {number} figure comparison panels: {args.output}")


if __name__ == "__main__":
    main()
