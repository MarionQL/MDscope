#!/usr/bin/env python3
"""Build a local HTML gallery comparing preview PNGs from two code revisions."""
from __future__ import annotations

import argparse
from html import escape
import json
import os
from pathlib import Path

from PIL import Image, ImageChops


def build_gallery(base: Path, head: Path, output: Path):
    base, head, output = Path(base).resolve(), Path(head).resolve(), Path(output).resolve()
    for label, directory in (("base", base), ("head", head)):
        if not (directory / "preview_manifest.json").is_file():
            raise FileNotFoundError(f"{label} missing preview_manifest.json: {directory}")
    base_paths = json.loads((base / "preview_manifest.json").read_text())["plots"]
    head_paths = json.loads((head / "preview_manifest.json").read_text())["plots"]
    all_paths = sorted(set(base_paths) | set(head_paths))
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
        for label, path in (("Base", left), ("PR head", right)):
            if path.is_file():
                url = relative(path)
                panels.append(f'<div class="panel"><strong>{label}</strong>'
                              f'<a href="{url}"><img loading="lazy" src="{url}" '
                              f'alt="{escape(label + " " + name, quote=True)}"></a></div>')
            else:
                panels.append(f'<div class="panel"><strong>{label}</strong>'
                              '<p>Not generated</p></div>')
        diff = ""
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
        rows.append(f'<section><h2>{escape(name)}</h2><div class="row">'
                    f'{"".join(panels)}{diff}</div></section>')

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
</style></head><body><h1>MDscope figure review — base vs PR head</h1>
<p class="notice">SYNTHETIC EXAMPLES — NOT MOLECULAR DYNAMICS RESULTS.</p>
<p>Both revisions render the same deterministic fake Parquet inputs through
the MDscope plotting functions. Pixel differences are diagnostic, not
automatic scientific validation. Inspect labels, units, axes and legends.</p>
''' + "\n".join(rows) + "\n</body></html>\n"
    output.write_text(html)
    return len(all_paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--head", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    number = build_gallery(args.base, args.head, args.output)
    print(f"Wrote {number} figure comparison panels: {args.output}")


if __name__ == "__main__":
    main()
