"""Assemble the GitHub Pages site from the galleries this repository already carries.

The site is the two rendered galleries, unchanged, behind one landing page - the same shape as
layout-benchmark's site (an orphan ``gh-pages`` branch with ``.nojekyll``, a small ``index.html``
at the root, and one folder per report):

    index.html          the landing page
    .nojekyll           so GitHub Pages serves the folders as they are, without Jekyll
    grouping/           from examples/        - the grouping gallery, 52 pictures
    layout/             from layout/examples/ - the layout gallery, 288 pictures

Nothing is re-rendered: the pages and their SVGs are copied as committed, so the site always shows
exactly what the repository holds. Rebuild the galleries first (README, steps 5 and 6) to change it.

Run:
    python scripts/build_site.py                 # writes site/
    python scripts/build_site.py -o build/site   # somewhere else
"""
from __future__ import annotations

import argparse
import html
import pathlib
import shutil

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = "https://github.com/elixir2026-project22/layout-benchmark-big-maps"

#: (folder on the site, folder in the repository, title, one line about it)
REPORTS = [
    ("grouping", "examples", "grouping",
     "52 pictures: every model under every grouping strategy, drawn by one renderer"),
    ("layout", "layout/examples", "layout",
     "288 pictures: 8 layout options x 9 groupings x 4 big models"),
]


def copy_report(source: pathlib.Path, target: pathlib.Path) -> tuple[int, int]:
    """Copy one gallery (its page, its table and its SVGs). Returns (files, bytes)."""
    target.mkdir(parents=True, exist_ok=True)
    files = size = 0
    for name in ("index.html", "gallery.tsv"):
        if (source / name).exists():
            shutil.copy2(source / name, target / name)
            files += 1
            size += (target / name).stat().st_size
    svg_out = target / "svg"
    if svg_out.exists():
        shutil.rmtree(svg_out)
    svg_out.mkdir()
    for svg in sorted((source / "svg").glob("*.svg")):
        shutil.copy2(svg, svg_out / svg.name)
        files += 1
        size += svg.stat().st_size
    return files, size


def landing(entries) -> str:
    """The landing page, in the same minimal style as layout-benchmark's site."""
    items = "".join(
        f'<li><a href="{folder}/index.html">{html.escape(title)}</a> - {html.escape(about)} '
        f'<span class=n>({files} files, {size / 1e6:.0f} MB)</span></li>'
        for folder, title, about, files, size in entries
    )
    return (
        "<!doctype html><meta charset=utf-8><title>Layout benchmark - big maps</title>"
        "<style>body{font-family:system-ui,sans-serif;margin:24px;max-width:56em;line-height:1.5}"
        "h1{font-size:22px}.n{color:#666}li{margin:.4em 0}</style>"
        "<h1>Layout benchmark - big maps</h1>"
        "<p>Grouping big pathway maps, and distributing the elements inside each group and the groups "
        "themselves. BioHackathon Europe 2026, project 22 (Open Pathway Maps). "
        f'<a href="{REPO}">Code</a>.</p>'
        f"<ul>{items}</ul>"
        "<p class=n>Each cell links the full-size SVG. The pictures are the repository's own "
        "<code>examples/</code> and <code>layout/examples/</code>, copied unchanged.</p>"
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--out", default="site")
    a = ap.parse_args(argv)
    site = pathlib.Path(a.out)
    if not site.is_absolute():
        site = ROOT / site
    site.mkdir(parents=True, exist_ok=True)
    entries = []
    for folder, rel, title, about in REPORTS:
        source = ROOT / rel
        if not (source / "index.html").exists():
            raise SystemExit(f"{source}/index.html is missing - build the gallery first (README, step 5 or 6)")
        files, size = copy_report(source, site / folder)
        entries.append((folder, title, about, files, size))
        print(f"  {folder:10s} {files:4d} files  {size / 1e6:7.1f} MB  from {rel}/")
    (site / ".nojekyll").write_text("", encoding="utf-8")
    (site / "index.html").write_text(landing(entries), encoding="utf-8")
    total = sum(e[4] for e in entries)
    print(f"{site}: {sum(e[3] for e in entries) + 2} files, {total / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
