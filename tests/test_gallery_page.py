"""The committed gallery page is the page the committed table produces.

``layout/build_gallery.py --no-render`` rebuilds ``layout/examples/index.html``
from ``layout/examples/gallery.tsv`` without laying anything out, so this runs
in a clean clone and needs neither the models nor Graphviz.
"""
from __future__ import annotations

import importlib.util


def _module(path, name):
    """Load a script by path: ``examples/`` and ``layout/`` both hold a build_gallery.py."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_index_html_matches_gallery_tsv(root):
    BG = _module(root / "layout" / "build_gallery.py", "layout_build_gallery")

    page = root / "layout" / "examples" / "index.html"
    before = page.read_bytes()
    try:
        BG.main(["--no-render"])
        after = page.read_bytes()
    finally:
        if page.read_bytes() != before:
            page.write_bytes(before)
    assert after.replace(b"\r\n", b"\n") == before.replace(b"\r\n", b"\n"), (
        "layout/examples/index.html is not what build_gallery.py --no-render writes "
        "from layout/examples/gallery.tsv"
    )
