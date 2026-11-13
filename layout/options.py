"""The layout options: one algorithm per family of ../../layout-benchmark-latest's registry.

Every option is a function ``(items, sizes, edges, seed) -> {item: (x, y)}``, run on an arbitrary graph
of boxes and nodes. The same function is used at every level of a grouping: on the quotient graph of
the top-level groups, on the groups inside a group, and on the species and reactions inside a leaf
group. ``sizes`` are box sizes in SVG points; engines that honour node sizes (Graphviz) get them,
the others are fixed up afterwards by the overlap removal in ``placement.py``.

Names and families follow ``layout_benchmark.layouts.REGISTRY``; this module re-implements the eight
picked options on top of the same libraries (igraph, Graphviz) because the benchmark's own runners
take a whole SBML ``Pathway``, not a sub-graph of one group.

Not runnable on this machine, hence not offered: every cytoscape.js / ELK / graphology / Reagraph
option (they need Node), and the pathway-specific tools (``rules-claude``, ``ecmapper``, ``metacarto``,
``sbml2map``, ``sbmlnetwork``) which need their own checkouts under ``.tools/``.

Run:
    python layout/options.py            # list the options, and whether each one can run here
"""
from __future__ import annotations

import math
import pathlib
import shutil
import subprocess
import tempfile
import threading

PT_PER_INCH = 72.0
GV_TIMEOUT = 300  # seconds, as in layout-benchmark
_IGRAPH_LOCK = threading.Lock()


# --------------------------------------------------------------------------------------- Graphviz

def _graphviz(engine: str, extra: dict[str, str]):
    """Run one Graphviz engine over the item graph, with the real box sizes (``-Tplain``, inches)."""

    def run(items, sizes, edges, seed):
        if not items:
            return {}
        if len(items) == 1:
            return {items[0]: (0.0, 0.0)}
        ids = {it: f"n{i}" for i, it in enumerate(items)}
        attrs = dict(extra)
        if engine in ("neato", "sfdp", "fdp"):
            attrs["start"] = f"random{seed}"
        lines = ["digraph G {", "  graph [" + ",".join(f'{k}="{v}"' for k, v in sorted(attrs.items())) + "];",
                 "  node [shape=box,fixedsize=true];"]
        for it in items:
            w, h = sizes[it]
            lines.append(f'  {ids[it]} [width={w / PT_PER_INCH:.4f},height={h / PT_PER_INCH:.4f},label=""];')
        for u, v, wt in edges:
            lines.append(f"  {ids[u]} -> {ids[v]} [weight={wt:g}];")
        lines.append("}")
        with tempfile.TemporaryDirectory() as td:
            f = pathlib.Path(td) / "g.gv"
            f.write_text("\n".join(lines), encoding="utf-8")
            proc = subprocess.run([engine, "-Tplain", str(f)], capture_output=True, text=True, timeout=GV_TIMEOUT)
        if proc.returncode != 0:
            raise RuntimeError(f"{engine}: {proc.stderr.strip()[:300]}")
        back = {v: k for k, v in ids.items()}
        pos = {}
        for line in proc.stdout.splitlines():
            part = line.split()
            if part and part[0] == "node" and part[1] in back:
                pos[back[part[1]]] = (float(part[2]) * PT_PER_INCH, float(part[3]) * PT_PER_INCH)
        missing = [it for it in items if it not in pos]
        if missing:
            raise RuntimeError(f"{engine}: {len(missing)} items unplaced")
        return pos

    return run


# ----------------------------------------------------------------------------------------- igraph

def _igraph(method: str, weighted: bool = False, **kwargs):
    """Run one igraph layout. igraph knows no node sizes; placement.py separates the boxes after."""

    def run(items, sizes, edges, seed):
        import random

        import igraph as ig

        if not items:
            return {}
        if len(items) == 1:
            return {items[0]: (0.0, 0.0)}
        index = {it: i for i, it in enumerate(items)}
        # directed: only Sugiyama uses the direction (it layers along it), the others ignore it
        g = ig.Graph(n=len(items), edges=[(index[u], index[v]) for u, v, _ in edges], directed=True)
        if edges:
            g.es["weight"] = [float(w) for _, _, w in edges]
        kw = dict(kwargs)
        if weighted and edges:
            kw["weights"] = "weight"
        # igraph's random number generator is global: one layout at a time, so a seed still fixes it
        with _IGRAPH_LOCK:
            ig.set_random_number_generator(random.Random(seed))
            try:
                coords = getattr(g, method)(**kw)
            finally:
                ig.set_random_number_generator(random)
        # a node cell is 74 pt wide in the grid renderer: igraph lays out in unit space
        return {it: (float(coords[i][0]) * 74.0, float(coords[i][1]) * 74.0) for it, i in index.items()}

    return run


# ---------------------------------------------------------------------------------------- options

#: name -> (function, family in layout-benchmark's registry, documentation link)
OPTIONS = {
    "graphviz-dot": (_graphviz("dot", {"rankdir": "TB", "ranksep": "0.45", "nodesep": "0.3"}),
                     "layered", "https://graphviz.org/docs/layouts/dot/"),
    "igraph-sugiyama": (_igraph("layout_sugiyama", hgap=1.4, vgap=1.4),
                        "layered (stands in for tree)",
                        "https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_sugiyama"),
    "graphviz-neato": (_graphviz("neato", {"overlap": "prism", "sep": "+12", "pack": "true"}),
                       "stress", "https://graphviz.org/docs/layouts/neato/"),
    "graphviz-sfdp": (_graphviz("sfdp", {"overlap": "prism", "sep": "+12", "pack": "true"}),
                      "force-directed", "https://graphviz.org/docs/layouts/sfdp/"),
    "igraph-drl": (_igraph("layout_drl", weighted=True), "for large graphs",
                   "https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_drl"),
    "graphviz-twopi": (_graphviz("twopi", {"overlap": "prism", "sep": "+12", "pack": "true"}),
                       "radial", "https://graphviz.org/docs/layouts/twopi/"),
    "graphviz-circo": (_graphviz("circo", {"overlap": "prism", "sep": "+12", "pack": "true"}),
                       "circular", "https://graphviz.org/docs/layouts/circo/"),
    "igraph-mds": (_igraph("layout_mds"), "baseline",
                   "https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_mds"),
}

#: registry families that cannot run on a machine without Node or the tool checkouts
UNAVAILABLE = {
    "pathway-specific": "rules-claude, ecmapper, metacarto, sbml2map - need their checkouts under .tools/",
    "SBML-specific": "sbmlnetwork - needs its checkout under .tools/",
    "compartments as groups": "fcose / cose-bilkent / reagraph compound layouts - need Node",
    "compound-capable force-directed": "fcose, cose-bilkent - need Node",
    "stress, with constraints": "cola - needs Node",
    "tree": "elk-mrtree - needs Node; igraph-sugiyama stands in as the hierarchical option",
}


def available(name: str) -> bool:
    """Whether the tool this option needs is on this machine."""
    if name.startswith("graphviz-"):
        return shutil.which(name.split("-", 1)[1]) is not None
    try:
        import igraph  # noqa: F401

        return True
    except ImportError:
        return False


def run(name: str, items, sizes, edges, seed: int = 0):
    """Lay out one item graph with one option. ``edges`` is ``[(u, v, weight), ...]``, undirected."""
    if name not in OPTIONS:
        raise ValueError(f"unknown option {name!r}, one of: {', '.join(OPTIONS)}")
    pos = OPTIONS[name][0](list(items), sizes, list(edges), seed)
    for it, (x, y) in pos.items():
        if not (math.isfinite(x) and math.isfinite(y)):
            raise RuntimeError(f"{name}: non-finite position for {it}")
    return pos


def main():
    print(f"{'option':20s}{'family':30s}{'runs here':10s}link")
    for name, (_f, family, link) in OPTIONS.items():
        print(f"{name:20s}{family:30s}{'yes' if available(name) else 'NO':10s}{link}")
    print("\nnot offered:")
    for family, why in UNAVAILABLE.items():
        print(f"  {family:34s}{why}")


if __name__ == "__main__":
    main()
