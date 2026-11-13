"""Place a grouped model with one layout option, the same algorithm at every level.

The grouping gives a tree of boxes (``examples/render_grouping.py`` builds it: nested group boxes,
localisation levels drawn as a cell from ``compartments.tsv``, boundary nodes in a band). This module
replaces that renderer's fixed grid and shelf packing by a real layout:

1. a group with exactly one entity and no sub-group is not a group - its node is moved up into the
   parent and the empty box dropped (``collapse_singletons``);
2. inside every box, its direct children (sub-boxes, already sized) and its own species and reactions
   are one item graph - two items are joined when the model has an edge between any member of the one
   and any member of the other, weighted by how many such edges there are - and the layout option
   places that graph (``lay``);
3. the same call at the top level places the top-level groups themselves, so groups are distributed
   by the same logic as the elements inside them;
4. items are then pushed apart until no two boxes overlap (``separate``), which is what gives a box
   its size for the level above.

Boxes keep their meaning: a space is still drawn as a compartment, a band still sits on the boundary
of its space (the option lays out the band's contents, the band itself is pinned to the top border).

Run:
    python layout/placement.py panel/features/e_coli_core strategies/hybrid/out/e_coli_core \
        --option graphviz-dot -o layout/out/graphviz-dot/hybrid/e_coli_core
"""
from __future__ import annotations

import argparse
import collections
import csv
import math
import pathlib
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "examples"))
import options as OPT  # noqa: E402
import render_grouping as RG  # noqa: E402

# margin added around a node so that labels and edge ends have room
NODE_MARGIN_X, NODE_MARGIN_Y = 14.0, 14.0
PAD, HEAD, GAP = RG.PAD, RG.HEAD, RG.GAP
BOX_GAP = 16.0      # free space kept between two boxes of one level
SPACING = 0.9       # nearest neighbours this many node widths apart, before separation
FILL = 0.10         # the share of a level's bounding box the items cover, before separation


def ikey(item) -> str:
    """A stable name for an item of a level: a drawn node, or a child box. Never a memory address."""
    return f"n:{item[1]}" if item[0] == "n" else f"b:{item[1].key}"


def node_size(node_row) -> tuple[float, float]:
    if node_row["kind"] == "reaction":
        return RG.RX + NODE_MARGIN_X, RG.RX + NODE_MARGIN_Y
    return RG.SP_W + NODE_MARGIN_X, RG.SP_H + NODE_MARGIN_Y


# ------------------------------------------------------------------- the box tree after collapsing

def collapse_singletons(box, dropped):
    """A group of one entity is not a group: move its node up, drop the box. Spaces and bands stay."""
    for child in list(box.children):
        collapse_singletons(child, dropped)
    if box.band is not None:
        collapse_singletons(box.band, dropped)
    keep = []
    for child in box.children:
        single = child.kind == "group" and not child.children and child.band is None and len(child.nodes) == 1
        empty = child.kind == "group" and not child.children and child.band is None and not child.nodes
        if single:
            box.nodes.extend(child.nodes)
            dropped.append(child.key)
        elif empty:
            dropped.append(child.key)
        else:
            keep.append(child)
    box.children = keep


def walk(box, chain, out):
    box.chain = chain + [box]
    out.append(box)
    for child in ([box.band] if box.band is not None else []) + box.children:
        walk(child, box.chain, out)


def members(box, into=None):
    """Every drawn key under a box (own nodes, band, children)."""
    into = [] if into is None else into
    into.extend(box.nodes)
    for child in ([box.band] if box.band is not None else []) + box.children:
        members(child, into)
    return into


# ------------------------------------------------------------------------------ overlap separation

def separate(pos, sizes, keys, gap=BOX_GAP, pushes=200):
    """Push items apart until no two boxes overlap, then guarantee it with one uniform scale.

    Scaling the positions while the sizes stay fixed never brings two items closer, so the final
    scale resolves every remaining overlapping pair and can create none."""
    if len(keys) < 2:
        return pos
    p = np.array([pos[k] for k in keys], dtype=float)
    half = np.array([[sizes[k][0] / 2 + gap / 2, sizes[k][1] / 2 + gap / 2] for k in keys])
    # deterministic jitter so coincident centres can be pushed apart at all
    same = cKDTree(p).query_pairs(r=1e-9)
    for i, _j in sorted(same):
        p[i] += ((i % 7) - 3) * 1e-3, ((i % 5) - 2) * 1e-3

    def overlaps():
        radius = float(np.hypot(*(half.max(axis=0) * 2)))
        pairs = sorted(cKDTree(p).query_pairs(r=radius))
        out = []
        for i, j in pairs:
            d = np.abs(p[i] - p[j])
            need = half[i] + half[j]
            if d[0] < need[0] and d[1] < need[1]:
                out.append((i, j, d, need))
        return out

    for _ in range(pushes):
        bad = overlaps()
        if not bad:
            break
        shift = np.zeros_like(p)
        for i, j, d, need in bad:
            over = need - d
            axis = 0 if over[0] <= over[1] else 1
            step = over[axis] / 2 + 0.5
            direction = 1.0 if p[i][axis] >= p[j][axis] else -1.0
            shift[i][axis] += direction * step
            shift[j][axis] -= direction * step
        p += shift
    bad = overlaps()
    if bad:
        scale = 1.0
        for _i, _j, d, need in bad:
            with np.errstate(divide="ignore", invalid="ignore"):
                s = np.where(d > 1e-9, need / np.maximum(d, 1e-9), np.inf)
            scale = max(scale, float(min(s)))
        p *= min(scale, 50.0)
    return {k: (float(p[i][0]), float(p[i][1])) for i, k in enumerate(keys)}


def components(keys, edges):
    """The connected components of one level's item graph, in the order the items were given.

    Every option is run on one component at a time and the components are then packed, so that an
    algorithm that throws unconnected pieces far apart (DrL does, by design) is judged on how it
    draws a connected piece. Graphviz does the same for its own layouts with ``pack=true``."""
    parent = {k: k for k in keys}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for u, v, _w in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
    groups = collections.OrderedDict()
    for k in keys:
        groups.setdefault(find(k), []).append(k)
    return list(groups.values())


def shelf_pack(parts, gap=BOX_GAP, aspect=1.4):
    """Place component bounding boxes in rows, tallest first. ``parts`` is [(width, height)],
    the result is each one's top-left corner."""
    order = sorted(range(len(parts)), key=lambda i: (-parts[i][1], -parts[i][0], i))
    area = sum(w * h for w, h in parts)
    limit = max([math.sqrt(area * aspect)] + [w for w, _h in parts]) if parts else 0.0
    out = [(0.0, 0.0)] * len(parts)
    x = y = row_h = 0.0
    for i in order:
        w, h = parts[i]
        if x > 0 and x + w > limit + 1e-6:
            x, y, row_h = 0.0, y + row_h + gap, 0.0
        out[i] = (x, y)
        x += w + gap
        row_h = max(row_h, h)
    return out


def normalise(pos, sizes, keys):
    """Scale a layout to a readable density, before the boxes are pushed apart.

    Two scales are computed and the smaller one is used: one that puts nearest neighbours SPACING
    node widths apart, and one that fills FILL of the bounding box with item area. The second is
    what rescues a degenerate result - MDS, for instance, returns coincident points for
    structurally equivalent nodes, and then the nearest-neighbour distance is ~0 and its scale
    useless. A layout with no extent at all (every item in one point) is packed in rows instead."""
    if len(keys) < 2:
        return pos
    p = np.array([pos[k] for k in keys], dtype=float)
    size = np.array([sizes[k] for k in keys], dtype=float)
    nn = cKDTree(p).query(p, k=2)[0][:, 1]
    median = float(np.median(nn[nn > 0])) if (nn > 0).any() else 0.0
    by_nn = SPACING * float(np.median(size[:, 0])) / median if median > 1e-9 else math.inf
    span = p.max(axis=0) - p.min(axis=0)
    box_area = float(span[0] * span[1])
    item_area = float((size[:, 0] * size[:, 1]).sum())
    by_area = math.sqrt(item_area / (FILL * box_area)) if box_area > 1e-9 else math.inf
    scale = min(by_nn, by_area)
    if not math.isfinite(scale):  # no extent: fall back to packing the items in rows
        corners = shelf_pack([tuple(s) for s in size])
        return {k: (corners[i][0] + size[i][0] / 2, corners[i][1] + size[i][1] / 2)
                for i, k in enumerate(keys)}
    return {k: (pos[k][0] * scale, pos[k][1] * scale) for k in keys}


# ----------------------------------------------------------------------------------- the placement

class Placement:
    """One model, one grouping, one layout option."""

    def __init__(self, features_dir, out_dir, option, seed=0):
        self.option, self.seed = option, seed
        root, nodes, edges, comp, tree, home, clones, meta = RG.build(features_dir, out_dir)
        self.root, self.nodes, self.edges, self.comp, self.meta = root, nodes, edges, comp, meta
        self.clones = clones
        self.dropped_singletons = []
        collapse_singletons(root, self.dropped_singletons)
        self.boxes = []
        walk(root, [], self.boxes)
        self.boxes = [b for b in self.boxes if b.kind != "root"]

        self.home = {}  # drawn key -> its box
        for b in [root] + self.boxes:
            for k in b.nodes:
                self.home[k] = b
        self.item_of = {}  # box -> {drawn key -> item under that box}
        for b in [root] + self.boxes:
            table = {k: ("n", k) for k in b.nodes}
            for child in ([b.band] if b.band is not None else []) + b.children:
                for k in members(child):
                    table[k] = ("b", child)
            self.item_of[id(b)] = table
        self.drawn_edges = self._drawn_edges()
        self.pos = {}
        self._lay(root)
        self._position(root, 0.0, 0.0)

    # the copy of a species a reaction connects to: a clone in the reaction's own boxes, else its home
    def _drawn(self, species, reaction):
        if reaction not in self.home:
            return species
        chain = {id(x) for x in self.home[reaction].chain}
        best = None
        for b, key in self.clones.get(species, ()):
            if id(b) in chain and (best is None or b.depth_chain > best[0]):
                best = (b.depth_chain, key)
        return best[1] if best else species

    def _drawn_edges(self):
        for b in self.boxes:
            b.depth_chain = len(b.chain)
        out = []
        for e in self.edges:
            r, s = e["reaction_id"], self._drawn(e["species_id"], e["reaction_id"])
            if r in self.home and s in self.home:
                out.append((s, r) if e["role"] != "product" else (r, s))
        return out

    def _item_edges(self, box, allowed):
        """Edges of the item graph of one box: aggregated model edges between its direct items.

        A boundary band is not laid out with the rest (it is pinned to the border of its space), so
        edges that end in it are left out of this level."""
        table = self.item_of[id(box)]
        allowed = set(allowed)
        weight = collections.Counter()
        for u, v in self.drawn_edges:
            iu, iv = table.get(u), table.get(v)
            if iu in allowed and iv in allowed and iu != iv:
                weight[(iu, iv)] += 1
        return [(u, v, w) for (u, v), w in sorted(weight.items(), key=lambda kv: (ikey(kv[0][0]), ikey(kv[0][1])))]

    def _lay(self, box):
        """Bottom-up: size every child, then place children and own nodes with the layout option."""
        for child in ([box.band] if box.band is not None else []) + box.children:
            self._lay(child)
        box.children.sort(key=lambda c: c.key)
        items = [("n", k) for k in sorted(box.nodes)] + [("b", c) for c in box.children]
        sizes = {}
        for it in items:
            sizes[it] = node_size(self.nodes[it[1].split("#")[0]]) if it[0] == "n" else (it[1].w, it[1].h)
        if items:
            edges = self._item_edges(box, items)
            parts, spans = [], []
            for comp in components(items, edges):
                inside = set(comp)
                sub = [(u, v, w) for u, v, w in edges if u in inside and v in inside]
                pos = OPT.run(self.option, comp, sizes, sub, self.seed)
                pos = normalise(pos, sizes, comp)
                pos = separate(pos, sizes, comp)
                x0 = min(pos[i][0] - sizes[i][0] / 2 for i in comp)
                y0 = min(pos[i][1] - sizes[i][1] / 2 for i in comp)
                parts.append({i: (pos[i][0] - x0, pos[i][1] - y0) for i in comp})
                spans.append((max(pos[i][0] + sizes[i][0] / 2 for i in comp) - x0,
                              max(pos[i][1] + sizes[i][1] / 2 for i in comp) - y0))
            corners = shelf_pack(spans)
            box.ipos = {}
            for part, (cx, cy) in zip(parts, corners):
                for i, (dx, dy) in part.items():
                    box.ipos[i] = (cx + dx, cy + dy)
            box.ipos = {i: box.ipos[i] for i in items}
            content_w = max(cx + w for (cx, _cy), (w, _h) in zip(corners, spans))
            content_h = max(cy + h for (_cx, cy), (_w, h) in zip(corners, spans))
        else:
            content_w = content_h = 0.0
            box.ipos = {}
        box.isize = sizes
        band_h = 0.0
        if box.band is not None:
            box.band.w = max(box.band.w, content_w)
            band_h = box.band.h + GAP / 2
        head = HEAD if box.kind != "root" else 0.0
        box.content = (content_w, content_h)
        box.band_h = band_h
        box.w = max(content_w, box.band.w if box.band is not None else 0.0) + 2 * PAD
        box.h = content_h + band_h + 2 * PAD + head

    def _position(self, box, x, y):
        box.x, box.y = x, y
        head = HEAD if box.kind != "root" else 0.0
        top = y + PAD + head
        if box.band is not None:
            self._position(box.band, x + PAD, top)
            top += box.band_h
        for item, (dx, dy) in box.ipos.items():
            if item[0] == "n":
                self.pos[item[1]] = (x + PAD + dx, top + dy)
            else:
                child = item[1]
                self._position(child, x + PAD + dx - child.w / 2, top + dy - child.h / 2)

    # ------------------------------------------------------------------------------------- tables

    def positions_rows(self):
        rows = []
        for key, (x, y) in sorted(self.pos.items()):
            node = self.nodes[key.split("#")[0]]
            w, h = node_size(node)
            box = self.home[key]
            rows.append({"node_id": key.split("#")[0], "drawn_key": key, "clone": int("#" in key),
                         "kind": node["kind"], "box": box.key, "x": round(x, 2), "y": round(y, 2),
                         "width": round(w - NODE_MARGIN_X, 2), "height": round(h - NODE_MARGIN_Y, 2)})
        return rows

    def boxes_rows(self):
        rows = []
        for b in sorted(self.boxes, key=lambda b: (len(b.chain), b.key)):
            parent = b.chain[-2] if len(b.chain) > 1 else None
            rows.append({"box_id": b.key, "kind": b.kind, "basis": b.basis, "label": b.label,
                         "level": len(b.chain) - 1, "parent": parent.key if parent is not None else "",
                         "x": round(b.x, 2), "y": round(b.y, 2), "width": round(b.w, 2),
                         "height": round(b.h, 2), "members": len(members(b))})
        return rows

    def stats(self):
        w, h = self.root.w, self.root.h
        return {"nodes": len(self.pos), "boxes": len(self.boxes),
                "singletons_ungrouped": len(self.dropped_singletons),
                "width": round(w), "height": round(h),
                "aspect": round(w / h, 2) if h else 0.0,
                "area_per_node": round(w * h / max(1, len(self.pos)))}


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("features_dir")
    ap.add_argument("out_dir", help="a strategy's groups.tsv / groups_meta.tsv directory")
    ap.add_argument("--option", required=True, choices=sorted(OPT.OPTIONS))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    place = Placement(a.features_dir, a.out_dir, a.option, a.seed)
    out = pathlib.Path(a.out)
    write_tsv(out / "positions.tsv", place.positions_rows())
    write_tsv(out / "boxes.tsv", place.boxes_rows())
    print(f"{out}: " + " · ".join(f"{k} {v}" for k, v in place.stats().items()))


if __name__ == "__main__":
    main()
