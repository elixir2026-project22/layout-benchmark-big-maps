"""Draw one model under one grouping: nested group boxes, localisation as a cell.

One renderer for every strategy, so pictures differ only by the grouping. Function, type, topology and
flow groups are plain nested boxes. Wherever a level is localisation, its groups are drawn as a cell
from compartments.tsv: the extracellular region around the cytosol, organelles inside the cytosol.
Nodes on a boundary (transport across it, membrane participants) sit in a band on the border of
their space. Species are coloured by compartment, reactions by SBO term.

Run:
    python examples/render_grouping.py panel/features/skm_ja strategies/<name>/out/skm_ja \
        -o examples/<name>/skm_ja.svg --title "<name> - skm_ja"
"""
from __future__ import annotations

import argparse
import collections
import csv
import math
import pathlib
import sys
from xml.sax.saxutils import escape

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
from compartments import Tree  # noqa: E402

CELL_W, CELL_H = 74.0, 34.0          # grid cell of a node
SP_W, SP_H, RX = 62.0, 20.0, 11.0    # species box, reaction square
PAD, HEAD, GAP = 12.0, 20.0, 14.0
ASPECT = 1.5
SPACE_FILL = {"extracellular": "#f1f1f1", "apoplast": "#e9f0df", "cell_wall": "#e5e9d8", "periplasm": "#eef0e4",
              "cytosol": "#fff8e1", "nucleus": "#e1ebfb", "nucleoplasm": "#e1ebfb", "mitochondrion": "#fbe4dd",
              "mito_matrix": "#fbe4dd", "chloroplast": "#dff3df", "er": "#f0e4f8", "er_lumen": "#f0e4f8",
              "golgi": "#f8ead3", "peroxisome": "#def2f2", "vacuole": "#e9e9f8", "lysosome": "#e9e9f8"}
RX_COLOURS = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1", "#76b7b2", "#edc948", "#9c755f", "#ff9da7",
              "#bab0ac", "#1b9e77", "#7570b3"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def darker(hex_colour, f=0.82):
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return f"#{int(r * f):02x}{int(g * f):02x}{int(b * f):02x}"


class Box:
    def __init__(self, key, label, kind, basis=""):
        self.key, self.label, self.kind, self.basis = key, label, kind, basis  # kind: root | group | space | band
        self.children, self.nodes, self.band = [], [], None
        self.canonical, self.frame = "", False
        self.x = self.y = self.w = self.h = 0.0
        self.grid = {}  # drawn node -> (dx, dy) inside the box
        self.depth = 0


def build(features_dir, out_dir):
    f, o = pathlib.Path(features_dir), pathlib.Path(out_dir)
    nodes = {r["node_id"]: r for r in read_tsv(f / "nodes.tsv")}
    edges = read_tsv(f / "edges.tsv")
    comp_rows = read_tsv(f / "compartments.tsv")
    tree, comp = Tree(comp_rows), {r["compartment_id"]: r for r in comp_rows}
    groups = read_tsv(o / "groups.tsv")
    meta = {r["group_id"]: r for r in read_tsv(o / "groups_meta.tsv")}

    root = Box("", "", "root")
    box = {"": root}
    for gid in sorted(meta, key=lambda g: (int(meta[g]["level"]), g)):
        m = meta[gid]
        b = Box(gid, m["label"], "group", m["basis"])
        box[gid] = b
        box[m["parent"]].children.append(b)
    home = {}
    for g in groups:
        chain = [g[lv] for lv in ("L1", "L2", "L3") if g[lv]]
        box[chain[-1]].nodes.append(g["node_id"])
        home[g["node_id"]] = box[chain[-1]]
    clones = collections.defaultdict(list)  # species -> [(box, drawn key)]
    for g in groups:
        for c in filter(None, g["clone_in"].split("|")):
            key = f"{g['node_id']}#{c}"
            box[c].nodes.append(key)
            clones[g["node_id"]].append((box[c], key))

    def to_cell(parent):
        """Replace the localisation groups under `parent` by spaces nested as compartments.tsv says."""
        loc = [c for c in parent.children if c.basis == "localisation"]
        for c in parent.children:
            if c.basis != "localisation":
                to_cell(c)
        if not loc:
            return
        inside, edge = {}, {}
        for b in loc:
            cid, _, on_edge = meta[b.key]["compartment_id"].partition("@")
            table = edge if on_edge else inside
            if cid in table:  # several groups of one compartment under one parent: they nest inside its space
                first = table[cid]
                if first.kind != "space":
                    holder = Box(f"{parent.key}|space:{cid}{'@boundary' if on_edge else ''}", "", "space")
                    holder.frame, holder.children = True, [first]
                    table[cid] = holder
                table[cid].children.append(b)
            else:
                table[cid] = b
            to_cell(b)
        used = sorted(set(inside) | set(edge))
        paths = [tree.ancestors(c)[::-1] for c in used]  # root ... c
        common = 0
        while all(len(p) > common for p in paths) and len({p[common] for p in paths}) == 1:
            common += 1
        top = paths[0][common - 1]
        needed = {c for p in paths for c in p[common - 1:]}
        space = {}
        for cid in sorted(needed, key=lambda c: (len(tree.ancestors(c)), c)):
            b = inside.get(cid) or Box(f"{parent.key}|space:{cid}", "", "space")
            b.frame = cid not in inside or b.key.endswith(f"|space:{cid}")
            b.kind, b.label, b.canonical = "space", comp[cid]["name"], comp[cid]["canonical"]
            space[cid] = b
            if cid in edge:
                band = edge[cid]
                band.kind, band.label = "band", "boundary"
                b.band = band
            elif cid in inside and not b.children:
                on = [n for n in b.nodes if nodes[n.split("#")[0]]["loc_kind"] == "boundary"]
                if on and len(on) < len(b.nodes):
                    band = Box(f"{b.key}|band", "boundary", "band")
                    band.nodes = on
                    b.nodes = [n for n in b.nodes if n not in set(on)]
                    for n in on:
                        if "#" not in n:
                            home[n] = band
                    b.band = band
        parent.children = [c for c in parent.children if c.basis != "localisation"]
        for cid, b in space.items():
            if cid == top:
                parent.children.append(b)
            else:
                up = next(a for a in tree.ancestors(cid)[1:] if a in space)
                space[up].children.append(b)

    to_cell(root)
    return root, nodes, edges, comp, tree, home, clones, meta


def order_nodes(keys, adj):
    """Breadth-first over the edges inside the box, so neighbours sit next to each other."""
    inside = set(keys)
    seen, out = set(), []
    deg = {k: sum(1 for x in adj.get(k.split("#")[0], ()) if x in inside) for k in keys}
    for seed in sorted(keys, key=lambda k: (-deg[k], k)):
        if seed in seen:
            continue
        queue = collections.deque([seed])
        seen.add(seed)
        while queue:
            k = queue.popleft()
            out.append(k)
            for x in sorted(adj.get(k.split("#")[0], ())):
                if x in inside and x not in seen:
                    seen.add(x)
                    queue.append(x)
    return out


def grid_size(n, cols=None):
    cols = cols or max(1, math.ceil(math.sqrt(n * ASPECT * CELL_H / CELL_W)))
    return cols, math.ceil(n / cols)


def measure(b, adj, depth=0):
    b.depth = depth
    for c in b.children:
        measure(c, adj, depth + 1)
    items = []  # (key, w, h)
    if b.nodes:
        cols, rows = grid_size(len(b.nodes))
        if b.kind == "band":
            cols, rows = len(b.nodes), 1
        items.append(("__grid__", cols * CELL_W, rows * CELL_H, cols))
    for c in sorted(b.children, key=lambda c: (-(c.w * c.h), c.key)):
        items.append((c, c.w, c.h, 0))
    band_h = 0.0
    if b.band is not None:
        measure(b.band, adj, depth + 1)
    area = sum(w * h for _, w, h, _ in items)
    limit = max([math.sqrt(area * ASPECT)] + [w for _, w, _, _ in items]) if items else 0.0
    if b.band is not None:  # the band runs along the top border: wrap it to the content width
        limit = max(limit, min(b.band.w, 6 * CELL_W))
        cols = max(1, int(limit // CELL_W))
        if b.band.nodes and not b.band.children:
            rows = math.ceil(len(b.band.nodes) / cols)
            b.band.w, b.band.h = cols * CELL_W + 2 * PAD, rows * CELL_H + PAD + HEAD * 0.8
            b.band.cols = cols
        band_h = b.band.h + GAP / 2
    x = y = row_h = width = 0.0
    place = []
    for key, w, h, cols in items:
        if x > 0 and x + w > limit + 1e-6:
            x, y, row_h = 0.0, y + row_h + GAP, 0.0
        place.append((key, x, y, cols))
        x += w + GAP
        row_h = max(row_h, h)
        width = max(width, x - GAP)
    height = y + row_h
    if b.band is not None:
        width = max(width, b.band.w - 2 * PAD if b.band.nodes and not b.band.children else b.band.w)
    head = HEAD if b.kind != "root" else 0.0
    b.place, b.band_h = place, band_h
    b.w = width + 2 * PAD
    b.h = height + 2 * PAD + head + band_h
    if b.kind == "band" and b.nodes and not b.children:
        b.cols = getattr(b, "cols", len(b.nodes))
        rows = math.ceil(len(b.nodes) / b.cols)
        b.w, b.h = b.cols * CELL_W + 2 * PAD, rows * CELL_H + PAD + HEAD * 0.8
    b.nodes = order_nodes(b.nodes, adj) if b.nodes else []


def position(b, x, y, pos):
    b.x, b.y = x, y
    head = HEAD if b.kind != "root" else 0.0
    if b.kind == "band" and b.nodes and not b.children:
        for i, n in enumerate(b.nodes):
            pos[n] = (x + PAD + (i % b.cols) * CELL_W + CELL_W / 2, y + HEAD * 0.8 + (i // b.cols) * CELL_H + CELL_H / 2)
        return
    top = y + PAD + head
    if b.band is not None:
        b.band.w = b.w - 2 * PAD if not b.band.children else b.band.w
        if b.band.nodes and not b.band.children:
            b.band.cols = max(1, int((b.band.w - 2 * PAD) // CELL_W))
            b.band.h = math.ceil(len(b.band.nodes) / b.band.cols) * CELL_H + PAD + HEAD * 0.8
        position(b.band, x + PAD, top, pos)
        top += b.band_h
    for key, dx, dy, cols in b.place:
        if key == "__grid__":
            for i, n in enumerate(b.nodes):
                pos[n] = (x + PAD + dx + (i % cols) * CELL_W + CELL_W / 2, top + dy + (i // cols) * CELL_H + CELL_H / 2)
        else:
            position(key, x + PAD + dx, top + dy, pos)


def render(features_dir, out_dir, title=""):
    root, nodes, edges, comp, tree, home, clones, meta = build(features_dir, out_dir)
    adj = collections.defaultdict(set)
    for e in edges:
        adj[e["reaction_id"]].add(e["species_id"])
        adj[e["species_id"]].add(e["reaction_id"])
    measure(root, adj)
    pos = {}
    position(root, 0.0, 0.0, pos)

    boxes = []

    def walk(b, chain):
        if b.kind != "root":
            boxes.append(b)
        b.chain = chain + [b]
        for c in ([b.band] if b.band is not None else []) + b.children:
            walk(c, b.chain)

    walk(root, [])
    sbo_terms = sorted({n["sbo"] for n in nodes.values() if n["kind"] == "reaction" and n["sbo"]})
    rx_fill = {t: RX_COLOURS[i % len(RX_COLOURS)] for i, t in enumerate(sbo_terms)}
    sbo_name = {n["sbo"]: n["sbo_name"] for n in nodes.values() if n["sbo"]}
    spaces = [c for c in comp.values() if c["kind"] == "space"]

    def space_fill(cid):
        return SPACE_FILL.get(comp.get(cid, {}).get("canonical", ""), "#f6f6f6")

    legend = [(darker(space_fill(c["compartment_id"]), 0.9), c["name"]) for c in spaces] if len(spaces) > 1 else []
    legend += [(rx_fill[t], f"{sbo_name.get(t, t)}") for t in sbo_terms] if len(sbo_terms) > 1 else []
    leg_cols = max(1, min(4, int(max(root.w, 600) // 230)))
    leg_h = math.ceil(len(legend) / leg_cols) * 16 + (30 if title else 8) if legend or title else 0
    W, H = max(root.w, 600.0), root.h + leg_h
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" '
           f'font-family="Helvetica, Arial, sans-serif">', f"<title>{escape(title)}</title>",
           '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" '
           'orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#666"/></marker></defs>',
           '<rect width="100%" height="100%" fill="#fff"/>']
    if title:
        out.append(f'<text x="10" y="18" font-size="14" font-weight="bold" fill="#222">{escape(title)}</text>')
    for i, (colour, text) in enumerate(legend):
        lx, ly = 10 + (i % leg_cols) * 230, (26 if title else 4) + (i // leg_cols) * 16
        out.append(f'<rect x="{lx}" y="{ly}" width="12" height="10" fill="{colour}" stroke="#666" stroke-width="0.5"/>'
                   f'<text x="{lx + 17}" y="{ly + 9}" font-size="10" fill="#333">{escape(text[:36])}</text>')
    out.append(f'<g transform="translate(0,{leg_h:.0f})">')
    for b in boxes:
        name = b.label if b.label else "[unnamed]"
        if b.kind == "space":
            cid = (b.key.split("space:")[-1].split("@")[0] if "space:" in b.key
                   else meta[b.key]["compartment_id"].split("@")[0])
            fill = space_fill(cid)
            dash = ' stroke-dasharray="5 3"' if b.frame else ""
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="16" fill="{fill}" '
                       f'stroke="{darker(fill, 0.6)}" stroke-width="1.8"{dash}/>')
            style = 'font-style="italic" fill="#444"'
        elif b.kind == "band":
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="4" fill="#ffffff" '
                       f'fill-opacity="0.55" stroke="#555" stroke-width="1" stroke-dasharray="2 2"/>')
            style = 'fill="#666"'
        else:
            shade = ["#ffffff", "#f7f7f7", "#efefef"][min(b.depth - 1, 2) % 3]
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="3" fill="{shade}" '
                       f'stroke="#333" stroke-width="{1.6 if b.depth == 1 else 1.0}"/>')
            style = ('fill="#111" font-weight="bold"' if b.label else 'fill="#999" font-style="italic"')
        size = 12 if b.depth == 1 else 10
        fit = max(3, int((b.w - 10) / (size * 0.56)))
        text = name if len(name) <= fit else name[: fit - 1] + "…"
        out.append(f'<text x="{b.x + 6:.1f}" y="{b.y + size + 3:.1f}" font-size="{size}" {style}>{escape(text)}'
                   f"<title>{escape(name)} [{escape(b.basis or b.kind)}]</title></text>")

    def drawn(species, reaction):
        """The copy of a species a reaction connects to: a clone in the reaction's own boxes, else its home."""
        chain = {id(x) for x in home[reaction].chain}
        best = None
        for b, key in clones.get(species, ()):
            if id(b) in chain and (best is None or b.depth > best[0]):
                best = (b.depth, key)
        return best[1] if best else species

    top_of = {}
    for b in boxes:
        for n in b.nodes:
            top_of[n] = b.chain[1].key if len(b.chain) > 1 else ""
    for e in edges:
        r, s = e["reaction_id"], drawn(e["species_id"], e["reaction_id"])
        if r not in pos or s not in pos:
            continue
        (x1, y1), (x2, y2) = (pos[s], pos[r]) if e["role"] != "product" else (pos[r], pos[s])
        cross = top_of.get(r) != top_of.get(s)
        style = f'stroke="{"#c0392b" if cross else "#777"}" stroke-width="0.7" stroke-opacity="{0.35 if cross else 0.7}"'
        if e["role"] == "modifier":
            style += ' stroke-dasharray="3 2"'
        marker = ' marker-end="url(#arr)"' if e["role"] == "product" else ""
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {style}{marker}/>')
    for key, (x, y) in sorted(pos.items()):
        n = nodes[key.split("#")[0]]
        tip = escape(f"{n['name'] or n['node_id']} ({n['node_id']}) {n['sbo_name']}")
        if n["kind"] == "reaction":
            out.append(f'<rect x="{x - RX / 2:.1f}" y="{y - RX / 2:.1f}" width="{RX}" height="{RX}" '
                       f'fill="{rx_fill.get(n["sbo"], "#fff")}" stroke="#222" stroke-width="0.8"><title>{tip}</title></rect>')
            continue
        fill = darker(space_fill(n["loc_space"]), 0.93) if len(spaces) > 1 else "#ffffff"
        rx = SP_H / 2 if n["sbo"] == "SBO:0000247" else 3
        out.append(f'<rect x="{x - SP_W / 2:.1f}" y="{y - SP_H / 2:.1f}" width="{SP_W}" height="{SP_H}" rx="{rx}" '
                   f'fill="{fill}" stroke="#222" stroke-width="0.8"><title>{tip}</title></rect>')
        if "#" in key:  # SBGN clone marker
            out.append(f'<rect x="{x - SP_W / 2 + 2:.1f}" y="{y + SP_H * 0.2:.1f}" width="{SP_W - 4}" height="{SP_H * 0.26:.1f}" '
                       f'fill="#555"/>')
        label = n["name"] or n["node_id"]
        label = label if len(label) <= 11 else label[:10] + "…"
        out.append(f'<text x="{x:.1f}" y="{y + 3:.1f}" font-size="8.5" text-anchor="middle" fill="#111">{escape(label)}</text>')
    out += ["</g>", "</svg>"]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("features_dir")
    ap.add_argument("out_dir")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--title", default="")
    a = ap.parse_args(argv)
    svg = render(a.features_dir, a.out_dir, a.title)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(a.out).write_text(svg, encoding="utf-8")
    print(f"{a.out}: {len(svg) / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
