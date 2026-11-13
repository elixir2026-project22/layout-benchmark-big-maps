"""Draw one model, one grouping and one layout option: the same picture as examples/render_grouping.py,
with the grid and shelf packing replaced by the layout option.

Colours, shapes and the cell drawing are taken from ``examples/render_grouping.py`` so that a picture
here and a picture there differ only by where things are put: species boxes coloured by compartment,
reactions as squares coloured by SBO term, spaces as rounded compartment boxes, boundary bands dashed,
edges red when they leave their top-level group, SBGN clone markers on cloned species.

Run:
    python layout/render_layout.py panel/features/skm_pss strategies/hybrid/out/skm_pss \
        --option graphviz-sfdp -o layout/examples/svg/graphviz-sfdp__hybrid__skm_pss.svg
"""
from __future__ import annotations

import argparse
import math
import pathlib
import sys
from xml.sax.saxutils import escape

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "examples"))
import placement as PL  # noqa: E402
import render_grouping as RG  # noqa: E402


def render(place: PL.Placement, title: str = "") -> str:
    nodes, comp, meta = place.nodes, place.comp, place.meta
    sbo_terms = sorted({n["sbo"] for n in nodes.values() if n["kind"] == "reaction" and n["sbo"]})
    rx_fill = {t: RG.RX_COLOURS[i % len(RG.RX_COLOURS)] for i, t in enumerate(sbo_terms)}
    sbo_name = {n["sbo"]: n["sbo_name"] for n in nodes.values() if n["sbo"]}
    spaces = [c for c in comp.values() if c["kind"] == "space"]

    def space_fill(cid):
        return RG.SPACE_FILL.get(comp.get(cid, {}).get("canonical", ""), "#f6f6f6")

    legend = [(RG.darker(space_fill(c["compartment_id"]), 0.9), c["name"]) for c in spaces] if len(spaces) > 1 else []
    legend += [(rx_fill[t], f"{sbo_name.get(t, t)}") for t in sbo_terms] if len(sbo_terms) > 1 else []
    root = place.root
    leg_cols = max(1, min(4, int(max(root.w, 600) // 230)))
    leg_h = math.ceil(len(legend) / leg_cols) * 16 + (30 if title else 8) if legend or title else 0
    W, H = max(root.w, 600.0), root.h + leg_h
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
           f'viewBox="0 0 {W:.0f} {H:.0f}" font-family="Helvetica, Arial, sans-serif">',
           f"<title>{escape(title)}</title>",
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

    for b in place.boxes:
        depth = len(b.chain) - 1
        name = b.label if b.label else "[unnamed]"
        if b.kind == "space":
            cid = (b.key.split("space:")[-1].split("@")[0] if "space:" in b.key
                   else meta[b.key]["compartment_id"].split("@")[0])
            fill = space_fill(cid)
            dash = ' stroke-dasharray="5 3"' if b.frame else ""
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="16" fill="{fill}" '
                       f'stroke="{RG.darker(fill, 0.6)}" stroke-width="1.8"{dash}/>')
            style = 'font-style="italic" fill="#444"'
        elif b.kind == "band":
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="4" fill="#ffffff" '
                       f'fill-opacity="0.55" stroke="#555" stroke-width="1" stroke-dasharray="2 2"/>')
            style = 'fill="#666"'
        else:
            shade = ["#ffffff", "#f7f7f7", "#efefef"][min(depth - 1, 2) % 3]
            out.append(f'<rect x="{b.x:.1f}" y="{b.y:.1f}" width="{b.w:.1f}" height="{b.h:.1f}" rx="3" fill="{shade}" '
                       f'stroke="#333" stroke-width="{1.6 if depth == 1 else 1.0}"/>')
            style = 'fill="#111" font-weight="bold"' if b.label else 'fill="#999" font-style="italic"'
        size = 12 if depth == 1 else 10
        fit = max(3, int((b.w - 10) / (size * 0.56)))
        text = name if len(name) <= fit else name[: fit - 1] + "…"
        out.append(f'<text x="{b.x + 6:.1f}" y="{b.y + size + 3:.1f}" font-size="{size}" {style}>{escape(text)}'
                   f"<title>{escape(name)} [{escape(b.basis or b.kind)}]</title></text>")

    pos = place.pos
    top_of = {k: (place.home[k].chain[1].key if len(place.home[k].chain) > 1 else "") for k in pos}
    for e in place.edges:
        r, s = e["reaction_id"], place._drawn(e["species_id"], e["reaction_id"])
        if r not in pos or s not in pos:
            continue
        (x1, y1), (x2, y2) = (pos[s], pos[r]) if e["role"] != "product" else (pos[r], pos[s])
        cross = top_of.get(r) != top_of.get(s)
        style = (f'stroke="{"#c0392b" if cross else "#777"}" stroke-width="0.7" '
                 f'stroke-opacity="{0.35 if cross else 0.7}"')
        if e["role"] == "modifier":
            style += ' stroke-dasharray="3 2"'
        marker = ' marker-end="url(#arr)"' if e["role"] == "product" else ""
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {style}{marker}/>')

    for key, (x, y) in sorted(pos.items()):
        n = nodes[key.split("#")[0]]
        tip = escape(f"{n['name'] or n['node_id']} ({n['node_id']}) {n['sbo_name']}")
        if n["kind"] == "reaction":
            out.append(f'<rect x="{x - RG.RX / 2:.1f}" y="{y - RG.RX / 2:.1f}" width="{RG.RX}" height="{RG.RX}" '
                       f'fill="{rx_fill.get(n["sbo"], "#fff")}" stroke="#222" stroke-width="0.8">'
                       f"<title>{tip}</title></rect>")
            continue
        fill = RG.darker(space_fill(n["loc_space"]), 0.93) if len(spaces) > 1 else "#ffffff"
        rx = RG.SP_H / 2 if n["sbo"] == "SBO:0000247" else 3
        out.append(f'<rect x="{x - RG.SP_W / 2:.1f}" y="{y - RG.SP_H / 2:.1f}" width="{RG.SP_W}" '
                   f'height="{RG.SP_H}" rx="{rx}" fill="{fill}" stroke="#222" stroke-width="0.8">'
                   f"<title>{tip}</title></rect>")
        if "#" in key:  # SBGN clone marker
            out.append(f'<rect x="{x - RG.SP_W / 2 + 2:.1f}" y="{y + RG.SP_H * 0.2:.1f}" width="{RG.SP_W - 4}" '
                       f'height="{RG.SP_H * 0.26:.1f}" fill="#555"/>')
        label = n["name"] or n["node_id"]
        label = label if len(label) <= 11 else label[:10] + "…"
        out.append(f'<text x="{x:.1f}" y="{y + 3:.1f}" font-size="8.5" text-anchor="middle" fill="#111">'
                   f"{escape(label)}</text>")
    out += ["</g>", "</svg>"]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("features_dir")
    ap.add_argument("out_dir", help="a strategy's groups.tsv / groups_meta.tsv directory")
    ap.add_argument("--option", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--title", default="")
    a = ap.parse_args(argv)
    place = PL.Placement(a.features_dir, a.out_dir, a.option, a.seed)
    svg = render(place, a.title)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(a.out).write_text(svg, encoding="utf-8")
    print(f"{a.out}: {len(svg) / 1e3:.0f} kB · " + " · ".join(f"{k} {v}" for k, v in place.stats().items()))


if __name__ == "__main__":
    main()
