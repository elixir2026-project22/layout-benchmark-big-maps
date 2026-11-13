"""Build the examples gallery: every model under every strategy, drawn by the one common renderer.

Run:
    python examples/build_gallery.py                     # re-renders every SVG and writes examples/index.html
    python examples/build_gallery.py --no-render         # only rebuild the page from the SVGs already there
"""
from __future__ import annotations

import argparse
import csv
import html
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import render_grouping as RG  # noqa: E402

STRATEGIES = [
    ("type_location", "strategies/type_location/out", "reaction type on top, localisation inside"),
    ("location_type", "strategies/location_type/out", "localisation on top, reaction type inside"),
    ("ontology", "strategies/ontology/out", "annotation trees (notes bins, EC, GO, xrefs), SBO fallback"),
    ("topology_leiden", "strategies/topology/topology_leiden/out", "Leiden communities (modularity)"),
    ("topology_infomap", "strategies/topology/topology_infomap/out", "Infomap communities (flow)"),
    ("topology_walktrap", "strategies/topology/topology_walktrap/out", "Walktrap communities (random walk)"),
    ("cascade", "strategies/cascade/out", "flow-anchored modules, stages inside"),
    ("hybrid", "strategies/hybrid/out", "function, then type, then localisation, size-bounded"),
    ("baseline", "baseline", "the project's current grouping (pathway-grouping)"),
]
MODELS = ["brenda_merged", "e_coli_core", "R-HSA-77289", "skm_ja", "skm_sa", "skm_pss"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def stats(out_dir):
    groups, meta = read_tsv(out_dir / "groups.tsv"), read_tsv(out_dir / "groups_meta.tsv")
    levels = {lv: {g[lv] for g in groups} - {""} for lv in ("L1", "L2", "L3")}
    size = {}
    for g in groups:
        size[g["L1"]] = size.get(g["L1"], 0) + 1
    named = sum(1 for m in meta if m["name_source"] not in ("none", ""))
    leaf = {}
    for g in groups:
        key = g["L3"] or g["L2"] or g["L1"]
        leaf[key] = leaf.get(key, 0) + 1
    return {"nodes": len(groups), "L1": len(levels["L1"]), "L2": len(levels["L2"]), "L3": len(levels["L3"]),
            "groups": len(meta), "named": f"{named}/{len(meta)}",
            "largest L1": f"{max(size.values()) / len(groups):.0%}",
            "singleton leaves": sum(1 for v in leaf.values() if v == 1),
            "clones": sum(1 for g in groups if g["clone_in"])}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-render", action="store_true")
    a = ap.parse_args(argv)
    (HERE / "svg").mkdir(exist_ok=True)
    scores = {}
    score_file = ROOT / "scoring" / "scores.csv"
    if score_file.exists():
        for r in csv.DictReader(open(score_file, encoding="utf-8", newline="")):
            scores[(r.get("strategy", ""), r.get("model", ""))] = r
    rows, cells = [], {}
    for name, rel, _desc in STRATEGIES:
        for model in MODELS:
            out = ROOT / rel / model
            if not (out / "groups.tsv").exists():
                cells[(name, model)] = None
                continue
            svg = HERE / "svg" / f"{name}__{model}.svg"
            if not a.no_render:
                svg.write_text(RG.render(ROOT / "panel" / "features" / model, out, f"{name} — {model}"), encoding="utf-8")
            st = stats(out)
            sc = scores.get((name, model), {})
            cells[(name, model)] = {"svg": svg.name, "stats": st, "score": sc.get("score", "")}
            rows.append({"strategy": name, "model": model, "svg": svg.name, **st, "score": sc.get("score", "")})
    with open(HERE / "gallery.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    def cell(name, model):
        c = cells[(name, model)]
        if c is None:
            return '<td class="none">not run</td>'
        st = " · ".join(f"{k} {v}" for k, v in c["stats"].items() if k in ("L1", "L2", "L3", "largest L1"))
        score = f'<div class="score">score {float(c["score"]):.3f}</div>' if c["score"] else ""
        return (f'<td><a href="svg/{c["svg"]}" target="_blank"><img src="svg/{c["svg"]}" loading="lazy" alt="{name} {model}"></a>'
                f'{score}<div class="st">{st}</div></td>')

    out = ["<!doctype html><html><head><meta charset='utf-8'><title>Grouping strategies — examples</title><style>",
           "body{font:13px/1.45 Helvetica,Arial,sans-serif;margin:24px;color:#222}",
           "h1{font-size:20px}h2{font-size:15px;margin:28px 0 6px}table{border-collapse:collapse}",
           "td,th{border:1px solid #ddd;padding:6px;vertical-align:top;text-align:left}",
           "th{background:#f5f5f5;font-weight:600}img{width:300px;border:1px solid #eee;background:#fff}",
           ".st{color:#666;font-size:11px;margin-top:3px}.score{font-weight:600;font-size:12px;margin-top:4px}",
           ".none{color:#999;font-style:italic}.desc{color:#555}</style></head><body>",
           "<h1>Grouping strategies — examples</h1>",
           "<p class='desc'>Every model drawn under every strategy by one renderer (<code>examples/render_grouping.py</code>): "
           "nested boxes for functional, type, topology and flow groups; wherever a level is localisation it is drawn as a cell "
           "from the compartment containment tree — extracellular region around the cytosol, organelles inside it, transport "
           "steps in a dashed band on the boundary they cross. Click a picture for the full-size SVG. "
           "L1/L2/L3 = number of groups per level, 'largest L1' = share of the model in the biggest top-level box.</p>",
           "<table><tr><th>strategy</th>" + "".join(f"<th>{html.escape(m)}</th>" for m in MODELS) + "</tr>"]
    for name, _rel, desc in STRATEGIES:
        out.append(f'<tr><th>{html.escape(name)}<div class="st">{html.escape(desc)}</div></th>'
                   + "".join(cell(name, m) for m in MODELS) + "</tr>")
    out.append("</table></body></html>")
    (HERE / "index.html").write_text("\n".join(out), encoding="utf-8")
    print(f"examples/index.html: {len(rows)} pictures")


if __name__ == "__main__":
    main()
