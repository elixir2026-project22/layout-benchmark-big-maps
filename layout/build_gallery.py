"""Build the layout gallery: every grouping of every big model, distributed by every layout option.

One picture per (layout option, grouping, model). The same algorithm is used at every level: it places
the species and reactions inside a group, the sub-groups inside a group, and the top-level groups
themselves. Only the four big models are used - the grouping work is about them.

Run:
    python layout/build_gallery.py                      # render everything, write layout/examples/index.html
    python layout/build_gallery.py --jobs 4             # same, on four processes
    python layout/build_gallery.py --no-render          # only rebuild the page from gallery.tsv
    python layout/build_gallery.py --option graphviz-dot --strategy hybrid      # a slice
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import html
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "examples"))
import options as OPT  # noqa: E402

# the four models the grouping work is about
MODELS = ["brenda_merged", "e_coli_core", "R-HSA-77289", "skm_pss"]

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
FIELDS = ["option", "family", "strategy", "model", "svg", "nodes", "boxes", "singletons_ungrouped",
          "width", "height", "aspect", "area_per_node", "seconds", "error"]


def one(job):
    """Place and draw one (option, strategy, model); returns a gallery row."""
    option, strategy, rel, model, seed = job
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(ROOT / "examples"))
    import placement as PL
    import render_layout as RL

    svg_name = f"{option}__{strategy}__{model}.svg"
    row = {"option": option, "family": OPT.OPTIONS[option][1], "strategy": strategy, "model": model,
           "svg": svg_name, "error": ""}
    start = time.time()
    try:
        place = PL.Placement(ROOT / "panel" / "features" / model, ROOT / rel / model, option, seed)
        svg = RL.render(place, f"{option} — {strategy} — {model}")
        (HERE / "examples" / "svg" / svg_name).write_text(svg, encoding="utf-8")
        out = HERE / "out" / option / strategy / model
        PL.write_tsv(out / "positions.tsv", place.positions_rows())
        PL.write_tsv(out / "boxes.tsv", place.boxes_rows())
        row.update(place.stats())
    except Exception as exc:  # a layout tool can fail or time out: the cell says so
        row["error"] = f"{type(exc).__name__}: {exc}"[:200]
    row["seconds"] = round(time.time() - start, 1)
    return {k: row.get(k, "") for k in FIELDS}


def page(rows):
    by = {(r["option"], r["strategy"], r["model"]): r for r in rows}
    opts = [o for o in OPT.OPTIONS if any(r["option"] == o for r in rows)]
    strategies = [s for s in STRATEGIES if any(r["strategy"] == s[0] for r in rows)]
    models = [m for m in MODELS if any(r["model"] == m for r in rows)]

    def cell(option, strategy, model):
        r = by.get((option, strategy, model))
        if r is None:
            return '<td class="none">not run</td>'
        if r["error"]:
            return f'<td class="none">failed<div class="st">{html.escape(str(r["error"])[:120])}</div></td>'
        st = (f'boxes {r["boxes"]} · ungrouped singletons {r["singletons_ungrouped"]}<br>'
              f'{r["width"]}×{r["height"]} · aspect {r["aspect"]} · {r["area_per_node"]} px²/node')
        return (f'<td><a href="svg/{r["svg"]}" target="_blank">'
                f'<img src="svg/{r["svg"]}" loading="lazy" alt="{option} {strategy} {model}"></a>'
                f'<div class="st">{st}</div></td>')

    out = ["<!doctype html><html><head><meta charset='utf-8'><title>Grouping layouts — examples</title><style>",
           "body{font:13px/1.45 Helvetica,Arial,sans-serif;margin:24px;color:#222}",
           "h1{font-size:20px}h2{font-size:15px;margin:30px 0 4px}table{border-collapse:collapse;margin-bottom:8px}",
           "td,th{border:1px solid #ddd;padding:6px;vertical-align:top;text-align:left}",
           "th{background:#f5f5f5;font-weight:600}img{width:300px;border:1px solid #eee;background:#fff}",
           ".st{color:#666;font-size:11px;margin-top:3px}.none{color:#999;font-style:italic}",
           ".desc{color:#555;max-width:62em}nav a{margin-right:12px}code{background:#f5f5f5;padding:0 3px}",
           "</style></head><body>",
           "<h1>Grouping layouts — examples</h1>",
           "<p class='desc'>Each picture is one grouping of one model, laid out by one layout option. The option is "
           "run at every level of the grouping: it places the species and reactions inside a leaf group, the "
           "sub-groups inside a group, and the top-level groups themselves (<code>layout/placement.py</code>), and "
           "the picture is drawn by the same renderer as the grouping gallery (<code>layout/render_layout.py</code>): "
           "species coloured by compartment, reactions by SBO term, compartment levels drawn as a cell, boundary "
           "nodes in a dashed band, edges that leave their top-level group in red. A group holding one entity and "
           "nothing else is not drawn as a group - its node is placed in the parent instead, counted per cell as "
           "<em>ungrouped singletons</em>. Click a picture for the full-size SVG.</p>",
           "<p class='desc'><b>Options</b> (one per family of layout-benchmark's registry): " +
           ", ".join(f'<a href="{OPT.OPTIONS[o][2]}"><code>{o}</code></a> ({OPT.OPTIONS[o][1]})' for o in opts) +
           ". Not available on the machine that built this page: " +
           "; ".join(f"{html.escape(f)} ({html.escape(w)})" for f, w in OPT.UNAVAILABLE.items()) + ".</p>",
           "<nav>" + " ".join(f'<a href="#{s[0]}">{s[0]}</a>' for s in strategies) + "</nav>"]
    for name, _rel, desc in strategies:
        out.append(f'<h2 id="{name}">{html.escape(name)}</h2><p class="desc">{html.escape(desc)}</p>')
        out.append("<table><tr><th>layout option</th>" + "".join(f"<th>{html.escape(m)}</th>" for m in models)
                   + "</tr>")
        for option in opts:
            out.append(f'<tr><th>{html.escape(option)}<div class="st">{html.escape(OPT.OPTIONS[option][1])}</div></th>'
                       + "".join(cell(option, name, m) for m in models) + "</tr>")
        out.append("</table>")
    out.append("</body></html>")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-render", action="store_true", help="rebuild index.html from gallery.tsv")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--option", action="append", help="only these options (repeatable)")
    ap.add_argument("--strategy", action="append", help="only these strategies (repeatable)")
    ap.add_argument("--model", action="append", help="only these models (repeatable)")
    a = ap.parse_args(argv)
    (HERE / "examples" / "svg").mkdir(parents=True, exist_ok=True)
    tsv = HERE / "examples" / "gallery.tsv"

    if a.no_render:
        with open(tsv, encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    else:
        jobs = [(o, s, rel, m, a.seed)
                for o in OPT.OPTIONS if not a.option or o in a.option
                for s, rel, _d in STRATEGIES if not a.strategy or s in a.strategy
                for m in MODELS if not a.model or m in a.model
                if (ROOT / rel / m / "groups.tsv").exists()]
        print(f"{len(jobs)} pictures, {a.jobs} job(s)")
        rows = []
        if a.jobs > 1:
            # threads, not processes: the work is in Graphviz subprocesses and in igraph/numpy,
            # which all release the GIL, and a sandboxed Windows host may refuse the process pool
            with cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
                for row in pool.map(one, jobs):
                    rows.append(row)
                    print(f"  {row['option']:18s}{row['strategy']:20s}{row['model']:16s}"
                          f"{row['seconds']:6}s {row['error']}")
        else:
            for job in jobs:
                row = one(job)
                rows.append(row)
                print(f"  {row['option']:18s}{row['strategy']:20s}{row['model']:16s}"
                      f"{row['seconds']:6}s {row['error']}")
        rows.sort(key=lambda r: (list(OPT.OPTIONS).index(r["option"]),
                                 [s[0] for s in STRATEGIES].index(r["strategy"]),
                                 MODELS.index(r["model"])))
        with open(tsv, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    (HERE / "examples" / "index.html").write_text(page(rows), encoding="utf-8")
    failed = sum(1 for r in rows if r["error"])
    print(f"layout/examples/index.html: {len(rows)} pictures, {failed} failed")


if __name__ == "__main__":
    main()
