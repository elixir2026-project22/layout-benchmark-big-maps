# layout

Distributing a grouped model: the elements inside each group, and the groups themselves.

A grouping (`../strategies/<name>/out/<model>/`) says *what belongs together*; it says nothing about
*where things go*. This module takes a grouping and a layout option from
[layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark) and places every species,
reaction and group with it. BioHackathon Europe 2026, project 22 (*Open Pathway Maps*).

The options are reimplemented here against the Graphviz binaries and igraph, so nothing of
layout-benchmark is imported; what is taken from it is which algorithm stands for which family. Nothing
here changes `../strategies/`, `../panel/` or `../examples/`.

## The protocol

1. **Take the grouping as given.** `../examples/render_grouping.py` already turns `groups.tsv` and
   `groups_meta.tsv` into a tree of boxes, with localisation levels drawn as a cell from
   `compartments.tsv` and boundary nodes in a band. That tree is reused unchanged, so a picture here
   and a picture in `../examples/` differ only by where things are put.
2. **A group of one entity is not a group.** A box holding exactly one species or reaction and no
   sub-group is dropped and its node placed in the parent instead. Each picture reports how many
   (`ungrouped singletons`). Compartment boxes are never dropped - a compartment is biology, not a
   grouping choice.
3. **Distribute the elements within each group.** Inside a box, its own species and reactions and its
   already-sized sub-boxes are one *item graph*: two items are joined when the model has an edge
   between any member of one and any member of the other, weighted by how many such edges there are.
   The layout option places that graph.
4. **Distribute the groups with the same logic.** The top level is the same call on the quotient graph
   of the top-level groups. One option therefore means one algorithm, run at every level.
5. **No box overlaps another.** Each connected component of a level is laid out on its own and the
   components are packed in rows (what Graphviz's `pack=true` does for its own layouts); items are
   then pushed apart, and finally scaled until nothing overlaps. Scaling positions while the sizes
   stay fixed can never bring two items closer, so that last step always terminates and never
   creates a new overlap.

## The options

One algorithm per family of `layout_benchmark.layouts.REGISTRY`, in the registry's own order:

| option | family |
|---|---|
| `graphviz-dot` | layered |
| `igraph-sugiyama` | layered (stands in for *tree*, which needs Node) |
| `graphviz-neato` | stress |
| `graphviz-sfdp` | force-directed |
| `igraph-drl` | for large graphs |
| `graphviz-twopi` | radial |
| `graphviz-circo` | circular |
| `igraph-mds` | baseline |

`python layout/options.py` prints this table and says which tools are present. The registry's other
families need software this module cannot assume: every cytoscape.js / ELK / graphology / Reagraph
option needs Node, and the pathway-specific tools (`rules-claude`, `ecmapper`, `metacarto`,
`sbml2map`) and `sbmlnetwork` need their own checkouts under `.tools/`. Adding one is a new entry in
`OPTIONS` in `options.py`; nothing else changes.

## Files

| file | contents |
|---|---|
| `options.py` | the eight options, each `(items, sizes, edges, seed) -> positions` |
| `placement.py` | the recursive placement, and the `positions.tsv` / `boxes.tsv` tables |
| `render_layout.py` | one SVG, drawn exactly as `../examples/render_grouping.py` draws a grouping |
| `build_gallery.py` | every option × grouping × model, and `examples/index.html` |
| `out/<option>/<strategy>/<model>/` | `positions.tsv` (where each node is drawn) and `boxes.tsv` (where each group box is) |
| `examples/` | `index.html`, `gallery.tsv` and the SVGs |

`positions.tsv`: `node_id`, `drawn_key` (a clone is `<node_id>#<group_id>`), `clone`, `kind`, `box`,
`x`, `y` (centre, SVG coordinates, y down), `width`, `height`.
`boxes.tsv`: `box_id`, `kind` (`group`, `space`, `band`), `basis`, `label`, `level`, `parent`, `x`,
`y` (top-left corner), `width`, `height`, `members`.

## Run

```bash
pip install python-libsbml numpy scipy igraph          # plus Graphviz on the PATH

# which options can run here
python layout/options.py

# one model, one grouping, one option: the tables
python layout/placement.py panel/features/skm_pss strategies/hybrid/out/skm_pss \
    --option graphviz-neato -o layout/out/graphviz-neato/hybrid/skm_pss

# the same, as a picture
python layout/render_layout.py panel/features/skm_pss strategies/hybrid/out/skm_pss \
    --option graphviz-neato -o layout/examples/svg/graphviz-neato__hybrid__skm_pss.svg \
    --title "graphviz-neato — hybrid — skm_pss"

# the whole gallery (8 options × 9 groupings × 4 models)
python layout/build_gallery.py --jobs 6
python layout/build_gallery.py --no-render                 # rebuild the page from gallery.tsv
python layout/build_gallery.py --option graphviz-dot --strategy hybrid --model skm_pss
```

Every run is seeded (`--seed`, default 0) and two runs with one seed give the same tables.

## Models

Only the four models the grouping work is about: `brenda_merged`, `e_coli_core`, `R-HSA-77289`,
`skm_pss`. The two small ones (`skm_ja`, `skm_sa`) stay in `../examples/`.
