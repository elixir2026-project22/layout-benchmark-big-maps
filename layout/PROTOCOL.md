# Protocol: from a grouping to a drawn map

What turns `groups.tsv` into coordinates. Four steps, each one a rule that can be checked against the
tables the run writes (`out/<option>/<strategy>/<model>/positions.tsv` and `boxes.tsv`).

## 1. The grouping is an input, never changed

`../strategies/<name>/out/<model>/groups.tsv` and `groups_meta.tsv` are read exactly as the contract
(`../CONTRACT.md`) defines them, through the same reader the grouping gallery uses
(`../examples/render_grouping.py`, function `build`). That reader gives a tree of boxes:

- a group of `groups_meta.tsv` is a box, nested by its `parent`;
- wherever a level's basis is `localisation`, its boxes are replaced by the compartment containment
  tree of `compartments.tsv` - extracellular around cytosol, organelles inside cytosol - so the
  nesting of compartments is never something a layout decides;
- a node with `loc_kind` = `boundary` goes into a band on the border of its space;
- `clone_in` puts a further copy of a species in each group named there.

No node is added, removed or moved between groups by this module, with the single exception of step 2.

## 2. A group of one entity is not a group

A box that holds exactly one species or reaction and has no sub-group and no band is dropped; its
node is placed in the parent instead. The count is reported per picture as *ungrouped singletons*
and the dropped ids are available from `Placement.dropped_singletons`.

A compartment box is never dropped, even when one species is alone in it: `../CONTRACT.md` rule 5
says localisation is given rather than chosen, and a compartment with one member is still that
compartment. In the panel this affects `type_location` (1, 1 and 12 compartments on `e_coli_core`,
`R-HSA-77289` and `skm_pss`), `location_type` (1 on `skm_pss`) and `hybrid` (12 on `skm_pss`); the
groups actually dropped are `baseline`'s (1 on `e_coli_core`, 2 on `R-HSA-77289`, 167 on `skm_pss`).

## 3. The elements within a group are distributed by the layout option

Inside one box, the items of that level are its own species and reactions together with its
sub-boxes, each sub-box already laid out and therefore of known size. Two items are joined when the
model has an edge between any member of the one and any member of the other; the weight is how many
such edges there are, and the direction is the model's (reactant and modifier into the reaction,
reaction out to the product), which only the layered options use. The option places that graph.

Edges that end in a boundary band are left out of its level: a band is pinned to the border of its
space, so no layout decides where it goes.

## 4. The groups are distributed by the same logic

The top level is the same call on the same kind of graph, where the items are the top-level groups.
One option therefore means one algorithm, run at every level of the grouping - inside a leaf group,
between the sub-groups of a group, and between the top-level groups.

## 5. Sizes are respected, and no box overlaps another

Layout tools differ in whether they know node sizes, so every level is finished the same way:

1. the level's item graph is split into connected components, each component laid out on its own and
   the components packed in rows (Graphviz does the same for its own layouts with `pack=true`);
   without this, an algorithm built for large graphs throws unconnected pieces far apart and is
   judged on that rather than on how it draws a connected piece;
2. the component is scaled to a readable density: nearest neighbours about `SPACING` node widths
   apart, or, when that number is meaningless because the tool returned coincident points, the scale
   that makes the items cover `FILL` of the bounding box - a layout with no extent at all is packed
   in rows instead;
3. items are pushed apart along their smaller overlap, and whatever still overlaps is resolved by one
   uniform scale of the positions. Scaling positions while the sizes stay fixed can never bring two
   items closer, so this terminates and creates no new overlap.

The box's own size follows from the result: the bounding box of its items, plus padding, its title
bar, and its band.

## Determinism

One seed gives one drawing. Every level sorts its items and its edges by a stable name before the
tool is called (never by object identity), Graphviz's random starts are seeded per call, and igraph's
global random generator is seeded under a lock so that the thread pool cannot interleave two layouts.
`--seed` changes the drawing; the default is 0.

## What a run writes

| file | one row per |
|---|---|
| `out/<option>/<strategy>/<model>/positions.tsv` | drawn node: its id, whether it is a clone, its box, its centre and size |
| `out/<option>/<strategy>/<model>/boxes.tsv` | box: its id, kind, basis, label, level, parent, corner, size, member count |
| `examples/gallery.tsv` | picture: the option, its family, the grouping, the model, the canvas, the density and the time |
| `examples/svg/<option>__<strategy>__<model>.svg` | picture |
| `examples/index.html` | the gallery |
