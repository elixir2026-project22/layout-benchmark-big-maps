# layout-benchmark-big-maps

Big pathway maps in SBML, drawn in two stages: first **grouping** - what belongs together - then
**joining** - where everything goes. From [BioHackathon Europe 2026](https://biohackathon-europe.org/),
project 22, *Open Pathway Maps: Scalable and FAIR Visualization of Metabolic Networks*.

A map of a thousand nodes is unreadable as one flat drawing, and a layout algorithm alone does not make it
readable. This repository splits the problem:

1. **Grouping.** Eight strategies, each a program that reads nothing but the SBML, turn a model into a
   nested grouping of all its species and reactions. One written rulebook of 27 rules scores them.
   `CONTRACT.md` defines exactly what a strategy reads and writes, so the strategies differ in grouping
   logic only.
2. **Joining (layout).** The grouping says what belongs together and nothing about position. The layout
   stage takes a grouping and one layout algorithm from
   [layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark)'s registry and runs **the
   same algorithm at every level** - inside a leaf group, between the sub-groups of a group, and between
   the top-level groups. `layout/PROTOCOL.md` states the five rules of that placement.

Both stages are shown as galleries of rendered SVG, committed here: 52 pictures for the grouping stage
and 288 (8 options x 9 groupings x 4 models) for the layout stage. They are online at
**<https://elixir2026-project22.github.io/layout-benchmark-big-maps/>** ([grouping](https://elixir2026-project22.github.io/layout-benchmark-big-maps/grouping/index.html),
[layout](https://elixir2026-project22.github.io/layout-benchmark-big-maps/layout/index.html)); see
[Website](#website).

> Status: the grouping comparison is finished and reported in [`doc/REPORT.md`](doc/REPORT.md); the layout
> stage is complete for the four big models and has no combined score yet - the gallery reports canvas
> size, density and time per picture. The models themselves are not in this repository; see
> [Not in this repository](#not-in-this-repository).

## Input

The work targets four models, with two small Stress Knowledge Map pathways carried along in the grouping
stage. `panel/panel.tsv` is the provenance record: for each model the source path it was taken from, its
SHA-256, and what was done to it.

| model | species | reactions | compartments | source | in the layout stage |
|---|---|---|---|---|---|
| `brenda_merged` | 178 | 218 | 1 | merged BRENDA pathways (groups and model-level notes removed: both name the source pathways) | yes |
| `e_coli_core` | 72 | 95 | 2 | [BiGG Models](http://bigg.ucsd.edu/) | yes |
| `R-HSA-77289` | 101 | 47 | 4 | [Reactome](https://reactome.org/content/detail/R-HSA-77289) | yes |
| `skm_pss` | 961 | 856 | 10 | [Stress Knowledge Map](https://skm.nib.si), PSS | yes |
| `skm_ja` | 117 | 82 | 5 | Stress Knowledge Map, JA pathway | no |
| `skm_sa` | 193 | 159 | 6 | Stress Knowledge Map, SA pathway | no |

Every strategy reads the same three tables per model, written by `common/sbml_features.py` into
`panel/features/<model>/`: `nodes.tsv` (every species and reaction with its compartment, SBO term and
ancestors, EC, GO, cross-references, notes, degree and its place in the cell), `edges.tsv` (one row per
participant with its role and stoichiometry) and `compartments.tsv` (the compartment containment tree).
`big-maps prepare` writes them; see [Reproduce it](#reproduce-it).

## Stage 1: grouping

Each strategy writes `groups.tsv` (every node, its group at levels 1-3, the rule that placed it, clones)
and `groups_meta.tsv` (every group: level, parent, label, basis, where the name came from). The format,
and the seven rules every strategy keeps - one program for all models, deterministic, from the SBML only,
names from identifiers rather than from words, localisation given rather than chosen, every node placed, a
declared fallback - are in [`CONTRACT.md`](CONTRACT.md).

| strategy | top level | inside | rules |
|---|---|---|---|
| `type_location` | reaction / signalling type from the SBO hierarchy, EC class where SBO does not divide | localisation per compartment, with boundary groups | [RULES.md](strategies/type_location/RULES.md) |
| `location_type` | localisation per compartment | reaction type, EC class as fallback | [RULES.md](strategies/location_type/RULES.md) |
| `ontology` | annotation trees: notes bins (MapMan), EC, pathway cross-references, GO; SBO fallback | the next level of the same tree, then localisation | [RULES.md](strategies/ontology/RULES.md) |
| `topology_leiden` | Leiden communities of the species-reaction graph (modularity) | the same algorithm again inside oversized groups | [RULES.md](strategies/topology/RULES.md) |
| `topology_infomap` | Infomap communities (flow) | as above | [RULES.md](strategies/topology/RULES.md) |
| `topology_walktrap` | Walktrap communities (random walk) | as above | [RULES.md](strategies/topology/RULES.md) |
| `cascade` | flow anchors (signalling), routes and cycles (metabolism) | stage, then localisation | [RULES.md](strategies/cascade/RULES.md) |
| `hybrid` | function term where annotated, connectivity where not, size-bounded | finer function, reaction type, localisation | [RULES.md](strategies/hybrid/RULES.md) |
| `baseline` | the project's current grouping (`pathway-grouping`), converted to this format | flat | [baseline/README.md](baseline/README.md) |

The scorer (`scoring/score_grouping.py`) implements the merged rulebook in [`scoring/RULES.md`](scoring/RULES.md),
which was written from three independent rule proposals (`rules/rules_graph.md`, `rules/rules_biology.md`,
`rules/rules_readability.md`, each with its metrics module and tests). Results:
[`doc/REPORT.md`](doc/REPORT.md), `scoring/scores.csv` (per model and rule) and `scoring/scores_summary.csv`.
`baseline` is last on every model; ranks 1-2 (`hybrid`, `location_type`) are a tie the evidence does not
separate.

## Stage 2: joining the groups

One option per runnable family of `layout_benchmark.layouts.REGISTRY`, in the registry's own order. The
median is over the 36 pictures of that option in the committed gallery; *area per node* is the canvas area
divided by the number of drawn nodes, so lower is denser.

| option | family | median px² per node | median s |
|---|---|---|---|
| [`graphviz-dot`](https://graphviz.org/docs/layouts/dot/) | layered | 52,450 | 6.3 |
| [`igraph-sugiyama`](https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_sugiyama) | layered (stands in for *tree*) | 76,963 | 3.2 |
| [`graphviz-neato`](https://graphviz.org/docs/layouts/neato/) | stress | 35,224 | 8.7 |
| [`graphviz-sfdp`](https://graphviz.org/docs/layouts/sfdp/) | force-directed | 40,420 | 10.1 |
| [`igraph-drl`](https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_drl) | for large graphs | 52,765 | 17.2 |
| [`graphviz-twopi`](https://graphviz.org/docs/layouts/twopi/) | radial | 38,766 | 13.6 |
| [`graphviz-circo`](https://graphviz.org/docs/layouts/circo/) | circular | 38,633 | 10.3 |
| [`igraph-mds`](https://python.igraph.org/en/stable/api/igraph.Graph.html#layout_mds) | baseline | 92,603 | 6.2 |

`big-maps options` prints this list and says which tools are present on the machine. The registry's other
families are not offered because they need software this repository cannot assume: every cytoscape.js /
ELK / graphology / Reagraph option needs Node, and the pathway-specific tools (`rules-claude`, `ecmapper`,
`metacarto`, `sbml2map`) and `sbmlnetwork` need their own checkouts under `.tools/`. Adding one is a new
entry in `OPTIONS` in `layout/options.py` and nothing else.

What the placement does, in short (the full statement is [`layout/PROTOCOL.md`](layout/PROTOCOL.md)):

1. the grouping is an input and is never changed, with one exception;
2. that exception: a box holding exactly one species or reaction and no sub-group is **not a group** - its
   node moves into the parent and the box is dropped. A compartment box is never dropped, because
   localisation is given rather than chosen. In this panel that drops groups only in `baseline` (1 on
   `e_coli_core`, 2 on `R-HSA-77289`, 167 on `skm_pss`);
3. inside a box, its own species and reactions and its already-sized sub-boxes form one item graph - two
   items are joined when the model has an edge between any member of one and any member of the other,
   weighted by how many - and the option places that graph;
4. the top level is the same call on the top-level groups, so groups are distributed by the same logic as
   the elements inside them;
5. components are laid out separately and packed in rows, then items are pushed apart and the positions
   scaled until nothing overlaps. Scaling positions while sizes stay fixed cannot bring two items closer,
   so this terminates and creates no new overlap.

A run writes `layout/out/<option>/<strategy>/<model>/positions.tsv` (every drawn node: id, clone, box,
centre, size) and `boxes.tsv` (every box: kind, basis, label, level, parent, corner, size, members).

### Determinism

One seed gives one drawing, byte for byte. Every level sorts its items and edges by a stable name before
the tool is called, never by object identity; Graphviz's random starts are seeded per call; igraph's global
random generator is seeded under a lock, so a thread pool cannot interleave two layouts. `--seed` changes
the drawing and defaults to 0. The grouping stage is deterministic by the same contract (rule 2), and the
scorer's rule M27 is set by a byte comparison of two runs.

## Setup

Requirements: Python 3.11 or newer, and the **Graphviz binaries on PATH** (`dot`, `neato`, `sfdp`, `twopi`,
`circo`) - they are called as subprocesses with `-Tplain`, so `pygraphviz` is not needed. Node is not
needed.

```bash
git clone https://github.com/elixir2026-project22/layout-benchmark-big-maps.git
cd layout-benchmark-big-maps

# with uv (as in layout-benchmark; uv.lock pins every version)
uv sync --extra tests
uv run pytest

# or with pip
python -m venv .venv && . .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e ".[tests]"
pytest
```

Graphviz itself: `sudo apt-get install graphviz` (Debian/Ubuntu), `brew install graphviz` (macOS),
`winget install graphviz` (Windows). Check with `big-maps options`: every option must say `yes`.

## Reproduce it

Everything below runs from the repository root in a clean clone. Only step 1 needs files that are not in
this repository.

```bash
# 1. the models. panel/panel.tsv records where each came from and its SHA-256.
#    --strip-groups removes groups a file already carries (brenda_merged needs it).
big-maps prepare /path/to/e_coli_core.sbml --name e_coli_core
big-maps prepare /path/to/brenda_merged.xml --name brenda_merged --strip-groups
big-maps prepare /path/to/R-HSA-77289.sbml --name R-HSA-77289
big-maps prepare /path/to/model-with-model-fixes.xml --name skm_pss
big-maps models                      # the panel and the SHA-256 of each source
# the same by hand, if you prefer the scripts:
#   python common/contract.py strip-groups SRC.xml -o panel/models/<model>.xml
#   python common/sbml_features.py panel/models/<model>.xml -o panel/features/<model>

# 2. run a grouping strategy on a model (identical command line for every strategy)
python strategies/hybrid/strategy.py --features panel/features/skm_pss \
    --sbml panel/models/skm_pss.xml --out strategies/hybrid/out/skm_pss
# the three topology variants differ only by --algorithm
python strategies/topology/strategy.py --features panel/features/skm_pss \
    --sbml panel/models/skm_pss.xml --algorithm leiden --seed 1 \
    --out strategies/topology/topology_leiden/out/skm_pss \
    --stats strategies/topology/topology_leiden/stats/skm_pss.json

# 3. check it against the contract, and export it as SBML Level 3 groups
python common/contract.py validate panel/features/skm_pss strategies/hybrid/out/skm_pss   # prints: valid
python common/contract.py write-sbml panel/models/skm_pss.xml strategies/hybrid/out/skm_pss \
    --features panel/features/skm_pss --strategy hybrid \
    -o strategies/hybrid/out/skm_pss/skm_pss_grouped.xml

# 4. score every strategy on every model
python scoring/collect_groupings.py          # strategies/*/out/* -> scoring/groupings/
python scoring/score_grouping.py --all --panel panel/features --strategies scoring/groupings \
    --baseline baseline --heldout heldout --out scoring/scores.csv --summary scoring/scores_summary.csv

# 5. the grouping gallery: every model under every strategy
python examples/build_gallery.py             # -> examples/svg/, gallery.tsv, index.html
python examples/build_gallery.py --no-render # the page only, from the committed gallery.tsv

# 6. the layout gallery: 8 options x 9 groupings x 4 models
python layout/options.py                                   # which options can run here
python layout/placement.py panel/features/skm_pss strategies/hybrid/out/skm_pss \
    --option graphviz-neato -o layout/out/graphviz-neato/hybrid/skm_pss
python layout/render_layout.py panel/features/skm_pss strategies/hybrid/out/skm_pss \
    --option graphviz-neato -o layout/examples/svg/graphviz-neato__hybrid__skm_pss.svg \
    --title "graphviz-neato — hybrid — skm_pss"
python layout/build_gallery.py --jobs 6                    # the whole gallery
python layout/build_gallery.py --no-render                 # the page only, from the committed gallery.tsv
python layout/build_gallery.py --option graphviz-dot --strategy hybrid --model skm_pss
```

Steps 5 and 6 need step 1 (the renderer reads `panel/features/<model>/`); `--no-render` does not, which is
why it is the end-to-end step CI runs. The committed gallery took 8,034 s of layout time in total, built
with `--jobs 6` (threads, not processes: the work is in Graphviz subprocesses and in igraph/numpy, which
release the GIL).

## What the example outputs show

| output | what it is |
|---|---|
| [`examples/index.html`](examples/index.html) | the grouping gallery: 52 pictures, every model under every strategy, with levels, group count and score per cell |
| `examples/gallery.tsv` | one row per picture: nodes, groups per level, named groups, largest level-1 group, singleton leaves, clones, score |
| [`layout/examples/index.html`](layout/examples/index.html) | the layout gallery: 288 pictures, every option x grouping x model |
| `layout/examples/gallery.tsv` | one row per picture: option, family, grouping, model, nodes, boxes, ungrouped singletons, canvas, aspect, area per node, seconds, error |
| `scoring/scores.csv`, `scoring/scores_summary.csv` | the rulebook score per model and per rule, and the ranking |
| `scoring/sensitivity/` | how the ranking moves under 11 weightings and with each rule dropped in turn |

Both galleries were built with 0 failures. The SVGs are the bulk of the repository: `layout/examples/svg/`
is 87.3 MB over 288 files (largest 0.83 MB, `skm_pss` accounts for 55.0 MB) and `examples/svg/` is 12.6 MB
over 52 files. They compress about 9-fold: the whole repository is 103.6 MB over 547 files in the working
tree but packs into 9.8 MiB, which is what a clone downloads.

## Website

Both galleries are published as a static site, so they can be looked at without cloning:

**<https://elixir2026-project22.github.io/layout-benchmark-big-maps/>**

| page | what it is |
|---|---|
| [`/grouping/`](https://elixir2026-project22.github.io/layout-benchmark-big-maps/grouping/index.html) | the grouping gallery, 52 pictures |
| [`/layout/`](https://elixir2026-project22.github.io/layout-benchmark-big-maps/layout/index.html) | the layout gallery, 288 pictures |

The site is the two galleries copied unchanged behind one landing page - nothing is re-rendered, so it
always shows exactly what the repository holds. It is published the way
[layout-benchmark](https://elixir2026-project22.github.io/layout-benchmark/) publishes its report: an
orphan `gh-pages` branch with `.nojekyll` at the root, which GitHub Pages serves as it is.

```bash
# assemble the site from the committed galleries (writes site/, which is gitignored)
python scripts/build_site.py

# publish it: build it inside a worktree on the gh-pages branch, commit, push
git worktree add --orphan -b gh-pages ../_gh-pages      # first time only
python scripts/build_site.py -o ../_gh-pages
git -C ../_gh-pages add -A
git -C ../_gh-pages commit -m "Site: galleries at $(git rev-parse --short main)"
git -C ../_gh-pages push -u origin gh-pages
```

Afterwards the worktree is reusable: `python scripts/build_site.py -o ../_gh-pages`, commit, push.
Drop it with `git worktree remove ../_gh-pages`; the branch stays.

In the repository's **Settings → Pages**, set *Source* to *Deploy from a branch*, branch `gh-pages`,
folder `/ (root)`. The site is about 100 MB over 346 files, well inside the 1 GB Pages limit; the
`skm_pss` pictures are the large ones, and every gallery page loads its thumbnails lazily.

## Not in this repository

The data is left out; the code, the documentation and the rendered results are kept. Everything left out is
regenerated by the commands above, and `.gitignore` lists each path so it cannot be committed by accident.

| left out | size | how to get it back |
|---|---|---|
| `panel/models/*.xml` | 5.8 MB, 6 files | your own copies of the models in `panel/panel.tsv`, then `big-maps prepare` (step 1) |
| `panel/features/**` | 2.2 MB, 18 files | `big-maps prepare` (step 1) |
| `heldout/**` | 66.2 MB | curated reference labels, used for scoring only, and a 50 MB Human-GEM copy; not public |
| `strategies/*/out/*/*_grouped.xml` | 49.5 MB, 49 files | `python common/contract.py write-sbml ...` (step 3) |
| `layout/out/**` | 17.2 MB, 576 files | `python layout/placement.py ...` (step 6) |
| `scoring/groupings/**` | 1.1 MB | `python scoring/collect_groupings.py` (step 4); it is a copy of `strategies/*/out/` |
| `scoring/heldout_round1/**` | 2.2 MB | a scoring round on the held-out model, which is not public |

`strategies/*/out/<model>/groups.tsv` and `groups_meta.tsv` **are** kept (1.2 MB over 49 groupings): they
are small, they are what the layout stage consumes, and keeping them makes the layout gallery reproducible
without rerunning every strategy first.

## Tests

`pytest` runs three sets, configured in `pyproject.toml`:

- `tests/` - the packaging and the layout stage, on a four-species model the tests build themselves with
  libSBML: the option registry, contract validity, placement with each of the eight options, no sibling
  box overlapping another, one seed giving one placement, the singleton rule, and that
  `layout/examples/index.html` is byte-for-byte what `layout/build_gallery.py --no-render` writes from the
  committed `gallery.tsv`;
- `rules/test_metrics_*.py` - the three rule proposals' metrics;
- `scoring/test_score_grouping.py` - the scorer.

The tests that need `panel/features/` skip when it is absent, which is the state of a clean clone: 23
pass in `tests/`, 102 pass and 104 skip in `rules/` and `scoring/`. Run step 1 and the 104 run too.

## Adding a strategy, an option or a rule

- **A grouping strategy**: a folder `strategies/<name>/` with `strategy.py` (the command line of
  `CONTRACT.md` rule 1), `RULES.md` (the numbered rules its `rule` column refers to) and `NOTES.md`. Then
  `python common/contract.py validate ...`, an entry in `STRATEGIES` in `examples/build_gallery.py` and
  `layout/build_gallery.py`, and in `SOURCES` in `scoring/collect_groupings.py`.
- **A layout option**: one entry in `OPTIONS` in `layout/options.py` - a function
  `(items, sizes, edges, seed) -> positions`, its family and a documentation link. Nothing else changes.
- **A rule**: a metric function in one of `rules/metrics_*.py` with its test, then the rule in
  `scoring/RULES.md` and its implementation in `scoring/score_grouping.py`.

## Licence

MIT (code), the licence of
[layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark), whose conventions this
repository follows; see [LICENSE.md](LICENSE.md) and [COPYING.md](COPYING.md). The repository the code was
assembled from carried no licence file of its own, so this is the project's licence applied to it.

The input models are not in this repository and keep their own licences: the Stress Knowledge Map data are
CC BY-NC 4.0 (<https://skm.nib.si>); BRENDA's pathway map is CC BY 4.0
(<https://www.brenda-enzymes.org/copy.php>); BiGG Models (`e_coli_core`) is © 2019 The Regents of the
University of California, free for educational, research and non-profit use (<https://bigg.ucsd.edu/>);
`R-HSA-77289` is from [Reactome](https://reactome.org). `common/SBO_OBO.obo` is the SBO release of
28 August 2021 from <https://github.com/EBI-BioModels/SBO>; the GO cellular-component ids in
`common/go_cellular_component.json` were looked up in the EMBL-EBI Ontology Lookup Service in November 2026.
