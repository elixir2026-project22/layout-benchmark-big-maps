# Readability rules for a nested grouping of an SBML model

Viewpoint: **what a human reader can take in** from a map drawn with one box per group.
Reference implementation: `metrics_readability.py` (`metrics(nodes, edges, comps, groups, meta, reference=None)`),
tests: `test_metrics_readability.py`. Written blind: no strategy output was seen or used.

```
python metrics_readability.py panel/features/MODEL strategies/NAME/out/MODEL            # metrics, verdicts, score
python metrics_readability.py panel/features/MODEL strategies/NAME/out/MODEL --json
python -m pytest test_metrics_readability.py -q                                        # 48 tests
```

## Notation (all from the contract tables)

- `N` nodes of `nodes.tsv` that are placed in `groups.tsv`; `E` rows of `edges.tsv` (one drawn arc each).
- `path(v)` = the non-empty `L1..L3` of node `v`; `home(v)` = its deepest group.
- `size(g)` = number of nodes with `g` on their path; `load(g)` = number of nodes whose home is `g`.
- **box** = a group as drawn. A space group `c` and a sibling `c@boundary` group (basis `localisation`, same
  parent) count as **one** box wherever boxes are counted: the membrane is the border of the same compartment.
- `boxes(p)` = number of distinct child boxes of parent `p` (+1 if `p` also holds loose nodes next to child
  groups); `boxes(map)` = number of distinct L1 boxes.
- `copies(s)` = `path(s)` plus the path of every group named in `clone_in` of species `s`.
- **Localisation levels are a backdrop, not a choice.** The contract fixes where a node is and draws
  compartments nested by `compartments.tsv`. Therefore lower limits on counts (RD01, RD08), the singleton
  rule (RD03) and depth (RD12) do not charge a strategy for what the compartments of the model impose;
  upper limits and every other rule apply to localisation groups like to any other group.
- NaN = the rule does not apply to this *model*. Where the *grouping* empties a denominator in a way that
  would hide a defect (no group is named), the metric takes its worst value, so omission is never rewarded.

Verdict and score of every rule come from three numbers: **pass** limit (score 1), **warn** limit (score
0.5), **zero** point (score 0), linear in between. `fail` = beyond the warn limit.

| id | metric | better | pass | warn | zero | weight |
|---|---|---|---|---|---|---|
| RD01_TOP | `top_groups` | inside range | 3..12 | 2 or 13..20 | 1 or >= 40 | 3 |
| RD02_BOXSIZE | `oversize_node_share` | lower | <= 0.10 | <= 0.30 | 0.70 | 3 |
| RD03_SINGLETON | `singleton_share` | lower | <= 0.10 | <= 0.25 | 0.50 | 3 |
| RD04_NAMED | `named_share` | higher | >= 0.90 | >= 0.60 | 0.20 | 3 |
| RD05_LONGEDGE | `long_edge_share` | lower | <= 0.35 | <= 0.60 | 0.90 | 2 |
| RD06_CLONE | `cloned_species_share` | lower | <= 0.10 | <= 0.25 | 0.50 | 2 |
| RD07_DISTINCT | `distinct_share` | higher | >= 0.95 | >= 0.80 | 0.50 | 2 |
| RD08_FAN | `fan_ok_share` | higher | >= 0.90 | >= 0.70 | 0.40 | 2 |
| RD09_BACKED | `backed_share` | higher | >= 0.80 | >= 0.50 | 0.20 | 1.5 |
| RD10_COMPDRAW | `comp_redraw` | lower | <= 1.5 | <= 3.0 | 6.0 | 1.5 |
| RD11_BALANCE | `size_balance` | higher | >= 0.50 | >= 0.25 | 0.10 | 1 |
| RD12_DEPTH | `depth_excess` | lower | 0 | <= 1.0 | 2.0 | 1 |

Sources are numbered [1]..[12] and listed at the end. "judgement" marks a number that is mine and has no
empirical source.

---

## The rules, most important first

### RD01_TOP - a handful of boxes at the top level
- **By eye:** the map shows between 3 and 12 top-level boxes.
- **Metric:** `top_groups = boxes(map)` = number of distinct L1 boxes (`c` and `c@boundary` merged).
- **Direction, range:** inside the range is better; integer >= 1.
- **Thresholds:** pass 3..12, warn 2 or 13..20, fail 1 or > 20.
  Source: the span of immediate memory is about seven items, 5 to 9 [1]; without rehearsal the limit is
  about four chunks, 3 to 5 [2] - hence the lower end 3. A map is scanned, not memorised, and in menu
  hierarchies more depth hurt search while a medium breadth/depth did best [3] - hence an upper limit
  somewhat above nine. The exact numbers 12 and 20 are **judgement**.
- **Not applicable:** `N <= 20` (a graph of about twenty nodes reads well as one node-link picture [4]).
  The lower limit is waived when every L1 group has basis `localisation`: the number of compartments is
  given by the model (flag `top_all_localisation`); the upper limit still holds.
- **Gaming:** (a) splitting everything into exactly 3..12 arbitrary groups passes - counter-rules RD05
  (edges) and RD04/RD07 (names). (b) 3..12 huge boxes - counter-rule RD02.

### RD02_BOXSIZE - no box is a hairball
- **By eye:** no box holds more than about 50 loose nodes (species + reactions) that are not in a sub-box.
- **Metric:** `oversize_node_share = sum(load(g) for load(g) > 50) / N`.
  Reported with it: `box_size_max`, `box_size_median`, `box_over100_node_share`, `largest_box_share`.
- **Direction, range:** lower is better; [0, 1].
- **Thresholds:** pass <= 0.10, warn <= 0.30, fail > 0.30. The box limit 50 has a source: readers had
  significant difficulty in dense node-link diagrams above 50 nodes and in sparse ones above 100 [5];
  above twenty vertices node-link diagrams already lose to matrices on most tasks in dense graphs [4].
  Pathway graphs are sparse (1.1 to 2.3 arcs per node in the panel), so 50 is a middle value and 100 is
  reported as the hard limit. The share cut-offs 0.10 / 0.30 are **judgement**.
- **Not applicable:** never (a model of <= 50 nodes passes trivially).
- **Gaming:** one node per box scores 0 - counter-rules RD03 and RD01. Cutting a large box into arbitrary
  chunks of 50 - counter-rule RD05.

### RD03_SINGLETON - boxes group something
- **By eye:** at most one box in ten holds a single node.
- **Metric:** `singleton_share = |{g : size(g) = 1, not forced}| / |groups|`. A singleton is *forced*, and
  not counted, when it has basis `localisation` and is the only group declared for its `compartment_id`
  (a compartment that holds one node). Reported with it: `singleton_node_share`.
- **Direction, range:** lower is better; [0, 1].
- **Thresholds:** pass <= 0.10, warn <= 0.25, fail > 0.25 - **judgement**. Anchor: the grouping that
  motivated this work has 167 singletons in 346 groups (0.48) and is described as unclear; it fails.
- **Not applicable:** never.
- **Gaming:** sweeping all leftovers into one unnamed bucket removes the singletons - counter-rules RD02
  (bucket size), RD04 (the bucket has no name), RD05 (its edges are long). Declaring leftovers as
  `localisation` boxes does not help: the exemption needs a unique `compartment_id`, and repeated boxes of
  one compartment are charged by RD10.

### RD04_NAMED - every box has a title
- **By eye:** nine boxes in ten carry a label.
- **Metric:** `named_share = |{g : label(g) != ""}| / |groups|`. Reported with it: `named_node_share`
  (share of nodes whose home box is labelled), `member_named_share`.
- **Direction, range:** higher is better; [0, 1].
- **Thresholds:** pass >= 0.90, warn >= 0.60, fail < 0.60. Source for the principle: compartments and
  other graphical annotations organise a pathway layout and should be labelled for identification [10].
  The cut-offs are **judgement**; anchor: 34 named groups of 346 (0.10) is the unclear baseline.
- **Not applicable:** never (the contract allows a `member` name for any group).
- **Gaming:** naming every box after one member gives 1.0 - counter-rules RD09 (a member name is not
  identifier-backed), RD03 (if the boxes are singletons). Giving every box the same coarse term gives 1.0 -
  counter-rule RD07.

### RD05_LONGEDGE - arrows stay inside a readable neighbourhood
- **By eye:** for at least two arrows in three, both ends (using the nearest clone of the species) lie in
  a common box of at most 100 nodes.
- **Metric:** an arc `(r, s)` is *long* when, for every `p` in `copies(s)`, `path(r)` and `p` share no
  group, or the deepest group they share has `size > 100`.
  `long_edge_share = |long arcs| / |E|`. Reported with it: `long_edge_share_raw` (clones ignored),
  `l1_cross_share` (ends in different L1 groups), `cut_edge_share` (ends not in the same box),
  `boundaries_crossed_mean`, `modularity_L1` (Newman-Girvan Q of the L1 partition).
- **Direction, range:** lower is better; [0, 1].
- **Thresholds:** pass <= 0.35, warn <= 0.60, fail > 0.60. Reasoning: arcs should be kept short and
  crossings limited [9]; edge crossings are by far the aesthetic with the largest effect on understanding
  [6], and arcs between distant boxes are the ones that cross. The span 100 is the size above which path
  following became hard even in sparse diagrams [5]. For a flat partition
  `Q = (1 - long share) - sum_g a_g^2`; networks with real community structure are reported to have Q
  between about 0.3 and 0.7 [7]. With about ten balanced groups (`sum a_g^2` about 0.1) Q = 0.3 is a long
  share of 0.60 and Q = 0.55 a long share of 0.35. The mapping from Q to the two limits is **judgement**.
- **Not applicable:** a model without edges (NaN).
- **Gaming:** "everything in one group" has no arc that leaves a box (`l1_cross_share` = 0) - this is why
  the rule uses the span and not the boundary: in a single box of more than 100 nodes every arc is long,
  and RD01/RD02 fail as well. Cloning every species next to every reaction removes all long arcs -
  counter-rule RD06, which carries the same weight so that the trade never pays (tested).

### RD06_CLONE - only a few hub species are drawn more than once
- **By eye:** at most one species in ten appears with a clone marker.
- **Metric:** `cloned_species_share = |{s : clone_in(s) names >= 1 known group}| / |species|`.
  Reported with it: `clone_load` (clone glyphs per species), `clone_useless_share` (clones placed in a
  group that holds none of the species' reactions; should be 0), and `clone_demand` = number of pairs
  (species, home box of a reaction it is long-connected to) / |species| - the clones the grouping would
  need to have no long arc, whether or not the strategy declared any.
- **Direction, range:** lower is better; [0, 1].
- **Thresholds:** pass <= 0.10, warn <= 0.25, fail > 0.25 - **judgement**. Sources for the principle: a
  layout should avoid multiple copies of a node unless they carry meaning, such as transport or common
  reactants like ATP/ADP [10]; clone markers are meant for duplicating hub nodes to limit crossing arcs
  [9]; SBGN allows duplicates only when marked with a clone marker [11]; node duplication helped
  community tasks but can interfere with other readability tasks [8]. The rule therefore counts *how many
  different species* are cloned (hubs are few), not how many copies one hub has.
- **Not applicable:** a model without species (NaN).
- **Gaming:** declaring no clone at all scores 0 - counter-rule RD05 (long arcs stay). Cloning one hub into
  every box is cheap here by design; `clone_load` shows it.

### RD07_DISTINCT - a title tells this box from its neighbours
- **By eye:** no two sibling boxes share a title, no box repeats the title of the box around it, and a box
  called X holds at least half of the X that its parent box contains.
- **Metric:** over named groups, `distinct_share` = share of groups that are (a) the only one with their
  label among their siblings (labels compared case- and whitespace-insensitively; `c` and `c@boundary` of
  one compartment do not count as duplicates), (b) not labelled like any ancestor, and (c) have
  `recall >= 0.5`, where `recall(g) = |carriers(x) in g| / |carriers(x) in parent(g)|` (parent = whole map
  for L1) and `carriers(x)` = nodes to which the identifier `name_xref` applies: SBO term on the node's
  `sbo_path`, EC number inside the EC class, GO id in `go`, the cross-reference (or a dotted child of it)
  in `xrefs`, `loc_space` = the compartment. (c) is skipped for `member` and `notes` names.
  Reported with it: `sibling_dup_share`, `same_as_ancestor_share`, `name_low_recall_share`,
  `name_recall_mean`, `global_dup_share`, `label_len_median`.
- **Direction, range:** higher is better; [0, 1]; 0 when no group is named.
- **Thresholds:** pass >= 0.95, warn >= 0.80, fail < 0.80 - **judgement**. Anchor: "EC 1" used twice among
  siblings fails (a) for both groups, and each of them holds at most half of the EC 1 reactions.
- **Not applicable:** never.
- **Gaming:** numbering the labels ("EC 1 (a)", "EC 1 (b)") passes (a) but still fails (c). Naming every
  box after a member passes - counter-rule RD09. Leaving boxes unnamed to avoid the check - RD04.

### RD08_FAN - every opened box shows a handful of sub-boxes
- **By eye:** a box that has sub-boxes has between 2 and 12 of them.
- **Metric:** over groups with at least one child group, `fan_ok_share` = share with
  `2 <= boxes(g) <= 12`. A parent whose children are all `localisation` groups has no lower limit.
  Reported with it: `fan_max`, `only_child_share`.
- **Direction, range:** higher is better; [0, 1].
- **Thresholds:** pass >= 0.90, warn >= 0.70, fail < 0.70. The limits 2..12 follow RD01 ([1], [2], [3]);
  an only child is a border without information. The share cut-offs are **judgement**.
- **Not applicable:** a flat grouping (no group has a child): NaN.
- **Gaming:** a flat grouping avoids the rule - counter-rules RD01 and RD02 then carry the load.

### RD09_BACKED - a title can be looked up
- **By eye:** four titles in five come from an ontology term, an EC class, a cross-reference or a
  compartment that the box's own members carry.
- **Metric:** over named groups, `backed_share` = share with `name_source` in {sbo, ec, go, xref,
  compartment}, a non-empty `name_xref`, and at least one member node in `carriers(name_xref)`.
  `member` and `notes` names are not counted as backed.
- **Direction, range:** higher is better; [0, 1]; 0 when no group is named.
- **Thresholds:** pass >= 0.80, warn >= 0.50, fail < 0.50 - **judgement**. Source for the principle:
  standardised names and identifiers are essential for users and tools to interpret a pathway model, in
  line with FAIR [10].
- **Not applicable:** the model carries no SBO term, EC number, GO id or cross-reference on any node and
  has a single compartment (NaN).
- **Gaming:** attaching a syntactically valid identifier that no member carries is not counted. Attaching
  the root-like term every node carries (e.g. the SBO term of all reactions) is backed but fails RD07.

### RD10_COMPDRAW - each compartment is drawn once
- **By eye:** a compartment appears as one region of the map, not as a patch inside many boxes.
- **Metric:** for node `v`, `container(v)` = the deepest group on `path(v)` whose members lie in more than
  one `loc_space`, or the map itself if there is none. `pieces(c)` = number of distinct containers of the
  nodes with `loc_space = c`. `comp_redraw = mean over compartments c with nodes of pieces(c)`.
  It needs no declared localisation group: a box that is pure in one compartment can be drawn inside that
  compartment, a box that mixes compartments needs its own patch of each. Reported with it:
  `comp_redraw_max`, `mixed_box_node_share`, and `loc_box_repeat` = the literal count asked for, the mean
  number of different parents under which a `localisation` group of the same compartment is declared.
- **Direction, range:** lower is better; >= 1 (1 = every compartment drawn once).
- **Thresholds:** pass <= 1.5, warn <= 3.0, fail > 3.0 - **judgement**. Source for the principle: in SBGN
  an entity pool node belongs to exactly one compartment [11], and curated pathway figures localise
  processes to cytoplasm, mitochondria or nucleus [10]. That one drawing per compartment is the
  readable form is my reading of both.
- **Not applicable:** all nodes share one `loc_space` (NaN).
- **Gaming:** both degenerate groupings score 1.0 (one box: one patch of each compartment; singletons:
  every box is pure) - counter-rules RD01, RD02, RD03.

### RD11_BALANCE - no giant box next to dust
- **By eye:** the boxes that hold nodes are of comparable size.
- **Metric:** over boxes with `load > 0` (k of them), `p_i = load_i / N`,
  `size_balance = exp(-sum p_i ln p_i) / k` = effective number of boxes / number of boxes.
- **Direction, range:** higher is better; (0, 1].
- **Thresholds:** pass >= 0.50, warn >= 0.25, fail < 0.25 - **judgement** (community sizes in real networks
  are heavy-tailed, so 1.0 is not expected).
- **Not applicable:** fewer than two boxes hold nodes (NaN).
- **Gaming:** singletons score 1.0 and equal random chunks score high - counter-rules RD03, RD01, RD05.

### RD12_DEPTH - no more nesting than the size of the map needs
- **By eye:** a map of up to about 175 nodes needs one level of the strategy's own boxes, up to about 1,225
  two, beyond that three; compartment boxes are not counted.
- **Metric:** `depth_mean` = mean over nodes of the number of groups on `path(v)` whose basis is not
  `localisation`; `depth_needed` = smallest `d >= 1` with `7^d * 25 >= N`;
  `depth_excess = max(0, depth_mean - depth_needed)`. Reported with it: `depth_levels_mean`, `depth_max`.
- **Direction, range:** lower is better; [0, 2].
- **Thresholds:** pass 0, warn <= 1, fail > 1. Source for the direction: more depth harmed search, a
  medium breadth and depth did best [3]. The capacities 7 (boxes per level, [1]) x 25 (half of the RD02
  limit) are **judgement**.
- **Not applicable:** never.
- **Gaming:** a flat grouping scores 0 for any model size - counter-rules RD01 and RD02 (a large model
  cannot be flat and pass both: 12 boxes x 50 nodes = 600 nodes).

---

## Which rule rewards which degenerate grouping, and what stops it

| grouping | rewarded by (pass) | stopped by (fail) |
|---|---|---|
| everything in one group | RD03 (no singleton), RD06 (no clone), RD10 (1.0), RD12 (flat) | RD01, RD02, RD04, RD05, RD07, RD09 |
| every node its own group | RD02 (no large box), RD06, RD10 (1.0), RD11 (1.0), RD12 | RD01, RD03, RD04, RD05, RD07, RD09 |
| one group per compartment | RD03, RD04, RD06, RD07, RD09, RD10, RD12 | RD02 (and RD05 where a compartment has > 100 nodes) |
| random, realistic sizes | RD03, RD06, RD11, RD12 | RD04, RD05, RD07, RD09, RD10 |

Measured on the six panel models with synthetic groupings that I built (not strategy outputs).
`louvain` = Louvain communities of the species-reaction graph, communities over 50 nodes split once, each
group named after its highest-degree member; `loc+louvain` = one L1 box per compartment, compartments over
50 nodes split by Louvain.

Combined score (0..1):

| model | one | single | comp | random | louvain | loc+louvain |
|---|---|---|---|---|---|---|
| brenda_merged | 0.15 | 0.16 | 0.38 | 0.35 | 0.89 | 0.49 |
| e_coli_core | 0.17 | 0.18 | 0.39 | 0.35 | 0.87 | 0.93 |
| R-HSA-77289 | 0.17 | 0.18 | 0.39 | 0.59 | 0.93 | 0.94 |
| skm_pss | 0.17 | 0.18 | 0.38 | 0.17 | 0.41 | 0.53 |
| skm_ja | 0.17 | 0.18 | 0.57 | 0.29 | 0.84 | 0.51 |
| skm_sa | 0.17 | 0.18 | 0.41 | 0.34 | 0.73 | 0.62 |

Verdicts in rule order RD01..RD12 (`P` pass, `w` warn, `F` fail, `-` not applicable):

| model | one | single | comp | random | louvain | loc+louvain |
|---|---|---|---|---|---|---|
| brenda_merged | `FFPFFPF-F--P` | `FPFFFPF-F-PP` | `PFPPFPP-P--P` | `wFPFFPF-F-PP` | `wPPPPPwPF-PP` | `PFPPPPPFF-PP` |
| e_coli_core | `FFPFFPF-FP-P` | `FPFFFPF-FPPP` | `PFPPFPP-PPPP` | `PFPFFPF-FFPP` | `PPPPwPP-FFPP` | `PPPPwPPPFPPP` |
| R-HSA-77289 | `FFPFFPF-FP-P` | `FPFFFPF-FPPP` | `PFPPFPP-PPwP` | `PPPFFPF-FFPP` | `PPPPPPP-FPPP` | `PPPPPPPPFPPP` |
| skm_pss | `FFPFFPF-FP-P` | `FPFFFPF-FPPP` | `PFPPFPP-PPwP` | `FFPFFPF-FFPP` | `FPPPPPwwFFPP` | `PFwPPPPFFPwP` |
| skm_ja | `FFPFFPF-FP-P` | `FPFFFPF-FPPP` | `PFPPPPP-PPPP` | `PFPFFPF-FFPP` | `wPPPPPP-FFPP` | `PPFPPPPFFPwP` |
| skm_sa | `FFPFFPF-FP-P` | `FPFFFPF-FPPP` | `PFPPFPP-PPPP` | `wFPFFPF-FFPP` | `wPPPPPP-FFPP` | `PPwPPPPFFPPP` |

Neither degenerate grouping is first on any model (asserted in the tests on all six models). Every rule
is passed by at least one baseline on at least one model, so no pass zone is empty; RD09 separates a
member-named topology grouping (fail) from a compartment-named one (pass). On `skm_pss` none of my baselines reaches
0.6: 37 top-level communities fail RD01, and a cytoplasm of 1,338 nodes needs a second and third
level that these baselines do not build.

## Combining the rules into one score

1. Per rule: `score_i` in [0, 1] from the pass / warn / zero table above (`rule_scores(m)`); NaN rules are
   left out and the weights renormalised.
2. `mean = sum(w_i * score_i) / sum(w_i)` with the weights of the table (3 for the four rules a reader
   meets first, 2 for edges, clones, distinct names and fan-out, 1 to 1.5 for the rest).
3. Gate: `readability_score = mean * (0.5 + 0.5 * min(score of RD01, RD02, RD03))`. A grouping that fails
   a core structural rule keeps at most about half of its mean, so good names cannot rescue a hairball or
   a field of singletons.
4. Across models: report the mean and the minimum over the six models; a strategy is only as readable as
   its worst map. Always show the verdict string next to the number.
5. With the other two viewpoints: use `readability_score` as one sub-score, or take `rule_scores` and
   re-weight; in the second case keep the gate (or an equivalent veto on RD01-RD03), and keep RD05 and
   RD06 at equal weight - otherwise cloning everything becomes profitable.
6. `ref_scatter` and `ref_mix` (only with a curated reference; see below) are diagnostics and are not in
   the score.

Optional reference (`reference=` table of `node_id, label`; several labels per node, unlabelled nodes
ignored, species-only labelling is fine): `ref_scatter` = label-size-weighted mean over reference labels
with >= 3 placed nodes of the effective number of L1 groups that hold the label (1 = a reader finds the
curated pathway in one place); `ref_mix` = node-weighted mean over boxes of the effective number of
reference labels in the box (1 = a box shows one curated pathway). One group gives `ref_scatter` = 1 and
singletons give `ref_mix` = 1, so the two must be read together.

## What I am unsure about

1. **Most cut-offs are judgement.** Only the box limits (20 / 50 / 100 nodes), the 3..9 item range and the
   Q range come from the literature; no user study exists for nested pathway boxes. Treat verdicts as
   ordinal; the ranking of strategies should be checked for stability when the pass limits move by a third.
2. **Localisation as a free backdrop** (no lower limit on compartment counts, compartment levels not
   counted as depth) favours localisation-first groupings in RD01, RD08 and RD12, while RD10 punishes
   localisation-below-function by construction. Whether a reader finds twelve patches of "cytosol" worse
   than one cytosol box with 1,338 nodes in sub-boxes is exactly what the gallery should be used to judge.
3. **Long edge = no common box of <= 100 nodes** is a proxy for arc length without a layout. A renderer
   that places related boxes side by side makes some "long" arcs short.
4. **Informative = recall >= 0.5 within the parent.** It does not test whether the term is too general in
   absolute terms: "EC 1" on the one box that holds all oxidoreductases passes. An ontology-depth test
   would catch it but needs depth tables for EC, GO and MapMan that the contract tables do not all carry.
5. **RD09 treats `notes` names as not backed** although a notes key is an allowed name source. If the
   curators count a notes value such as a pathway name as an identifier, move `notes` into `BACKED_SOURCES`.
6. **Counting unit.** Nodes are species plus reactions; a box of 50 nodes is about 25 reactions. If the
   renderer hides modifiers or draws complexes as one glyph, the limits should move.
7. **Clone accounting** assumes a clone is drawn as a loose glyph in the named group. A strategy that
   declares no clone is not penalised by RD06 and is judged on RD05 alone; `clone_demand` is reported so
   that such strategies can still be compared.
8. `skm_pss` (1,817 nodes) may be beyond what any three-level grouping can pass on RD01, RD02 and RD08 at
   once with boxes of 50; if every strategy fails there, compare the scores rather than the verdicts.

## Sources (only what was retrieved)

1. Miller GA (1956) The magical number seven, plus or minus two: some limits on our capacity for
   processing information. Psychological Review 63(2):81-97. doi:10.1037/h0043158. (Bibliographic record
   and the "about 5 to 9" summary retrieved; the 1956 text itself was not opened.)
2. Cowan N (2001) The magical number 4 in short-term memory: a reconsideration of mental storage capacity.
   Behavioral and Brain Sciences 24(1):87-114. doi:10.1017/S0140525X01003922. (Abstract.)
3. Larson K, Czerwinski M (1998) Web page design: implications of memory, structure and scent for
   information retrieval. Proc. CHI '98. doi:10.1145/274644.274649. (Abstract.)
4. Ghoniem M, Fekete J-D, Castagliola P (2004) A comparison of the readability of graphs using node-link
   and matrix-based representations. IEEE InfoVis 2004, 17-24. doi:10.1109/INFVIS.2004.1. (Abstract.)
5. Yoghourdjian V, Yang Y, Dwyer T, Lawrence L, Wybrow M, Marriott K (2021) Scalability of network
   visualisation from a cognitive load perspective. IEEE TVCG 27(2):1677-1687.
   doi:10.1109/TVCG.2020.3030459; arXiv:2008.07944. (Abstract.)
6. Purchase H (1997) Which aesthetic has the greatest effect on human understanding? Graph Drawing '97,
   LNCS 1353:248-261. doi:10.1007/3-540-63938-1_67. (Abstract.)
7. Newman MEJ, Girvan M (2004) Finding and evaluating community structure in networks. Physical Review E
   69:026113. (The 0.3-0.7 range is quoted here from papers citing it; the original was not opened.)
   https://www.researchgate.net/publication/7659762_Finding_and_Evaluating_Community_Structure_in_Networks
8. Henry N, Bezerianos A, Fekete J-D (2008) Improving the readability of clustered social networks using
   node duplication. IEEE TVCG 14(6):1317-1324. doi:10.1109/TVCG.2008.141. (Finding as summarised in Di
   Giacomo et al., A user study on hybrid graph visualizations, doi:10.1007/978-3-030-92931-2_2.)
9. Toure V, Le Novere N, Waltemath D, Wolkenhauer O (2018) Quick tips for creating effective and impactful
   biological pathways using the Systems Biology Graphical Notation. PLoS Comput Biol 14(2):e1005740.
   doi:10.1371/journal.pcbi.1005740. (Full text read; Tip 7.)
10. Hanspers K, Kutmon M, Coort SL, et al. (2021) Ten simple rules for creating reusable pathway models for
    computational analysis and visualization. PLoS Comput Biol 17(8):e1009226.
    doi:10.1371/journal.pcbi.1009226. (Full text read; Rules 2, 3 and 7.)
11. Rougny A, Toure V, Moodie S, et al. Systems Biology Graphical Notation: Process Description language
    Level 1 Version 2.0 (clone marker; an entity pool node belongs to one compartment).
    https://www.researchgate.net/publication/333779482_Systems_Biology_Graphical_Notation_Process_Description_language_Level_1_Version_20
    (Specification text as indexed; not read in full.)
12. Vehlow C, Beck F, Weiskopf D (2017) Visualizing group structures in graphs: a survey. Computer Graphics
    Forum 36(6):201-225. https://diglib.eg.org/items/3860e224-f912-4809-a648-c9446477bd75 (Background for
    drawing hierarchical groups as nested contours; no threshold taken from it.)
