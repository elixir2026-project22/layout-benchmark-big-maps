# RULES.md - merged rulebook for judging a grouping of an SBML pathway map

BioHackathon Europe 2026, project 22 (Open Pathway Maps). One rulebook merged from the three blind
proposals in `rules/`:

| proposal | rules | implementation |
|---|---|---|
| `rules_graph.md` (viewpoint: graph) | G1-G12 | `rules/metrics_graph.py` |
| `rules_readability.md` (viewpoint: what a reader takes in) | RD01-RD12 | `rules/metrics_readability.py` |
| `rules_biology.md` (viewpoint: biology and reuse/FAIR) | BIO-1, BIO-2, LOC-1..LOC-5, FAIR-1..FAIR-4 | `rules/metrics_biology.py` |

35 proposed rules in, **27 kept** (M01-M27) and **8 dropped** as duplicates of a kept rule. Every kept
rule keeps the metric, the direction and the thresholds of the proposal it comes from; where two
proposals measured the same thing differently the conflict is stated and decided in section 2, always on
grounds of what the metric measures and never on what it does to any ranking. No threshold was changed.
Thresholds the proposals marked as judgement are marked **judgement (provisional)** here; a threshold that
proved uninformative on the panel is reported as a finding in `SCORING_NOTES.md`, not retuned.

Scored by `score_grouping.py`, which computes nothing itself except one metric (M05, section 3.1) and
otherwise calls `metrics()` of the three reference modules. The "from" line of every rule names the
module and the key it reads, so each number is traceable to the track that defined it.

```
python score_grouping.py --features panel/features/skm_ja --grouping strategies/hybrid/skm_ja \
    --reference heldout/skm_ja.reference.tsv --json skm_ja_hybrid.json
python score_grouping.py --all --panel panel/features --strategies strategies --out scores.csv
python -m pytest test_score_grouping.py -q
```

## 1. What a grouping is judged on

A grouping is one nested partition (L1-L3) of all species and reactions of one model, in the format of
`CONTRACT.md`. The rules ask five questions, in the order of the blocks below:

| block | question | rules | weight |
|---|---|---|---|
| A. Structure and granularity | is the map made of a readable number of boxes of readable size? | M01-M08 | 0.25 |
| B. Graph structure | do the drawn arcs stay inside the boxes, and are the boxes the network's own clusters? | M09-M14 | 0.20 |
| C. Localisation fidelity | can the reader see where in the cell every node is? | M15-M19 | 0.20 |
| D. Naming and reuse | does every box carry a title that is true, distinguishing and resolvable? | M20-M24 | 0.15 |
| E. Biology | does the grouping agree with curated pathways, and are members of a box alike in annotation? | M25, M26 | 0.20 |
| gate | can a third party reproduce the files? | M27 | gate only |

Definitions used throughout are those of the source proposals and are not restated here:
the bipartite species/reaction graph, the species-normalised edge weight `w_rs = 1/d_s` and the
data-derived hub set `H` (`rules_graph.md` section 1); `path(v)`, `home(v)`, `size(g)`, `load(g)` and the
**box** convention that a localisation group `c` and its sibling `c@boundary` count as one drawn box
(`rules_readability.md`, notation); `chain(i)`, `leaf(i)`, `space(i)`, `kind(i)`, the four annotation
channels and their information content (`rules_biology.md`, notation and BIO-2).

**Verdicts.** Every rule returns `pass`, `warn`, `fail` or `na`. `na` means the rule does not apply to
this model or this grouping; an `na` metric is **NaN, never 0**, is left out of the score, and the weights
of the block are renormalised over the rules that do apply.

## 2. Merge decisions

### 2.1 Duplicates removed (8 rules)

| dropped | kept instead | why the kept rule wins |
|---|---|---|
| **G3** singleton share (`singleton_group_share_leaf`, pass 0.05 / warn 0.20) | **M02** = RD03 (`singleton_share`, pass 0.10 / warn 0.25) | Same quantity, same anchor (the motivating grouping, 167 singletons of 346 groups = 0.48, fails either way). RD03 exempts a singleton that is *forced*: a localisation group that is the only group declared for its `compartment_id`. The contract makes localisation given, not chosen, so a compartment holding one node must not be charged to the strategy. G3 has no such exemption. |
| **G8** readable box size (`readable_node_share_leaf`, 3 <= n_g <= 50) | **M01** = RD02 (`oversize_node_share`) for the upper half, **M02** for the lower half | Both take the 50-node limit from Yoghourdjian et al. 2021. RD02 reports the same bound as a node share over the limit and adds the 100-node hard limit as a second diagnostic; its lower end (boxes of 1) is already M02. Keeping both would charge a large box twice. `readable_node_share_leaf` stays as a diagnostic. |
| **G9** siblings at every level (`max_siblings`, pass 2-9, fail > 50) | **M03** = RD01 (top boxes 3-12) and **M05** = RD08 (child boxes per parent) | Same concern (how many boxes stand side by side), split over the two places a reader meets it. G9 is a worst-case over all sibling sets, so one wide parent in a deep model decides the verdict; M03/M05 separate the overview from the opened box and M05 is a share, not a maximum. The conflict over the upper bound (9 vs 12) is resolved in 2.2. `max_siblings` and `fan_max` stay as diagnostics. |
| **G12** reference agreement (`ref_bcubed_f1_best`, pass 0.60) | **M25** = BIO-1 (`ref_skill`) | Same reference, same BCubed family. BIO-1 uses the extension of Amigo et al. 2009 for *overlapping* assignments, which is what the held-out reference is (41 nodes of brenda_merged carry 2-8 labels), and subtracts the best trivial grouping (one group / singletons / one box per compartment) so the value is chance-corrected. G12's own text calls its 0.60 threshold uncalibrated. `ref_bcubed_f1_best`, `ref_cut_w` and the per-level F1 values stay as diagnostics. |
| **RD05_LONGEDGE** long arcs (`long_edge_share`: no common box of <= 100 nodes) | **M09** = G1 (`cut_w_leaf`) and **M13** = G7 (`edge_span_ratio`) | Three rules measured edge locality. G1 is the boundary-crossing part with the species-normalised weight, so a map is not judged mainly on where water and ATP sit; G7 is the span part, integrated over the whole hierarchy with Dasgupta's cost. RD05 mixes both and, by its own text, is removed completely by cloning every species next to every reaction - the gaming it needs RD06 to block. G1 honours a copy only for a species in the hub set `H`, so the same trick buys nothing. Dropping RD05 therefore also removes the reason to keep RD06. `long_edge_share`, `l1_cross_share`, `cut_edge_share` and `modularity_L1` stay as diagnostics. |
| **RD06_CLONE** cloned species share (<= 0.10 of species) | **M14** = G11 (`clone_misuse_share`) | Both bound cloning. G11 asks the question a curator asks - *is this copy earning its place* (is a reaction of the species in that box) and *is this species one that is used all over the map* (membership of the data-derived set `H`) - whereas RD06 bounds the count of cloned species by a flat 10 %, a number the contract cannot supply. RD06's role as the counter-rule of RD05 lapses with RD05. `cloned_species_share`, `clone_load` and `clone_demand` stay as diagnostics. |
| **RD09_BACKED** backed names (`backed_share`, share of *groups*) | **M23** = FAIR-3 (`name_grounded`, share of *nodes*) | Same property: a title that can be looked up outside the model. FAIR-3 weighs it by nodes, which is what a reader meets (a thousand nodes under an unnamed box outweigh twenty well-named small boxes), and gives half credit to a `member` or `notes` name, which the contract explicitly allows as a name source. `backed_share` and `name_resolvable` stay as diagnostics. |
| **RD10_COMPDRAW** compartment drawn once (`comp_redraw`, container-based) | **M16** = LOC-2 (`loc_repeat`, box-based) | Same property, two assumptions about the renderer. LOC-2 counts the *declared* localisation boxes, which is what `contract.py write-sbml` exports and what the common renderer nests by `compartments.tsv`; RD10 infers a drawing from boxes that happen to be pure in one compartment, which assumes a renderer behaviour the project has not fixed. Where a strategy declares no localisation box at all, LOC-2 is `na` and M15 (LOC-1) fails instead, so nothing escapes. `comp_redraw`, `comp_redraw_max`, `mixed_box_node_share` and `loc_box_repeat` stay as diagnostics. |

### 2.2 Conflicts between surviving rules, and the decision

1. **Upper bound on sibling boxes: 9 (G9, from Miller 1956) or 12 (RD01/RD08, judgement).**
   Decided for **12**. Both numbers extend a memory-span result to boxes on a map, which neither source
   supports; the graph track itself writes that 9 "may be too strict at L2", and the readability track
   derives its range from a lower bound of 3 (Cowan 2001, chunks without rehearsal) and an upper bound
   "somewhat above nine" (Larson & Czerwinski 1998, breadth beats depth in menu hierarchies). The looser
   bound is the one with the explicit reasoning chain, and a map is scanned, not memorised. Marked
   judgement (provisional) in M03 and M05.
2. **Singleton threshold: 0.05/0.20 (G3) or 0.10/0.25 (RD03).** Decided for **RD03** - see 2.1; both
   anchor on the same example, and only RD03 excuses a compartment that holds one node.
3. **Largest box: relative (G2, share of the map) or absolute (RD02, nodes over 50).** Decided to **keep
   both** (M04 and M01): they are not the same statement. A 300-node box in a 1,800-node model passes G2
   and is unreadable; twelve boxes of 40 nodes in a 100-node model pass RD02 and one of them still holds a
   third of the map. On the panel the two diverge exactly as expected (see `SCORING_NOTES.md`).
4. **Does reference agreement enter the score?** Graph track: no (it keeps models with and without a
   reference on one scale). Biology track: yes, with weight 0.35, the largest of all. Decided: **yes, as
   block E with weight 0.20**, and the scorer additionally reports `score_core`, the same score without
   M25, so that the six models can be compared on one scale after all. Reasons for 0.20 rather than 0.35:
   the reference exists for 4 of the 6 panel models, labels only species in 3 of them, and covers 56-100 %
   of their nodes; reason for not 0: it is the only evidence that is independent of the annotation a
   strategy may have grouped by.
5. **Name truth: lift (FAIR-2) or recall-within-parent plus distinctness (RD07).** Decided to **keep
   both** (M21, M22) with the overlap declared: condition (c) of RD07 (the box holds at least half of what
   its name denotes within its parent) is a recall statement that FAIR-2's support term partly repeats.
   RD07 is kept for its conditions (a) and (b) - no two sibling boxes with the same title, no box titled
   like an ancestor - which FAIR-2 cannot see at all. Because of the overlap M22 is given weight 1 in
   block D while M21 and M23 are given 2.
6. **How a metric becomes a number: three different mappings.** Graph: pass 1 / warn 0.5 / fail 0, a step
   function. Readability: linear from the pass limit (1) through the warn limit (0.5) to a zero point.
   Biology: linear through four knots, 0 at the worst value, 0.25 at the warn limit, 0.75 at the pass
   limit, 1 at the best value. Decided for the **biology four-knot mapping for every rule** (section 4.1):
   it is continuous, it distinguishes a grouping that just passes from one that passes comfortably, and it
   does not reward reaching a threshold by a hair. The discrete verdict is reported separately, so the
   step reading of the graph track is still available in `scores.csv`.
7. **The gate: a cap (graph: score capped at 0.40 if G2 or G3 fails) or a factor (readability: score
   multiplied by 0.5 + 0.5 * min over the three core rules).** Decided to **apply both** (section 4.3):
   the factor first, because it is continuous and charges a near-fail as well, then the cap, because it
   guarantees that a map with one dominant box, a hairball or a field of singletons cannot reach the
   middle of the field whatever the rest of its numbers do. The cap covers the four structural gate rules
   M01-M04, the union of the two proposals' gate sets.
8. **Determinism as a gate (FAIR-1).** Kept as a gate (score 0), as the biology track proposed and as
   contract rule 2 requires. See section 5 for how it is set in this run: **verified input, not measured
   here.**

## 3. The rules

Each rule gives: the one-sentence statement a curator can check; the metric with its formula or a pointer
to the formula in the source proposal; direction and range; pass / warn / fail thresholds and the *zero
point* used by the utility mapping (section 4.1); when it does not apply; the source of the threshold; and
the module and key the scorer reads.

### Block A - Structure and granularity (weight 0.25)

**M01 - No box is a hairball.** [RD02_BOXSIZE] *weight 3, gate*
- Curator: no box holds more than about 50 loose nodes (species plus reactions) that are not in a sub-box.
- Metric: `oversize_node_share = sum(load(g) for g with load(g) > 50) / N`.
- Lower is better; [0, 1]. Zero point 0.70.
- pass <= 0.10, warn <= 0.30, fail > 0.30.
- Applies always (a model of <= 50 nodes passes trivially).
- Source: the 50-node box limit is from Yoghourdjian et al. 2021 (readers struggle with dense node-link
  diagrams above 50 nodes, sparse ones above 100) with Ghoniem et al. 2004 for the 20-node point; the
  share limits 0.10 / 0.30 are **judgement (provisional)**.
- From: `metrics_readability.oversize_node_share`. Diagnostics: `box_size_max`, `box_over100_node_share`,
  `box_size_median`, and `readable_node_share_leaf` (the dropped G8).

**M02 - Boxes group something.** [RD03_SINGLETON, absorbs G3] *weight 3, gate*
- Curator: at most one box in ten holds a single node.
- Metric: `singleton_share = |{g : size(g) = 1, not forced}| / |groups|`. A singleton is *forced*, and not
  counted, when it has basis `localisation` and is the only group declared for its `compartment_id`.
- Lower is better; [0, 1]. Zero point 0.50.
- pass <= 0.10, warn <= 0.25, fail > 0.25.
- Applies always.
- Source: **judgement (provisional)**; anchored on the grouping that motivated the project (167 singletons
  of 346 groups = 0.48, described as unclear), which fails.
- From: `metrics_readability.singleton_share`. Diagnostics: `singleton_node_share`,
  `metrics_graph.singleton_group_share_leaf` (the dropped G3), `residual_node_share`.

**M03 - A handful of boxes at the top level.** [RD01_TOP] *weight 3, gate*
- Curator: the map shows between 3 and 12 top-level boxes.
- Metric: `top_groups` = number of distinct L1 boxes, a localisation group and its sibling `@boundary`
  group counted as one.
- Inside the range is better; integer >= 1. Zero points 1 (below) and 40 (above).
- pass 3..12, warn 2 or 13..20, fail 1 or > 20.
- Does not apply when `N <= 20` (such a model reads as one picture, Ghoniem et al. 2004). The **lower**
  limit is waived when every L1 group has basis `localisation`: then the number of top boxes is the number
  of compartments, which the model gives (flag `top_all_localisation`); the upper limit still holds.
- Source: span of immediate memory 5-9 (Miller 1956) and about 4 chunks without rehearsal (Cowan 2001) for
  the lower end; breadth beats depth in menu hierarchies (Larson & Czerwinski 1998) for an upper end above
  nine. The numbers 12 and 20 are **judgement (provisional)**; see conflict 2.2.1.
- From: `metrics_readability.top_groups`. Diagnostics: `n_groups_L1`, `n_groups_L2`, `n_groups_L3`,
  `metrics_graph.max_siblings` (the dropped G9).

**M04 - No dominant box.** [G2] *weight 3, gate*
- Curator: no single box holds more than a third of the map.
- Metric: `largest_share_leaf = max_g n_g / N` over leaf groups (`n_g` = home members).
- Lower is better; (0, 1]. Zero point 1.0.
- pass <= 1/3, warn <= 1/2, fail > 1/2.
- Applies always; for a model too small to split (N < 9) read the verdict with care. Judged at the leaf
  only: at L1 a localisation group may legitimately hold most of the model (cytosol), and localisation is
  given, not chosen.
- Source: **judgement (provisional)**; anchored on the brief's example of a 42-node model drawn as one box
  (share 1.0).
- From: `metrics_graph.largest_share_leaf`. Diagnostics: `largest_share_L1`, `largest_box_share`.

**M05 - An opened box shows a handful of sub-boxes.** [RD08_FAN, upper bound only] *weight 2*
- Curator: a box that has sub-boxes has at most 12 of them.
- Metric: `fan_max12_share` = share of groups with at least one child group whose number of child boxes is
  at most 12 (same box convention as M03; a parent that also holds loose nodes beside child groups counts
  one extra box). The **lower** bound of RD08 (at least 2 children) is not part of this metric: it is
  M06, which charges the same defect with zero tolerance, and charging it twice would double-count.
- Higher is better; [0, 1]. Zero point 0.40.
- pass >= 0.90, warn >= 0.70, fail < 0.70.
- Does not apply to a flat grouping (no group has a child): NaN.
- Source: same as M03; the share limits are **judgement (provisional)**.
- From: computed in `score_grouping.py` from `metrics_readability._Model` (the only metric the scorer
  computes itself), because the published `fan_ok_share` mixes the lower bound in. Diagnostics:
  `fan_ok_share`, `fan_max`, `only_child_share`.

**M06 - No box that holds just one box.** [G10] *weight 1*
- Curator: no box contains exactly one box and nothing else.
- Metric: `single_child_share` = share of parent groups with exactly one child group and no node of their
  own; a nested grouping with a single top-level group counts the undrawn root as such a parent.
- Lower is better; [0, 1]. Zero point 1.0.
- pass = 0, warn <= 0.20, fail > 0.20.
- Does not apply to a grouping with one level: NaN.
- Source: **judgement (provisional)**, but the pass value of exactly 0 is forced by the contract: L2 and L3
  are optional per node, so a group that would have one child can simply have none.
- From: `metrics_graph.single_child_share`.

**M07 - No giant box next to dust.** [RD11_BALANCE] *weight 1*
- Curator: the boxes that hold nodes are of comparable size.
- Metric: over the k boxes with `load > 0` and `p_i = load_i / N`,
  `size_balance = exp(-sum p_i ln p_i) / k` (effective number of boxes over number of boxes).
- Higher is better; (0, 1]. Zero point 0.10.
- pass >= 0.50, warn >= 0.25, fail < 0.25.
- Does not apply when fewer than two boxes hold nodes: NaN.
- Source: **judgement (provisional)** (community sizes in real networks are heavy-tailed, so 1.0 is not
  expected).
- From: `metrics_readability.size_balance`.

**M08 - No more nesting than the size of the map needs.** [RD12_DEPTH] *weight 1*
- Curator: a map of up to about 175 nodes needs one level of the strategy's own boxes, up to about 1,225
  two, beyond that three; compartment boxes do not count.
- Metric: `depth_mean` = mean over nodes of the number of groups on `path(v)` whose basis is not
  `localisation`; `depth_needed` = smallest d >= 1 with `7^d * 25 >= N`;
  `depth_excess = max(0, depth_mean - depth_needed)`.
- Lower is better; [0, 2]. Zero point 2.0.
- pass = 0, warn <= 1, fail > 1.
- Applies always.
- Source: direction from Larson & Czerwinski 1998 (more depth harmed search); the capacities 7 boxes per
  level (Miller 1956) times 25 nodes per box (half the M01 limit) are **judgement (provisional)**.
- From: `metrics_readability.depth_excess`. Diagnostics: `depth_mean`, `depth_max`, `depth_needed`.

### Block B - Graph structure (weight 0.20)

**M09 - Edges stay inside boxes.** [G1] *weight 3*
- Curator: pick any species - most of the reactions it takes part in are drawn in its own box.
- Metric: `cut_w_leaf = (1/W) * sum over (r,s) of w_rs * [P_leaf(r) != P_leaf(s) and no honoured copy of s
  in P_leaf(r)]` with `w_rs = 1/d_s`, i.e. the mean over species of the share of its reactions drawn
  outside its box. A `clone_in` copy is honoured only for a species in the data-derived hub set `H`.
- Lower is better; [0, 1]. Zero point 1.0.
- pass <= 0.25, warn <= 0.50, fail > 0.50.
- Does not apply to a model without edges (NaN). A leaf level with one group is a **fail**, not `na`.
- Source: the 0.50 limit follows from Radicchi et al. 2004 - summed over groups, a set of weak communities
  has internal degree above external, i.e. cut share below 1/2. The 0.25 limit is **judgement
  (provisional)**, placed between random groupings (0.92-0.99) and module-based ones (0.02-0.33) on the
  panel. The weight `1/d_s` follows Newman 2001 and Croes et al. 2006.
- From: `metrics_graph.cut_w_leaf`. Diagnostics: `cut_raw_leaf`, `cut_w_L1..L3`, `cut_w_mean`,
  `metrics_readability.long_edge_share` and `cut_edge_share` (the dropped RD05).

**M10 - Every box is one connected piece.** [G4] *weight 2*
- Curator: inside a box you can walk from any node to any other along edges without leaving it.
- Metric: `connected_node_share_leaf` = sum over leaf groups with `n_g >= 2` of the size of the largest
  connected component of the induced subgraph, over the total size of those groups. Copies never connect
  a group.
- Higher is better; (0, 1]. Zero point 0.
- pass >= 0.95, warn >= 0.80, fail < 0.80.
- Does not apply when no group has two nodes (NaN; M02 then fails).
- Source: the requirement that a community be internally connected is from Traag, Waltman & van Eck 2019;
  the tolerances 0.95 / 0.80 are **judgement (provisional)**.
- From: `metrics_graph.connected_node_share_leaf`. Diagnostics: `disconnected_group_share_leaf`,
  `connected_nohub_node_share_leaf` (hub species removed first - not judged, because in signalling models
  the hub is often the legitimate centre of the group).

**M11 - Boxes follow the network's own clusters.** [G5] *weight 2*
- Curator: the boxes agree with the clusters a dedicated clustering algorithm finds.
- Metric: Barber's bipartite modularity with the species-normalised weight,
  `Q_k = sum_g [in(g)/W - KR_g * KS_g / W^2]`; `Q_w_best` = max over levels with >= 2 groups;
  `Q_eff = Q_w_best / Q_ref`, where `Q_ref` is the Q of a deterministic Louvain-type optimiser of the same
  function. Judged on `Q_eff`.
- Higher is better; `Q_eff` from about 0 (random) to about 1, slightly above 1 if a strategy beats the
  greedy optimiser. Zero point 0.
- pass >= 0.70, warn >= 0.40, fail < 0.40.
- Does not apply to a model without edges (NaN). One group at the leaf: **fail**.
- Source: a relative threshold is used because absolute modularity is not comparable across networks
  (Guimera, Sales-Pardo & Amaral 2004); the null model is the bipartite configuration model (Barber 2007);
  the best level is taken because modularity cannot resolve modules below its intrinsic scale (Fortunato &
  Barthelemy 2007). The limits 0.70 / 0.40 are **judgement (provisional)**, checked on the panel.
- From: `metrics_graph.Q_eff`. Diagnostics: `Q_w_best`, `Q_raw_best`, `Q_ref`, `n_groups_ref`,
  `metrics_readability.modularity_L1`.

**M12 - Each box is a community.** [G6] *weight 1*
- Curator: each box has more edge weight inside it than leaving it.
- Metric: for a leaf group, conductance `phi(g) = cut(g) / min(vol(g), 2W - vol(g))`; g is a weak community
  when `2*in(g) > cut(g)`. `weak_node_share_leaf` = share of nodes in weak communities.
- Higher is better; [0, 1]. Zero point 0.
- pass >= 0.80, warn >= 0.50, fail < 0.50.
- Does not apply to a model without edges (NaN). One group at the leaf: **fail**.
- Source: the per-group criterion is the weak-community definition of Radicchi et al. 2004; conductance is
  among the structural scores that best recover ground-truth communities (Yang & Leskovec 2012). The
  shares 0.80 / 0.50 are **judgement (provisional)**.
- From: `metrics_graph.weak_node_share_leaf`. Diagnostic: `conductance_mean_leaf`.

**M13 - Crossing edges stay close.** [G7] *weight 2*
- Curator: an edge that leaves its box ends in a sibling box under the same parent, not on the far side of
  the map.
- Metric: Dasgupta's cost of the group tree, normalised. `box(e)` = the smallest group holding both ends
  of edge e (an honoured hub copy counts as an end), the whole model if none;
  `edge_span_w = sum_e w_e * |box(e)| / (N * W)`; `edge_span_null` = the same with every (reaction,
  species) pair as an edge; `edge_span_ratio` = their quotient. This is how the levels are combined into
  one number: a crossing costs the size of the lowest box that still holds both ends, so a dummy level
  gains nothing.
- Lower is better; (0, about 1]. Zero point 1.0 (1 = no better than a random filling).
- pass <= 1/3, warn <= 2/3, fail > 2/3.
- Does not apply to a model without edges (NaN).
- Source: the cost function is Dasgupta 2016; the limits are **judgement (provisional)**, checked on the
  panel. Note the proposal's own caveat: the ratio is easier to lower in large models.
- From: `metrics_graph.edge_span_ratio`. Diagnostics: `edge_span_w`, `edge_span_null`.

**M14 - Copies are spent on hubs and every copy is used.** [G11] *weight 1*
- Curator: only species used all over the map are drawn more than once, and every copy has one of its
  reactions in the same box.
- Metric: `clone_misuse_share` = share of `clone_in` entries (species, group) that are idle (no reaction of
  the species lies in that group) or belong to a species outside the hub set `H`. `H` is derived from the
  degree distribution: `d_s > c`, c the smallest integer with `P(Poisson(lambda) > c) <= 0.05/|S|`, lambda
  the mean degree of the non-hub species, iterated to a fixed point.
- Lower is better; [0, 1]. Zero point 1.0.
- pass <= 0.10, warn <= 0.30, fail > 0.30.
- Does not apply when `clone_in` is not used (NaN). Not using copies is not penalised.
- Source: **judgement (provisional)**. The reason to keep low-degree species single is Guimera & Amaral
  2005 (the low-degree metabolites that link modules are better conserved than hubs); deriving the hub set
  from the network rather than from a list of names follows Huss & Holme 2007 and is required by contract
  rule 1.
- From: `metrics_graph.clone_misuse_share`. Diagnostics: `clone_idle_share`, `clone_nonhub_share`,
  `clones_per_species`, `clone_absorbed_w`, `hub_degree_cutoff`, `n_hubs`,
  `metrics_readability.cloned_species_share` and `clone_demand` (the dropped RD06).

### Block C - Localisation fidelity (weight 0.20)

The contract fixes where every node is (`loc_space`, `loc_kind`, `loc_outer` of `nodes.tsv`) and fixes the
containment tree (`compartments.tsv`, with `parent_id`, `encloses_id` and the `source` of each link).
These five rules therefore test *compliance with given information*, which is why four of them have a pass
value of exactly 1.00: any strategy can reach it by construction.

**M15 - Localisation is expressed.** [LOC-1] *weight 0.25 of the block*
- Curator: in a model with several compartments, every node sits inside a box declared as its own
  compartment.
- Metric: `loc_expressed = |{i : some g in chain(i) is a localisation box whose space equals space(i)}| / N`.
- Higher is better; [0, 1]. Zero point 0.
- pass = 1.00, warn >= 0.90, fail < 0.90.
- Does not apply to models with fewer than two populated spaces (NaN) - of the panel, `brenda_merged`.
- Source: the pass value follows from contract rule 5 ("localisation is given, not chosen") and from SBGN
  PD, where an entity pool node belongs to a single compartment; the 0.90 warn limit is **judgement
  (provisional)**.
- From: `metrics_biology.loc_expressed`. Diagnostic: `loc_pure_leaf`.

**M16 - Each compartment is drawn as one box.** [LOC-2] *weight 0.25 of the block*
- Curator: each compartment, and each compartment boundary, appears as one box.
- Metric: for every node with a box, `share(i) = |{j : b(j) = b(i)}| / |{j : compartment_id(b(j)) =
  compartment_id(b(i))}|`, `loc_repeat = 1 / mean_i share(i)` - the effective number of boxes per
  compartment (inverse Simpson, node-weighted).
- Lower is better; from 1 upwards. Mapped as `1/loc_repeat` for the utility (section 4.1).
- pass = 1.00, warn <= 2.00, fail > 2.00.
- Does not apply to single-compartment models, or to a grouping with no localisation box at all (NaN; M15
  fails instead).
- Source: **judgement (provisional)**; one box per compartment is what the contract promises the reader
  ("organelles inside cytosol, extracellular around it") and what a compartment container is in a Process
  Description map. The limit 2 means "at worst two equal halves".
- From: `metrics_biology.loc_repeat`. Diagnostics: `loc_boxes_mean`, `loc_boxes_max`,
  `metrics_readability.comp_redraw` and `loc_box_repeat` (the dropped RD10).

**M17 - Transport sits on the boundary it crosses.** [LOC-3] *weight 0.20 of the block*
- Curator: a reaction that moves something across a compartment boundary, and a species that sits in a
  membrane, is drawn in the boundary box of that compartment, not among its interior reactions.
- Metric: with `B = {i : kind(i) = boundary}`,
  `loc_boundary = |{i in B : some g in chain(i) has compartment_id = space(i) + "@boundary"}| / |B|`.
- Higher is better; [0, 1]. Zero point 0.
- pass = 1.00, warn >= 0.80, fail < 0.80.
- Does not apply when no node has `loc_kind = boundary` (NaN), or to single-compartment models.
- Source: transport drawn as a process between the pools of two compartments is the Process Description
  convention (SBGN bricks dictionary); pass = 1.00 because `loc_kind` is given in `nodes.tsv`; the 0.80
  warn limit is **judgement (provisional)**.
- From: `metrics_biology.loc_boundary`. Diagnostic: `loc_boundary_in_interior`.

**M18 - Enclosing compartments are drawn around each node.** [LOC-4] *weight 0.15 of the block*
- Curator: around every node, each compartment that encloses its own compartment is drawn exactly once in
  the same part of the map: organelles inside cytosol, cytosol inside the extracellular space.
- Metric: `A(i)` = the populated proper ancestors of `space(i)` in `compartments.tsv`; `scope(i)` = the
  parent group of the outermost box of i (the top of the map if that box is at L1); an ancestor a is
  *drawn* for i when among the groups whose parent is `scope(i)` exactly one interior box has compartment
  a (or no interior box and exactly one boundary box).
  `loc_nesting = sum_i |{a in A(i) : a drawn for i}| / sum_i |A(i)|`.
- Higher is better; [0, 1]. Zero point 0.
- pass = 1.00, warn >= 0.90, fail < 0.90.
- Does not apply to single-compartment models or where no node has a populated enclosing space (NaN).
- Source: **judgement (provisional)**. The containment tree itself is not in question: SBML Level 3 dropped
  the `outside` attribute, so containment comes from `compartments.tsv`, the same for every strategy, and
  each link there carries its own `source` column.
- **Assumption carried from the proposal:** nesting is assumed to happen among boxes that share a parent.
  If the common renderer inserts a missing enclosing compartment itself, this rule is too strict for
  function-first hierarchies. Recorded as a limit, not corrected.
- From: `metrics_biology.loc_nesting`.

**M19 - A box spanning compartments holds their transport and is split by place.** [LOC-5] *weight 0.15 of the block*
- Curator: a group that holds nodes from several compartments also holds the transport reactions that
  connect those compartments, and is divided into compartment boxes one level down.
- Metric: a group *spans* when its members have two or more `loc_space` values; it is *linked* when
  joining, for every member reaction, the spaces of its reactants, products and its own `loc_space` (a
  membrane counting as the space it encloses) connects all spaces of the group.
  `bad(i)` = i lies in a spanning group that is not linked, or in a spanning group and has no box;
  `loc_span_ok = 1 - |{i : bad(i)}| / N`.
- Higher is better; [0, 1]. Zero point 0. Never spanning is allowed and scores 1.00: the rule forbids
  unjustified spans, it does not ask for spans.
- pass = 1.00, warn >= 0.90, fail < 0.90.
- Does not apply to single-compartment models (NaN).
- Source: **judgement (provisional)**. This is the rulebook's answer to "when is a functional group across
  compartments acceptable": its own transport steps connect its compartments (M19), every member still
  sits in a compartment box (M15, M19), and the resulting repetition stays within M16.
- From: `metrics_biology.loc_span_ok`. Diagnostics: `loc_span_share`, `loc_span_connected`.

### Block D - Naming and reuse (weight 0.15)

**M20 - Every box has a title.** [RD04_NAMED] *weight 2*
- Curator: nine boxes in ten carry a label.
- Metric: `named_share = |{g : label(g) != ""}| / |groups|`.
- Higher is better; [0, 1]. Zero point 0.20.
- pass >= 0.90, warn >= 0.60, fail < 0.60.
- Applies always (the contract allows a `member` name for any group).
- Source: the principle (graphical annotations organise a layout and should be labelled) is from Hanspers
  et al. 2021; the cut-offs are **judgement (provisional)**, anchored on 34 named groups of 346 (0.10) in
  the motivating grouping.
- From: `metrics_readability.named_share`. Diagnostics: `named_node_share`, `member_named_share`,
  `metrics_biology.name_coverage`.

**M21 - A title is true of the box's members and specific to them.** [FAIR-2] *weight 2*
- Curator: the identifier behind a group's name is carried by the members of that group, and by clearly
  fewer nodes outside it.
- Metric: for a named group with at least two members, `support_g = |E_g and Car| / |E_g|`,
  `prev = |Car| / |nodes of the kinds in Car|`, `lift_g = (support_g - prev) / (1 - prev)`;
  `name_lift = sum_g |E_g| * lift_g / sum_g |E_g|` over the groups with `prev < 1`. A name carried by
  every node of its kinds is vacuous and left out (counted in `name_vacuous_share`).
- Higher is better; at most 1, negative when the name fits the rest of the model better than the group.
  Zero point 0.
- pass >= 0.61, warn >= 0.21, fail < 0.21.
- Does not apply when no informative named group of two or more members exists (NaN). The statistic has
  the form (observed - expected) / (1 - expected).
- Source: the bands 0.21 / 0.61 borrow the "fair" and "substantial" limits of Landis & Koch 1977 by
  analogy of form; that scale is itself a convention, so **judgement (provisional)**.
- From: `metrics_biology.name_lift`. Diagnostics: `name_support`, `name_traceable`, `name_vacuous_share`.

**M22 - A title tells this box from its neighbours.** [RD07_DISTINCT] *weight 1 (overlap with M21 declared)*
- Curator: no two sibling boxes share a title, no box repeats the title of the box around it, and a box
  called X holds at least half of the X its parent box contains.
- Metric: over named groups, `distinct_share` = share that are (a) the only one with their label among
  their siblings (case- and whitespace-insensitive; a compartment box and its own `@boundary` box are not
  duplicates), (b) not labelled like any ancestor, and (c) have `recall(g) = |carriers(x) in g| /
  |carriers(x) in parent(g)| >= 0.5`; (c) is skipped for `member` and `notes` names.
- Higher is better; [0, 1]; 0 when no group is named. Zero point 0.50.
- pass >= 0.95, warn >= 0.80, fail < 0.80.
- Applies always.
- Source: **judgement (provisional)**. Weight reduced to 1 because condition (c) overlaps M21 (conflict
  2.2.5); conditions (a) and (b) are what M21 cannot see.
- From: `metrics_readability.distinct_share`. Diagnostics: `sibling_dup_share`, `same_as_ancestor_share`,
  `name_low_recall_share`, `name_recall_mean`, `global_dup_share`.

**M23 - Most nodes sit under a name backed by a resolvable identifier.** [FAIR-3] *weight 2*
- Curator: pick any node - at least one of the boxes around it has a name whose identifier can be looked up
  outside the model (SBO, EC, GO, a database cross-reference, or a compartment mapped to GO).
- Metric: a group is *well named* when it is named, has at least two members, its identifier is carried by
  a member, is not vacuous and resolves; *locally named* when all of that holds except resolution (a
  `member` or `notes` name). `name_grounded` = share of nodes with a well named group in their chain plus
  0.5 times the share of nodes with only locally named groups in their chain.
- Higher is better; [0, 1]. Zero point 0.
- pass >= 0.80, warn >= 0.50, fail < 0.50.
- Applies always; 0 when nothing is named.
- Source: FAIR principles (Wilkinson et al. 2016) for persistent identifiers and FAIR vocabularies, and
  SBML Level 3 Groups (Hucka & Smith 2016) for annotating a group with what it represents. The limits and
  the half credit are **judgement (provisional)**.
- From: `metrics_biology.name_grounded`. Diagnostics: `name_resolvable`, `named_group_share`,
  `metrics_readability.backed_share` (the dropped RD09).

**M24 - A node without annotation is still placed with a reaction partner.** [FAIR-4] *weight 1.5*
- Curator: a species or reaction with no ontology annotation is in the same box as something it reacts
  with, not alone and not in a rest group of unrelated leftovers.
- Metric: with `U` = nodes with no term in any annotation channel and at least one neighbour in
  `edges.tsv`, `stranded_poor = |{i in U : no neighbour of i is in the same innermost group (or in a
  `clone_in` group of i)}| / |U|`.
- Lower is better; [0, 1]. Zero point 1.0.
- pass <= 0.10, warn <= 0.30, fail > 0.30.
- Does not apply when no unannotated node has a neighbour (NaN).
- Source: **judgement (provisional)**, from the proposal's calibration (random groupings 0.58-0.94,
  all-singletons 1.00, non-trivial controls 0.00-0.33).
- From: `metrics_biology.stranded_poor`. Diagnostics: `stranded_annotated`, `poor_share`.
- Not implemented, carried from the proposal as a limit: the stronger test is to blank the annotation of a
  fixed 30 % of nodes and rerun the strategy. That needs the strategies to be rerun, which this scorer
  does not do.

### Block E - Biology (weight 0.20)

**M25 - Agreement with the curated reference.** [BIO-1, absorbs G12] *weight 0.6 of the block*
- Curator: nodes the curators put in one pathway share a box, and a box does not mix pathways.
- Metric: extended BCubed precision, recall and F1 over the labelled nodes (Amigo et al. 2009), with
  `C_d(i)` = the group of i at depth d plus its `clone_in` groups; `ref_f1` = the best of d in {1,2,3};
  `ref_skill = (ref_f1 - F0) / (1 - F0)` with `F0` the best F1 of the trivial groupings on the same labels
  (one group, singletons, one box per compartment with and without boundary boxes).
- Higher is better; <= 1; 0 means no better than the best trivial grouping, negative means worse. Zero
  point 0 (negative values are clipped to 0 for the utility, and reported unclipped).
- pass >= 0.25, warn >= 0.05, fail < 0.05.
- Does not apply when the model has no reference or fewer than two labelled nodes (NaN). On the SKM models
  only species are labelled, so the rule says nothing about where reactions are placed there.
- Source: **judgement (provisional)**; no published threshold for BCubed F1 was found. Calibrated by the
  proposal on a *planted* reference, not on the curated one, which it never saw.
- From: `metrics_biology.ref_skill`. Diagnostics: `ref_f1`, `ref_precision`, `ref_recall`,
  `ref_best_level`, `ref_f1_L1..L3`, `ref_baseline_f1`, `ref_f1_one_group`, `ref_f1_singletons`,
  `ref_f1_compartments`, `ref_labelled_share`, `metrics_graph.ref_bcubed_f1_best` and `ref_cut_w` (the
  dropped G12), `metrics_readability.ref_scatter` and `ref_mix`.

**M26 - Members of a box are alike in annotation.** [BIO-2] *weight 0.4 of the block*
- Curator: two annotated members of the same innermost box are more alike - same reaction type, enzyme
  class, GO term, pathway or bin code - than two annotated nodes picked from anywhere in the model.
- Metric: Resnik similarity over four channels (`sbo` for reactions, `ec`, `go`, `class`), each closed
  under ancestors where a hierarchy exists, with information content `IC(t) = ln(N_ch/n_t)/ln(N_ch)`;
  `hom = (mean_i S_in(i) - mean_i S_all(i)) / (1 - mean_i S_all(i))`, where `S_in` averages over the
  comparable members of the node's innermost group and `S_all` over the whole model.
- Higher is better; at most 1; 0 is chance; negative when boxes mix unlike nodes more than chance. Zero
  point 0.
- pass >= 0.21, warn >= 0.05, fail < 0.05.
- Does not apply when fewer than two comparable annotated nodes exist (NaN). A channel whose only term is
  carried by every node contributes nothing (all 218 `brenda_merged` reactions are "biochemical reaction").
- Source: the form (observed - expected)/(1 - expected) and the 0.21 band are the Landis & Koch 1977 "fair"
  limit by analogy - **judgement (provisional)**; the fail limit 0.05 is calibration (random groupings
  -0.032 to 0.047 over the six models, 20 seeds each). Information content and Resnik similarity follow
  Pesquita et al. 2009; the chance correction follows the reasoning of Rosenberg & Hirschberg 2007 and
  Vinh et al. 2010.
- Carried caveat: on the SKM models reactions carry only an SBO type and the best value the proposal could
  produce was 0.18, below the pass limit, so **0.21 may be out of reach there for honest reasons**. Read
  with `hom_gain` (the part not already given by the compartments).
- From: `metrics_biology.hom`. Diagnostics: `hom_L1`, `hom_sbo`, `hom_ec`, `hom_go`, `hom_class`,
  `hom_compartments`, `hom_gain`, `annotated_share`.

### Gate

**M27 - Two runs give identical files.** [FAIR-1] *gate*
- Curator: running the strategy twice on the same model gives byte-identical `groups.tsv` and
  `groups_meta.tsv`.
- Metric: `fair_deterministic` = 1 if both files have the same SHA-256 in the two runs, else 0.
- pass = 1, fail = 0. A fail sets the score of that model to 0.
- Does not apply when no second run is supplied (NaN, not counted).
- Source: contract rule 2; recording seeds is rule 6 of Sandve et al. 2013.
- From: `metrics_biology.fair_deterministic` when `--rerun` is given; otherwise from the provenance flag
  described in section 5. Diagnostics: `fair_rerun_same_rows`, `fair_sorted`, `fair_rule_share`.

## 4. How the rules combine into one score per model

Implemented in `score_grouping.py` as `rule_table()`, `block_scores()` and `score()`.

### 4.1 From a metric to a utility (one mapping for every rule)

Each rule declares four knots: its *worst* value, its *warn* threshold, its *pass* threshold and its *best*
value. The utility is piecewise linear through

```
worst -> 0.00      warn -> 0.25      pass -> 0.75      best -> 1.00
```

clipped to [worst, best]. Where the pass threshold equals the best value (M06, M08, M15, M16, M17, M18,
M19) the knots coincide and passing earns the full 1.00. `loc_repeat` (M16) is inverted to `1/loc_repeat`
first, so that its unbounded worst value becomes 0. For the one range rule (M03) a value inside 3..12 earns
1.00 and a value outside it is mapped one-sidedly to 0.25 at the warn edge and 0 at the zero point; the
lower bound is waived when every L1 group is a localisation group. Negative `ref_skill` (M25) and negative
`hom` (M26) are clipped to 0 in the utility and reported unclipped in `scores.csv`.

A rule's **verdict** is read from the metric and its own thresholds, not from the utility: `pass`, `warn`
(between the pass and warn limits), `fail` (beyond the warn limit), `na` (NaN metric).

### 4.2 Weights

Inside a block the applicable rules are averaged with the weights below (a rule that does not apply is
dropped and the remaining weights are rescaled). The blocks are then averaged with the block weights; a
block in which no rule applies is dropped and the block weights are rescaled.

| block | weight | rules (weight inside the block) |
|---|---|---|
| A. Structure and granularity | 0.25 | M01 3, M02 3, M03 3, M04 3, M05 2, M06 1, M07 1, M08 1 |
| B. Graph structure | 0.20 | M09 3, M10 2, M11 2, M12 1, M13 2, M14 1 |
| C. Localisation fidelity | 0.20 | M15 0.25, M16 0.25, M17 0.20, M18 0.15, M19 0.15 |
| D. Naming and reuse | 0.15 | M20 2, M21 2, M22 1, M23 2, M24 1.5 |
| E. Biology | 0.20 | M25 0.6, M26 0.4 |

The within-block weights are those of the source proposals, with two changes, both argued in section 2 and
both fixed before any strategy was scored: M22 drops to 1 because of its declared overlap with M21
(conflict 2.2.5), and M04 (the graph track's gate, weight 2 there) rises to 3 to match the weight of the
other three gate rules of block A.

Block weights are judgement. The reasoning, in the order of the blocks:

- **A 0.25** holds the three rules the readability track weighted highest (RD01-RD03, weight 3 each) and
  the graph track's own gate; it is also what the project's motivating failure was about (167 singletons of
  346 groups, 34 of 346 groups named).
- **B 0.20**: the only block that measures the arcs that will actually be drawn.
- **C 0.20**: compliance with information the contract declares given; a map that misplaces compartments is
  wrong, not merely unclear.
- **D 0.15**: the naming rules are close to pass/fail for a careful strategy and separate little - the
  biology track's own reason for weighting its FAIR block 0.15.
- **E 0.20**: the reference is the only evidence independent of the annotation a strategy grouped by, but
  it covers 4 of 6 panel models and labels only species in 3 of them (conflict 2.2.4). The biology track's
  0.35 is not carried over.

### 4.3 Gates

Applied in this order:

1. **Determinism (M27).** `fair_deterministic = 0` sets the model's score to 0. `na` leaves the score
   untouched.
2. **Core-structure factor.** `gated = mean * (0.5 + 0.5 * min(utility of M01, M02, M03, M04))` over the
   gate rules that apply. A grouping that fails a structural gate keeps at most about half of its mean, so
   good names cannot rescue a hairball or a field of singletons.
3. **Cap.** If any of M01-M04 has verdict `fail`, the score is capped at 0.40. A map with one dominant box,
   a hairball or mostly single-node boxes is not a grouping, whatever its edges do.

`score_core` - the same computation with M25 removed from block E - is reported next to the score so that
models with and without a reference stay comparable (conflict 2.2.4).

### 4.4 Across models

- Report the per-model score, its `score_core`, the block scores and the full 27-verdict string for every
  (strategy, model) pair: one row of `scores.csv`.
- Rank strategies by the **mean score over the models they were run on**; break ties by the minimum over
  models, then by `edge_span_ratio` (M13), as the graph track proposed.
- Also report the minimum: a strategy is only as good as its worst map.
- Do not average raw metric values across rules or models; their scales differ.
- A strategy that is missing for a model is not scored for it and is not imputed. The baseline exists for
  4 of the 6 panel models, so its mean is over 4 models and is marked as such in `scores_summary.csv`
  (`n_models`).
- **Contaminated reference.** For `brenda_merged` and `skm_pss` the baseline was built from the same
  curated labels that are held out here. Its M25 is therefore not a prediction. The scorer marks those
  rows `ref_contaminated = 1`, keeps the value, and reports the baseline's mean twice: with M25 and with
  M25 forced to `na` (`mean_score_refblind`). Strategy rows are never marked.

## 5. Determinism in this run: verified input, not measured

The bundle contains one run per strategy, so this scorer cannot observe a second run. Byte-identical
reruns were verified by the project lead for all 8 strategies x 6 models before the bundle was assembled.
M27 is therefore recorded as **passed on verified input**: every strategy row carries
`fair_deterministic = 1` with `determinism_source = verified_input`. The baseline is a converted artefact
of an earlier pipeline rather than a strategy run, and no rerun of it was verified, so its rows carry
`fair_deterministic = NaN`, `determinism_source = not_verified`, and M27 is `na` for them - it does not
gate the baseline either way.

Given `--rerun <dir2>`, the scorer compares the two directories byte for byte and sets the metric from the
comparison, with `determinism_source = measured`. Nothing in `scores.csv` was produced that way.

## 6. Dropped rules

| dropped rule | measured | kept rule that covers it | reason (short) |
|---|---|---|---|
| G3 | singleton share | M02 (RD03) | duplicate; RD03 exempts compartments that hold one node, as the contract requires |
| G8 | readable box size 3..50 | M01 (upper), M02 (lower) | duplicate of the 50-node bound and of the singleton rule; keeping it charges a large box twice |
| G9 | siblings at every level | M03 (top), M05 (fan) | duplicate; a worst-case over all sibling sets lets one wide parent decide the verdict |
| G12 | reference agreement (plain BCubed F1) | M25 (BIO-1) | duplicate; BIO-1 handles multi-label references and is chance-corrected |
| RD05_LONGEDGE | long arcs | M09 (cut), M13 (span) | duplicate of two kept rules, and removable by cloning everything |
| RD06_CLONE | share of cloned species | M14 (G11) | duplicate; G11 tests whether a copy earns its place, with a data-derived hub set |
| RD09_BACKED | backed names, per group | M23 (FAIR-3) | duplicate; FAIR-3 weighs the same property by nodes and credits contract-allowed local names |
| RD10_COMPDRAW | compartment drawn once | M16 (LOC-2) | duplicate; LOC-2 counts declared boxes rather than inferring the drawing from box purity |

All eight stay in `scores.csv` as diagnostics, so no information is lost - only the double counting.

## 7. Sources

Rules that carry a citation take it from the proposal that defined them; the full records are in
`rules/rules_graph.md` section 7, `rules/rules_readability.md` ("Sources") and `rules/rules_biology.md`
("Sources"). Short form, by rule:

- M01, M03 (upper end): Yoghourdjian et al. 2021; Ghoniem, Fekete & Castagliola 2004.
- M03, M05, M08 (lower end and depth): Miller 1956; Cowan 2001; Larson & Czerwinski 1998.
- M09, M12: Radicchi et al. 2004; Yang & Leskovec 2012; weight from Newman 2001 and Croes et al. 2006.
- M10: Traag, Waltman & van Eck 2019.
- M11: Barber 2007; Guimera, Sales-Pardo & Amaral 2004; Fortunato & Barthelemy 2007.
- M13: Dasgupta 2016.
- M14: Guimera & Amaral 2005; Huss & Holme 2007.
- M15-M19: SBGN PD Level 1 (entity pool node in one compartment); SBGN bricks dictionary (transport as a
  process between two pools); libSBML Compartment documentation (SBML L3 has no `outside` attribute).
- M20, M23: Hanspers et al. 2021; Wilkinson et al. 2016; Hucka & Smith 2016.
- M21, M26 (band limits, by analogy of form only): Landis & Koch 1977.
- M25: Amigo et al. 2009.
- M26: Pesquita et al. 2009; Rosenberg & Hirschberg 2007; Vinh et al. 2010.
- M27: Sandve et al. 2013.

Dating note: Yoghourdjian et al. is dated **2021** here, after the journal issue both source proposals
cite in their reference lists (IEEE TVCG 27(2):1677-1687, doi:10.1109/TVCG.2020.3030459). The inline
comment of `metrics_readability.py` (`BOX_MAX = 50`) dates the same work 2020, after its preprint
(arXiv:2008.07944). Same study, same 50-node limit; only the year differs.

Everything else in this rulebook is marked **judgement (provisional)** at the rule where it is used. Three
tracks wrote those marks honestly; none was upgraded to a citation here. M02, M04, M06, M07, M16, M17, M18,
M19, M22, M24 and all share thresholds rest on judgement alone.

## 8. Known limits carried from the proposals

1. No threshold in this rulebook comes from a user study of nested pathway maps. Verdicts are ordinal;
   check the ranking for stability when the pass limits move by a third.
2. M11, M13 and M09 favour topology-based strategies by construction; function-first groupings will warn or
   fail on them for honest reasons. That is why block B carries 0.20 and not more.
3. M16 and M18 assume the renderer nests compartment boxes among siblings. If the common renderer inserts
   the missing enclosing compartment, both are too strict for function-first hierarchies.
4. M26 rewards grouping by the very annotation it reads; M25 is the only counter-evidence, and it exists
   for 4 of 6 models.
5. The species-normalised weight of M09 makes the cut easier in models with many leaf species (the SKM
   models) than in metabolic models, under one threshold.
6. M25 takes the best of three levels, which is fair to a strategy that works blind and lenient to a
   hierarchy whose other levels are poor; the per-level values are reported.
7. M24's real test (blank 30 % of the annotation and rerun) is not implemented here.
8. `skm_pss` (1,817 nodes) may be beyond what any three-level grouping can pass on M01, M03 and M05 at
   once; where every strategy fails a rule there, compare the scores rather than the verdicts.
