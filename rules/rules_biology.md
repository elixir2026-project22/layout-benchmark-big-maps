# Rules for judging a grouping — viewpoint: biology and reuse (FAIR)

BioHackathon Europe 2026, project 22 (Open Pathway Maps). Written blind: no strategy output was seen or used.
Reference implementation: `metrics_biology.py`; tests: `test_metrics_biology.py`.

```
python metrics_biology.py panel/features/MODEL strategies/NAME/out/MODEL --reference REF.tsv --rerun SECOND_RUN_DIR
pytest -q test_metrics_biology.py
```

## Reading guide

A grouping is judged on four questions, in this order of weight.

1. Does it agree with what curators call a pathway (BIO-1, only where a held-out reference exists)?
2. Do members of a box share biology that can be read from the SBML annotation (BIO-2)?
3. Can the reader see where every node is in the cell (LOC-1 to LOC-5)?
4. Can a third party reproduce it and trace every name to an identifier (FAIR-1 to FAIR-4)?

Eleven rules follow, most important first. The contract validator (`common/contract.py`) already rejects a
localisation group that holds a node of another compartment, a node without a rule id, and a named group without
`name_xref`. No rule here repeats those checks; every rule covers something the validator lets through.

| # | id | one-line statement | metric | better | pass | warn | fail |
|---|---|---|---|---|---|---|---|
| 1 | BIO-1 | groups agree with the curated pathways | `ref_skill` | higher | >= 0.25 | 0.05 – 0.25 | < 0.05 |
| 2 | LOC-1 | every node sits in a box declared as its compartment | `loc_expressed` | higher | = 1.00 | 0.90 – 1.00 | < 0.90 |
| 3 | BIO-2 | members of a box are alike in annotation | `hom` | higher | >= 0.21 | 0.05 – 0.21 | < 0.05 |
| 4 | LOC-2 | each compartment is drawn as one box | `loc_repeat` | lower | = 1.00 | 1.00 – 2.00 | > 2.00 |
| 5 | LOC-3 | transport sits on the boundary it crosses | `loc_boundary` | higher | = 1.00 | 0.80 – 1.00 | < 0.80 |
| 6 | LOC-4 | enclosing compartments are drawn around each node | `loc_nesting` | higher | = 1.00 | 0.90 – 1.00 | < 0.90 |
| 7 | FAIR-1 | two runs give identical files | `fair_deterministic` | higher | = 1 | – | = 0 (gate) |
| 8 | FAIR-2 | a group's name is true of its members and tells them apart | `name_lift` | higher | >= 0.61 | 0.21 – 0.61 | < 0.21 |
| 9 | LOC-5 | a box spanning compartments holds their transport and is split by place | `loc_span_ok` | higher | = 1.00 | 0.90 – 1.00 | < 0.90 |
| 10 | FAIR-3 | most nodes sit under a name backed by a resolvable identifier | `name_grounded` | higher | >= 0.80 | 0.50 – 0.80 | < 0.50 |
| 11 | FAIR-4 | a node without annotation is still placed with a reaction partner | `stranded_poor` | lower | <= 0.10 | 0.10 – 0.30 | > 0.30 |

### Notation

- `N` nodes (species and reactions) of `nodes.tsv`. `chain(i)` = the groups of node `i` from L1 down; `leaf(i)` = the
  deepest of them (the innermost box drawn around the node); `members(g)` = all nodes with `g` in their chain.
- `space(i)` = `loc_space`, `kind(i)` = `loc_kind` of `nodes.tsv`. A **box** is a group with `basis = localisation`;
  its `compartment_id` is either a space `c` (interior box) or `c@boundary` (boundary box).
- A space is **populated** if at least one node has it as `loc_space`. Virtual spaces are never populated.
- A model is **multi-compartment** if it has two or more populated spaces. Of the panel, `brenda_merged` is not.
- "judgement" marks a threshold that is my own choice. "calibration" refers to the synthetic groupings of the
  section *Behaviour on synthetic groupings*, built on all six panel models.

---

## BIO-1 — Agreement with the curated reference

**Curator's check.** Nodes that the curators assigned to the same pathway sit in the same group, and nodes of
different pathways do not, at one of the levels of the grouping.

**Metric.** Extended BCubed precision, recall and F1 (Amigó et al. 2009), over the labelled nodes only.
`L(i)` = set of reference labels of node `i`. `C_d(i)` = set holding the group of `i` at depth `d` (its deepest group
if its chain is shorter), plus every `clone_in` group of `i`, lifted to depth `d` when it lies deeper.

```
precision_d = mean_i  mean_{j : C_d(i) ∩ C_d(j) ≠ ∅}  min(|C_d(i) ∩ C_d(j)|, |L(i) ∩ L(j)|) / |C_d(i) ∩ C_d(j)|
recall_d    = mean_i  mean_{j : L(i) ∩ L(j) ≠ ∅}      min(|C_d(i) ∩ C_d(j)|, |L(i) ∩ L(j)|) / |L(i) ∩ L(j)|
F_d         = 2 · precision_d · recall_d / (precision_d + recall_d)        (j = i is included, as in BCubed)
ref_f1      = max over d in {1, 2, 3} of F_d          (ref_best_level = the d that wins)
ref_skill   = (ref_f1 − F0) / (1 − F0),   F0 = max F of the trivial groupings on the same labels:
              one group, every node alone, one group per compartment (with and without boundary boxes)
```

**Why this measure.** It is defined per item, so a few large groups do not dominate. It has a published extension
to overlapping assignments, which is what a multi-label reference and `clone_in` are. Entropy-based and
pair-counting measures do not apply directly to overlapping clusters (Amigó et al. 2009, section 6). Unlabelled
nodes are left out of both averages, so a partially labelled reference neither rewards nor punishes them.
BCubed is not corrected for chance; `ref_skill` supplies the correction by subtracting the best trivial grouping.

**Which level is compared, and why.** All three depths are scored and the best is reported. The strategies work
blind and cannot know the granularity of the reference, and a nested map can be collapsed to any level by the
reader. Taking the best depth cannot reward a trivial level: a single L1 group is one of the baselines in `F0`.

**Direction and range.** Higher is better. `ref_skill` <= 1; 0 means "no better than the best trivial grouping";
negative means worse than that.

**Thresholds.** pass >= 0.25, warn 0.05 – 0.25, fail < 0.05. Judgement. I found no published threshold for BCubed
F1. Calibration with a *planted* reference (not the curated one, which I never saw): random groupings with realistic
sizes scored −0.58 to −0.10, the three trivial groupings <= 0 by construction, two sensible groupings that differ
from the planted one scored 0.00 to 0.68, exact recovery 1.00. The fail limit of 0.05 sits above everything trivial
or random; the pass limit of 0.25 is where the non-trivial controls begin to separate.

**Does not apply.** No reference for the model, or fewer than two labelled nodes: all `ref_*` are NaN. For the SKM
models only species are labelled, so this rule says nothing about where reactions are placed there.

**Gaming and counter-rule.** (a) Three levels are three attempts: bounded at three, and trivial levels are the
baseline. (b) Cloning a species into many groups raises recall; the precision term divides by the number of shared
groups, so surplus copies cost precision, and a pair that shares a group but no label scores 0. (c) Localisation-only
groupings score 0 by construction, so the localisation rules cannot be bought twice.

## LOC-1 — Localisation is expressed

**Curator's check.** In a model with several compartments, every node sits inside a box that is declared as its own
compartment.

**Metric.** `loc_expressed = |{ i : some g in chain(i) is a box whose space equals space(i) }| / N`.
Diagnostic `loc_pure_leaf` = share of nodes whose innermost group holds nodes of one space only (localisation that
is respected but not declared).

**Direction and range.** Higher, 0 to 1.

**Thresholds.** pass = 1.00, warn >= 0.90, fail < 0.90. The pass value follows from the contract ("localisation is
given, not chosen"): `nodes.tsv` places every node, so any strategy can reach 1.00 by construction, and anything
less is a choice. It also matches Process Description maps, where an entity pool node belongs to a single compartment
(SBGN PD user manual). The 0.90 warn limit is judgement.

**Does not apply.** Models with fewer than two populated spaces (NaN).

**Gaming and counter-rule.** Declaring every node its own compartment box gives 1.00. LOC-2 is the counter-rule: that
grouping has `loc_repeat` between 39.8 and 181.7 on the panel and fails it.

## BIO-2 — Members of a box are alike in annotation

**Curator's check.** Two annotated members of the same innermost box are more alike — same reaction type, enzyme
class, GO term, or pathway or bin code — than two annotated nodes picked from anywhere in the model.

**Metric.** Four channels of ontology terms per node, each closed under ancestors where a hierarchy is available:

| channel | nodes | terms |
|---|---|---|
| `sbo` | reactions | the SBO term and its ancestors (`sbo_path`) |
| `ec` | all | EC numbers as prefixes: `1`, `1.2`, `1.2.1`, `1.2.1.105` |
| `go` | all | GO ids, flat (the bundle carries no GO hierarchy) |
| `class` | all | cross-references with a class qualifier (`isVersionOf`, `isPartOf`, `occursIn`, `hasProperty`, `isPropertyOf`) other than EC and GO; hierarchical codes in the notes (`key: 14.5.7_text`) as prefixes |

```
IC(t)      = ln(N_ch / n_t) / ln(N_ch)        n_t nodes carry t; N_ch = reactions for sbo, all nodes otherwise
sim_ch(i,j)= max IC(t) over the terms i and j share in channel ch           (Resnik: most informative common ancestor)
sim(i,j)   = max over the channels in which both i and j have a term;  i, j are "comparable" if such a channel exists
S_in(i)    = mean sim(i,j) over the comparable j ≠ i in leaf(i);   if there is none, S_in(i) = S_all(i)
S_all(i)   = mean sim(i,j) over all comparable j ≠ i of the model
hom        = (mean_i S_in(i) − mean_i S_all(i)) / (1 − mean_i S_all(i))
```

The same value per channel is reported as `hom_sbo`, `hom_ec`, `hom_go`, `hom_class`; at L1 as `hom_L1`; for the
compartments alone as `hom_compartments`, and `hom_gain = hom − max(0, hom_compartments)`.

**Why this form.** Purity-type homogeneity is 1 for every singleton, so many small groups win for free; this is
why V-measure pairs its homogeneity term with completeness (Rosenberg and Hirschberg 2007), and why agreement
measures are corrected for chance (Vinh et al. 2010). Here a node with nobody to compare with gets exactly its chance value, so all-singletons and
one-group both score exactly 0. Sharing a term that most nodes carry is worth little because of the information
content weight (Pesquita et al. 2009): two reactions that are both merely "EC 1" are barely alike, two that share
`1.13.11` are. The species SBO term is left out on purpose: it is a material type, and a pathway box rightly mixes
metabolites and proteins.

**Level.** The innermost group of each node, because that is the box drawn directly around it. Upper levels are
meant to be broader.

**Direction and range.** Higher. At most 1; 0 is chance; negative when boxes mix unlike nodes more than chance.

**Thresholds.** pass >= 0.21, warn 0.05 – 0.21, fail < 0.05. The statistic has the form (observed − expected) /
(1 − expected) of a kappa coefficient; 0.21 is the lower end of the band Landis and Koch (1977) called "fair".
That scale is a convention, not a validated one, and the transfer to this statistic is by analogy: treat it as
judgement. The fail limit is calibration: random groupings with realistic sizes, 20 seeds per model, gave `hom`
between -0.032 and 0.047 on the six models.

**Does not apply.** Fewer than two comparable annotated nodes (NaN). A channel with no terms, or whose only term is
carried by every node (all 218 `brenda_merged` reactions are "biochemical reaction"), contributes nothing.

**Gaming and counter-rule.** (a) A grouping cut directly from one annotation scores high by construction (grouping
`brenda_merged` by EC sub-subclass: 0.62). BIO-1 is the counter-rule where a reference exists; where none exists
this is an open weakness (see *Unsure*). (b) Very small groups that each share a rare term score well; singletons
earn nothing, FAIR-4 penalises nodes separated from their reaction partners, and group-size rules belong to the
other viewpoints. (c) Coarse classes: grouping `brenda_merged` by EC class gives 0.40, by sub-subclass 0.62; grouping
the SKM models by SBO reaction type gives 0.12 to 0.18.

## LOC-2 — Each compartment is drawn as one box

**Curator's check.** Each compartment, and each compartment boundary, appears as one box; count the boxes that carry
the same compartment.

**Metric.** For every node that has a box, `b(i)` = its innermost box.

```
share(i)   = |{ j : b(j) = b(i) }| / |{ j : compartment_id(b(j)) = compartment_id(b(i)) }|
loc_repeat = 1 / mean_i share(i)
```

`mean share` is the chance that two nodes of the same compartment are drawn in the same box; its inverse is the
effective number of boxes per compartment (an inverse Simpson index, weighted by nodes). A compartment split into
`k` equal boxes gives `k`; a main box with a few stragglers stays close to 1. Raw counts are reported as
`loc_boxes_mean` and `loc_boxes_max`.

**Direction and range.** Lower is better, from 1 upwards.

**Thresholds.** pass = 1.00, warn <= 2.00, fail > 2.00. Judgement. One box per compartment is what the contract
promises the reader ("organelles inside cytosol, extracellular around it") and what a compartment container is in a
Process Description map. The limit of 2 means "at worst two equal halves".

**Does not apply.** Single-compartment models, and groupings with no box at all (NaN; LOC-1 fails instead).

**Gaming and counter-rule.** Declaring no localisation avoids this rule and fails LOC-1. The two rules must be read
together: LOC-1 asks that place is shown, LOC-2 that it is shown once.

## LOC-3 — Transport sits on the boundary it crosses

**Curator's check.** A reaction that moves something across a compartment boundary, and a species that sits in a
membrane, is drawn in the boundary box of that compartment and not among its interior reactions.

**Metric.** `B = { i : kind(i) = boundary }`.
`loc_boundary = |{ i in B : some g in chain(i) has compartment_id = space(i) + "@boundary" }| / |B|`.
Diagnostic `loc_boundary_in_interior` = share of `B` that sits in the interior box of its space and in no boundary box.
Nodes with `loc_kind = between` have no boundary to sit on; they belong to the interior of their `loc_space` and are
covered by LOC-1.

**Direction and range.** Higher, 0 to 1.

**Thresholds.** pass = 1.00, warn >= 0.80, fail < 0.80. Pass is 1.00 because `loc_kind` is given in `nodes.tsv`; the
warn limit is judgement. Transport drawn as a process between the pools of two compartments is the Process
Description convention (SBGN bricks dictionary).

**Does not apply.** No node with `loc_kind = boundary` (e.g. `brenda_merged`), or a single-compartment model.

**Gaming and counter-rule.** A boundary box can hold only boundary nodes (validator), so it cannot be padded. Making
one boundary box per transport reaction satisfies this rule and fails LOC-2.

## LOC-4 — Enclosing compartments are drawn around each node

**Curator's check.** Around every node, each compartment that encloses its own compartment is drawn exactly once in
the same part of the map: organelles inside cytosol, cytosol inside the extracellular space.

**Metric.** `A(i)` = the populated proper ancestors of `space(i)` in `compartments.tsv`. `scope(i)` = the parent group
of the outermost box of `i` (the top of the map if that box is at L1). An ancestor `a` is *drawn* for `i` if, among
the groups whose parent is `scope(i)`, exactly one interior box has compartment `a` (or no interior box and exactly
one boundary box).

```
loc_nesting = Σ_i |{ a in A(i) : a drawn for i }| / Σ_i |A(i)|          (a node without a box contributes 0)
```

**Direction and range.** Higher, 0 to 1.

**Thresholds.** pass = 1.00, warn >= 0.90, fail < 0.90. Judgement. The tree itself is not in question: SBML Level 3
dropped the `outside` attribute of compartments, so containment comes from `compartments.tsv`, the same for all.

**Does not apply.** Single-compartment models; models in which no node has a populated enclosing space.

**Assumption to check with the renderer.** The contract says localisation groups are drawn nested by
`compartments.tsv`. I assume this nesting happens among boxes that share a parent. A pathway group that holds a
peroxisome box and no cytosol box then shows a peroxisome with no cytosol around it. If the common renderer adds the
missing enclosing compartment itself, this rule is too strict for function-first hierarchies.

**Gaming and counter-rule.** An empty enclosing box cannot be added (a group needs a member). Splitting an enclosing
compartment in two makes the nesting ambiguous and counts as not drawn.

## FAIR-1 — Two runs give identical files

**Curator's check.** Running the strategy twice on the same model gives byte-identical `groups.tsv` and
`groups_meta.tsv`.

**Metric.** `fair_deterministic` = 1 if both files have the same SHA-256 in the two runs (CLI `--rerun`), else 0. In
the Python API, `rerun=(groups, meta)` compares the tables including row order. Diagnostics: `fair_rerun_same_rows`
(share of nodes with the same row), `fair_sorted` (rows sorted by id, as the contract asks).

**Protocol.** The scorer makes the second run itself, in a fresh process with a different `PYTHONHASHSEED`, which
exposes dependence on set or dict order. A third run on input tables with shuffled rows is a useful extra.

**Thresholds.** pass = 1, fail = 0; a fail sets the combined score of that model to 0. Source: contract rule 2;
recording random seeds is rule 6 of Sandve et al. (2013).

**Does not apply.** No second run supplied (NaN, reported as "na", not counted).

**Gaming and counter-rule.** Copying the first output; hence the scorer runs the strategy, not the strategy's author.

## FAIR-2 — A group's name is true of its members and tells them apart

**Curator's check.** The identifier behind a group's name is carried by the members of that group, and by clearly
fewer nodes outside it.

**Metric.** For a named group `g` with at least two members and `name_source` in sbo, ec, go, xref, notes,
compartment: `Car` = the nodes of the model that carry `name_xref` (SBO: the term is on the node's `sbo_path`; EC:
the code is a prefix of one of its EC numbers; GO: exact id; xref and notes: the `db:id` or `key:code` equals or is a
dotted prefix of one of the node's identifiers; compartment: the node is located there). `K` = the kinds (species,
reaction) of `Car`; `E_g` = the members of `g` of those kinds.

```
support_g = |E_g ∩ Car| / |E_g|
prev      = |Car| / |nodes of kinds K|
lift_g    = (support_g − prev) / (1 − prev)
name_lift = Σ_g |E_g| · lift_g / Σ_g |E_g|          over the groups with prev < 1
```

A name carried by every such node of the model (`prev = 1`: the compartment of a one-compartment model) is true
and says nothing; it is left out and counted in `name_vacuous_share`. A name whose identifier is on no member has
support 0. Groups named after a member are only checked for that member being in the group. Diagnostics:
`name_support`, `name_traceable`.

**Direction and range.** Higher; at most 1; 0 when the group is no richer in the term than the model; negative when
the name fits the rest of the model better than the group.

**Thresholds.** pass >= 0.61, warn 0.21 – 0.61, fail < 0.21: the "substantial" and "fair" bands of Landis and Koch
(1977), by analogy of form. Judgement.

**Does not apply.** No informative named group of two or more members (NaN; the combined score then counts it as 0).
GO names are matched by exact id only, because the bundle has no GO hierarchy.

**Gaming and counter-rule.** Naming one small perfect group and nothing else gives 1.00; FAIR-3 measures how many
nodes are under such a name. Names that hold for everything are excluded as vacuous.

## LOC-5 — A box spanning compartments holds their transport and is split by place

**Curator's check.** A group that holds nodes from several compartments also holds the transport reactions that
connect those compartments, and is divided into compartment boxes one level down.

**Metric.** A group `g` *spans* if `members(g)` have two or more `loc_space` values. It is *linked* if joining, for
every member reaction, the spaces of its reactants and products (a membrane counts as the space it encloses) and its
own `loc_space`, connects all spaces of `g`.

```
bad(i)      = i lies in a spanning group that is not linked, or lies in a spanning group and has no box
loc_span_ok = 1 − |{ i : bad(i) }| / N
```

Diagnostics: `loc_span_share` (share of nodes in spanning groups), `loc_span_connected`.

**This is the answer to "when is a functional group across compartments acceptable".** Three conditions: its own
transport steps connect its compartments (this rule); every member still sits in a compartment box (this rule and
LOC-1); the repetition of compartments that results stays within LOC-2. A pathway such as jasmonate synthesis, which
runs from chloroplast through peroxisome to cytosol, meets the first two when its transport reactions are members.

**Direction and range.** Higher, 0 to 1. 1.00 also when no group spans: the rule forbids unjustified spans, it does
not ask for spans. What a span gains shows in BIO-1.

**Thresholds.** pass = 1.00, warn >= 0.90, fail < 0.90. Judgement.

**Does not apply.** Single-compartment models.

**Gaming and counter-rule.** Never spanning is allowed and scores 1.00. One group for the whole model scores 0.00 on
every multi-compartment panel model.

## FAIR-3 — Most nodes sit under a name backed by a resolvable identifier

**Curator's check.** Pick any node: at least one of the boxes around it has a name whose identifier can be looked up
outside the model — an SBO, EC or GO term, a database cross-reference, or a compartment mapped to GO.

**Metric.** A group is *well named* if it is named, has at least two members, its identifier is carried by a
member, is not vacuous (FAIR-2), and resolves: `SBO:nnnnnnn`, `GO:nnnnnnn`, an EC code, a cross-reference whose
prefix is a database used in the model's annotation, or a compartment with a `go_id` in `compartments.tsv`. It is
*locally named* if all of that holds except resolution (named after a member node, or after a notes key).

```
name_grounded = share of nodes with a well named group in their chain
                + 0.5 · share of nodes with only locally named groups in their chain
```

Diagnostics: `name_resolvable` (the first term), `name_coverage` (any name at all), `named_group_share`.

**Direction and range.** Higher, 0 to 1.

**Thresholds.** pass >= 0.80, warn 0.50 – 0.80, fail < 0.50. Judgement. The half credit for local names is
judgement: the contract allows them, and they are honest, but they cannot be looked up elsewhere. Background: the FAIR
principles ask for persistent identifiers and for vocabularies that are themselves FAIR (Wilkinson et al. 2016); SBML
groups can be annotated to say what they represent (Hucka and Smith 2016); `contract.py write-sbml` writes such an
annotation only for SBO, GO and EC identifiers.

**Does not apply.** Always applies; 0 when nothing is named.

**Gaming and counter-rule.** Naming every box after its compartment gives 1.00 and is in fact well grounded; BIO-1
and BIO-2 show that such a grouping says nothing about function.

## FAIR-4 — A node without annotation is still placed with a reaction partner

**Curator's check.** A species or reaction that has no ontology annotation is in the same box as something it reacts
with, not alone and not in a rest group of unrelated leftovers.

**Metric.** `U` = nodes with no term in any BIO-2 channel and at least one neighbour in `edges.tsv`. Two nodes are
*together* if they have the same innermost group, or a `clone_in` group of one lies in the chain of the other.

```
stranded_poor = |{ i in U : no neighbour of i is together with i }| / |U|
```

Diagnostics: `stranded_annotated` (the same for annotated nodes), `poor_share` (= |nodes without annotation| / N).

**Direction and range.** Lower, 0 to 1.

**Thresholds.** pass <= 0.10, warn <= 0.30, fail > 0.30. Judgement, from calibration: random groupings 0.58 to 0.94,
all-singletons 1.00, the non-trivial controls 0.00 to 0.33.

**Does not apply.** No node without annotation that has a neighbour (NaN).

**Gaming and counter-rule.** One group for everything gives 0.00; BIO-1, BIO-2 and LOC-1 fail it.

**Stronger test, not implemented here.** Blank the annotation columns of a fixed 30 % of nodes (chosen by sorting
node ids, so it is reproducible), rerun the strategy, and require that it still validates and that the grouping of
the untouched nodes changes little. This needs the strategies to be rerun and is left to the scoring script.

---

## Combining the rules into one score

`metrics_biology.score(m)` implements this.

1. **Gate.** FAIR-1 = 0 sets the score of that model to 0. A grouping that fails the contract validator is not scored.
2. **Utility of each rule.** The metric is mapped to 0..1 by straight lines through four points: 0 at its worst
   value, 0.25 at the warn threshold, 0.75 at the pass threshold, 1 at its best value. A rule that passes therefore
   earns at least 0.75, one that only warns between 0.25 and 0.75. `loc_repeat` is first turned into `1 / loc_repeat`.
3. **Blocks.** Rules are averaged inside four blocks with these weights; a rule that does not apply is left out and
   the remaining weights are rescaled.

   | block | rules and weights inside the block | block weight |
   |---|---|---|
   | reference | BIO-1 (1.0) | 0.35 |
   | homogeneity | BIO-2 (1.0) | 0.25 |
   | localisation | LOC-1 0.25, LOC-2 0.25, LOC-3 0.20, LOC-4 0.15, LOC-5 0.15 | 0.25 |
   | fair | FAIR-2 0.40, FAIR-3 0.30, FAIR-4 0.30 | 0.15 |

4. **Score.** Weighted mean of the blocks that apply (no reference: the reference block is left out; one compartment:
   the localisation block is left out). `name_lift` counts as 0 when no group carries an informative name.

The weights are judgement. The reasoning: the held-out reference is the only evidence that is independent of the
annotation a strategy may have grouped by, so it weighs most; localisation and homogeneity are equal; the FAIR
rules are close to pass/fail for a careful strategy and separate little, so they weigh least.

Advice for the overall ranking across the three viewpoints:

- Report the eleven verdicts side by side per strategy and model; use the score only to order them.
- Average over models, not over nodes, so that `skm_pss` does not decide everything; show models with and without a
  reference separately, because the score means different things in the two cases (next section).
- Read LOC-1 and LOC-2 as a pair, BIO-2 together with `hom_gain`, and FAIR-2 together with FAIR-3. Each first rule
  has a cheap way to be satisfied which the second one closes.
- Group count, group size and singleton share are not judged here. They belong to the clarity viewpoint; my rules
  only make sure singletons earn nothing.

## Behaviour on synthetic groupings

All six panel models; no model dropped or subsampled. Groupings: one group; every node its own group; every node
its own compartment box; random with the group sizes of a Louvain partition; one box per compartment (with boundary
boxes); compartments then Louvain communities; Louvain communities then compartment boxes; compartments then
reaction type (EC sub-subclass, else SBO term); and three one-level groupings by annotation. The reference used here
is **planted** from Louvain communities with multi-label hub species and 30 % of nodes unlabelled (species only in
the SKM models); "planted" is the grouping that equals it. It stands in for the curated reference, which I did not
see, so the reference columns show how the metric behaves, not how any real strategy will score.

Combined score with the planted reference:

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | 0.06 | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| singletons | 0.00 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| singletons_as_localisation | 0.00 | 0.10 | 0.10 | 0.10 | 0.10 | 0.10 |
| random | 0.01 | 0.04 | 0.01 | 0.00 | 0.00 | 0.00 |
| compartments | 0.06 | 0.56 | 0.44 | 0.43 | 0.45 | 0.44 |
| compartments>communities | 0.61 | 0.55 | 0.66 | 0.77 | 0.69 | 0.71 |
| communities>compartments | 0.60 | 0.80 | 0.77 | 0.71 | 0.72 | 0.72 |
| compartments>type | 0.49 | 0.56 | 0.81 | 0.66 | 0.76 | 0.70 |
| planted | 0.65 | 0.48 | 0.62 | 0.48 | 0.49 | 0.48 |
| by_sbo_type | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_class | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_subsubclass | n/a | n/a | n/a | n/a | n/a | n/a |

Combined score without a reference:

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | 0.11 | 0.07 | 0.07 | 0.07 | 0.07 | 0.07 |
| singletons | 0.00 | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| singletons_as_localisation | 0.00 | 0.15 | 0.15 | 0.15 | 0.16 | 0.15 |
| random | 0.02 | 0.06 | 0.01 | 0.00 | 0.00 | 0.01 |
| compartments | 0.11 | 0.86 | 0.68 | 0.66 | 0.70 | 0.68 |
| compartments>communities | 0.37 | 0.80 | 0.84 | 0.75 | 0.69 | 0.71 |
| communities>compartments | 0.34 | 0.76 | 0.73 | 0.57 | 0.59 | 0.60 |
| compartments>type | 0.92 | 0.87 | 0.92 | 0.78 | 0.83 | 0.79 |
| planted | 0.34 | 0.21 | 0.42 | 0.20 | 0.21 | 0.20 |
| by_sbo_type | 0.11 | 0.50 | 0.07 | 0.41 | 0.48 | 0.43 |
| by_ec_class | 0.88 | 0.50 | 0.57 | 0.11 | 0.07 | 0.12 |
| by_ec_subsubclass | 0.92 | 0.50 | 0.56 | 0.13 | 0.07 | 0.15 |

`hom` (BIO-2):

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| singletons | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| singletons_as_localisation | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| random | 0.00 | 0.02 | -0.02 | -0.00 | -0.00 | 0.00 |
| compartments | 0.00 | 0.18 | 0.04 | 0.03 | 0.07 | 0.05 |
| compartments>communities | 0.09 | 0.12 | 0.16 | 0.10 | 0.06 | 0.07 |
| communities>compartments | 0.09 | 0.22 | 0.19 | 0.11 | 0.13 | 0.11 |
| compartments>type | 0.62 | 0.38 | 0.40 | 0.13 | 0.18 | 0.14 |
| planted | 0.09 | 0.08 | 0.20 | 0.07 | 0.09 | 0.08 |
| by_sbo_type | 0.00 | 0.18 | 0.00 | 0.12 | 0.18 | 0.14 |
| by_ec_class | 0.40 | 0.29 | 0.50 | 0.00 | 0.00 | 0.00 |
| by_ec_subsubclass | 0.62 | 0.31 | 0.44 | 0.00 | 0.00 | 0.00 |

`ref_skill` (BIO-1, planted reference):

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | 0.00 | -0.03 | -0.32 | -0.10 | -0.29 | -0.24 |
| singletons | -0.47 | -0.40 | -0.74 | -0.09 | -0.26 | -0.24 |
| singletons_as_localisation | -0.47 | -0.40 | -0.74 | -0.09 | -0.26 | -0.24 |
| random | -0.24 | -0.25 | -0.58 | -0.10 | -0.26 | -0.23 |
| compartments | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| compartments>communities | 0.68 | 0.02 | 0.09 | 0.43 | 0.23 | 0.23 |
| communities>compartments | 0.68 | 0.66 | 0.59 | 0.89 | 0.86 | 0.81 |
| compartments>type | 0.00 | 0.00 | 0.19 | 0.12 | 0.20 | 0.16 |
| planted | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| by_sbo_type | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_class | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_subsubclass | n/a | n/a | n/a | n/a | n/a | n/a |

`loc_repeat` (LOC-2), `loc_nesting` (LOC-4):

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | n/a | n/a | n/a | n/a | n/a | n/a |
| singletons | n/a | n/a | n/a | n/a | n/a | n/a |
| singletons_as_localisation | n/a | 83.50 | 49.33 | 181.70 | 39.80 | 58.67 |
| random | n/a | n/a | n/a | n/a | n/a | n/a |
| compartments | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| compartments>communities | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| communities>compartments | n/a | 5.89 | 4.71 | 12.36 | 3.81 | 5.95 |
| compartments>type | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| planted | n/a | n/a | n/a | n/a | n/a | n/a |
| by_sbo_type | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_class | n/a | n/a | n/a | n/a | n/a | n/a |
| by_ec_subsubclass | n/a | n/a | n/a | n/a | n/a | n/a |

| grouping | brenda_merged | e_coli_core | R-HSA-77289 | skm_pss | skm_ja | skm_sa |
|---|---|---|---|---|---|---|
| one_group | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| singletons | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| singletons_as_localisation | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| random | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| compartments | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| compartments>communities | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| communities>compartments | n/a | 0.76 | 0.00 | 0.17 | 0.85 | 1.00 |
| compartments>type | n/a | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| planted | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| by_sbo_type | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| by_ec_class | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| by_ec_subsubclass | n/a | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

What the tables show:

- One group, all singletons, singletons declared as compartment boxes, and random never exceed 0.10 with a reference
  and 0.16 without; the weakest non-trivial control is 0.49 and 0.34. Neither degenerate grouping comes out best on
  any model. This is asserted in the tests.
- One box per compartment is not degenerate, and it is the grouping any strategy gets for free. With a reference it
  scores 0.43 to 0.56 (0.06 on the one-compartment model). Without a reference it scores 0.66 to 0.86 (0.11 on the
  one-compartment model), because it satisfies every localisation and FAIR rule and only BIO-2 is left to tell it
  from finer groupings.
- Function-first hierarchies (communities, then compartment boxes) pay on LOC-2 and LOC-4 and gain on BIO-1.

## What I am unsure about

1. **How the renderer nests boxes.** LOC-2 and LOC-4 assume that compartment boxes are nested only among siblings.
   If the renderer draws one cell and places function groups inside its compartments, function-first hierarchies are
   penalised too much and both metrics should be recomputed on what is actually drawn.
2. **BIO-2 rewards grouping by the annotation itself**, and without a reference nothing in my set contradicts it. On
   the SKM models reactions carry only an SBO type; the best value I could produce was 0.18, below the pass limit, so
   0.21 may be out of reach there for honest reasons. A per-model limit, or `hom_gain` > 0, may serve better.
3. **Without a reference the score mostly measures compliance.** One box per compartment reaches up to 0.86. For
   models without a reference, granularity has to be judged by the other viewpoints.
4. **The compartment baseline inside `ref_skill`.** It stops localisation from being rewarded twice, and it would be
   too harsh if the curated labels themselves follow compartments. `ref_f1` and the three baseline values are
   reported so that this can be seen.
5. **Best-of-three levels in BIO-1.** Fair to blind strategies, lenient to hierarchies whose other levels are poor.
   `ref_f1_L1`, `ref_f1_L2`, `ref_f1_L3` are reported for a stricter reading.
6. **No GO hierarchy in the bundle.** GO terms are compared and matched by exact id; a group named by a GO ancestor
   of its members' terms counts as untraceable.
7. **The `class` channel is a heuristic.** It reads class-type qualifiers and dotted numeric codes in the notes. It
   fits the MapMan bins and KEGG pathway ids of the SKM models; other models may keep such information elsewhere.
8. **Thresholds.** None comes from a study of pathway maps. The 0.21 and 0.61 limits borrow the Landis–Koch bands by
   analogy, and that scale is itself only a convention. The rest is judgement, checked against synthetic groupings.
9. **FAIR-4 is a single-run proxy.** The ablation test described under FAIR-4 is the real test of behaviour under
   missing annotation.
10. **Species-only labels.** On the SKM models BIO-1 cannot see reactions; a strategy that places species well and
    reactions badly is not caught by it.

## Deviations from the brief

- Eleven rules, not twelve.
- `metrics()` has one optional argument beyond the requested signature, `rerun=(groups, meta)`, because determinism
  cannot be read from one run. Without it FAIR-1 is NaN.
- The annotation-ablation test for missing annotation is described, not implemented.
- SBO homogeneity uses reactions only; GO is flat. Both are stated under BIO-2.
- The curated reference was not available to me. BIO-1 was tested on hand-made cases with known values and on a
  planted reference; its behaviour on the real reference is untested.
- Literature was read through search results (titles, abstracts, excerpts), not full texts.
- Part 2 of the tests needs the bundle (`GROUPING_PANEL` or a `panel/features` folder nearby) and networkx; it is
  skipped with a message if the panel is not found.

## Sources

Only records retrieved by web search during this task are listed.

- Amigó E, Gonzalo J, Artiles J, Verdejo F (2009). A comparison of extrinsic clustering evaluation metrics based on
  formal constraints. Information Retrieval 12:461–486. doi:10.1007/s10791-008-9066-8 — BCubed and its extension to
  overlapping clusters (BIO-1).
- Rosenberg A, Hirschberg J (2007). V-Measure: a conditional entropy-based external cluster evaluation measure.
  EMNLP-CoNLL 2007, 410–420. https://aclweb.org/anthology/D/D07/D07-1043.pdf — homogeneity and completeness (BIO-2).
- Vinh NX, Epps J, Bailey J (2010). Information theoretic measures for clusterings comparison: variants, properties,
  normalization and correction for chance. JMLR 11:2837–2854.
  https://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf — why agreement measures need a chance correction.
- Pesquita C, Faria D, Falcão AO, Lord P, Couto FM (2009). Semantic similarity in biomedical ontologies. PLoS Comput
  Biol 5(7):e1000443. doi:10.1371/journal.pcbi.1000443 — information content and Resnik similarity (BIO-2).
- Landis JR, Koch GG (1977). The measurement of observer agreement for categorical data. Biometrics 33(1):159–174.
  Bands as reported at https://en.wikipedia.org/wiki/Cohen's_kappa, which also notes that they rest on opinion.
- Wilkinson MD et al. (2016). The FAIR Guiding Principles for scientific data management and stewardship. Scientific
  Data 3:160018. doi:10.1038/sdata.2016.18 — identifiers, vocabularies, provenance (FAIR-2, FAIR-3).
- Sandve GK, Nekrutenko A, Taylor J, Hovig E (2013). Ten simple rules for reproducible computational research. PLoS
  Comput Biol 9(10):e1003285. doi:10.1371/journal.pcbi.1003285 — rule 6, random seeds (FAIR-1).
- Hucka M, Smith LP (2016). SBML Level 3 package: Groups, Version 1 Release 1. J Integr Bioinform 13(3):8–29.
  doi:10.1515/jib-2016-290 — nested, annotated groups (FAIR-3).
- Rougny A et al. (2019). Systems Biology Graphical Notation: Process Description language Level 1 Version 2.0.
  J Integr Bioinform. doi:10.1515/jib-2019-0022; and the SBGN PD Level 1 user manual,
  https://sbgn.github.io/downloads/usermanual/sbgn_PD-level1-user-public.pdf — an entity pool node belongs to a
  single compartment (LOC-1, LOC-2).
- SBGN bricks dictionary, https://sbgnbricks.sourceforge.net/sbgnbricks_dictionary.html — transport as a process
  between pools in two compartments (LOC-3).
- libSBML API documentation, class Compartment,
  https://sbml.org/software/libsbml/5.18.0/docs/formatted/python-api/classlibsbml_1_1_compartment.html — the
  `outside` attribute exists in SBML Level 1 and 2 and not in Level 3 (LOC-4).
- BioModels.net qualifiers, COMBINE specification,
  https://github.com/combine-org/combine-specifications/blob/main/specifications/qualifiers-1.1.md — meaning of
  `isVersionOf`, `isPartOf`, `occursIn` (the `class` channel of BIO-2).
