# Graph rules for judging a grouping (viewpoint: graph)

BioHackathon Europe 2026, project 22. Twelve rules, most important first, computed by `metrics_graph.py`
from the contract tables only (`nodes.tsv`, `edges.tsv`, `groups.tsv`; `compartments.tsv` and
`groups_meta.tsv` are accepted but not needed). Written blind: no strategy output was seen. Thresholds
were checked against synthetic groupings built here (section 5).

Run:

```
python metrics_graph.py score --features panel/features/skm_ja --grouping strategies/NAME/out/skm_ja
python metrics_graph.py score --features panel/features/skm_ja --grouping strategies/NAME/out/skm_ja --reference reference/skm_ja.tsv --json skm_ja_graph.json
python metrics_graph.py demo  --features panel/features/skm_ja
python -m pytest test_metrics_graph.py -q
```

## 1. Definitions used by every rule

**Graph.** Bipartite: species S on one side, reactions R on the other, V = S ∪ R, N = |V|. One edge per
distinct (reaction, species) pair of `edges.tsv`; reactant, product and modifier count alike; a species
with two roles in one reaction is one edge. d_s = number of reactions of species s.

**Hub treatment (data-derived, two parts).**

1. *Species-normalised weight.* Edge (r, s) weighs w_rs = 1 / d_s. Every species then has total weight 1
   ("one vote"), whether it takes part in 1 reaction or in 65; W = Σ w = number of species with an edge.
   The strength of a node is k_v = Σ of w over its edges (species: 1; reaction: Σ_s 1/d_s). All `_w`
   metrics use this weight. The `_raw` twins (every edge 1) are returned as diagnostics, never judged.
   Precedent: weighting a tie inversely to the size of the set it comes from (Newman 2001), and
   weighting compounds by their degree instead of deleting a hand-made list of them (Croes et al. 2006).
2. *Hub set H.* A species is a hub when it has significantly more reactions than uniform attachment
   gives: d_s > c, with c the smallest integer such that P(Poisson(λ) > c) ≤ 0.05 / |S| and λ the mean
   degree of the non-hub species (iterated until H stops changing). H decides which copies (`clone_in`)
   are honoured (G1, G7, G11) and feeds one diagnostic (`connected_nohub_node_share`). It replaces a list
   of molecule names, which the contract forbids. Precedent for deriving the set from the network itself:
   Huss & Holme 2007 (as described in "Pathway identification by network pruning ...", section 7:
   remove the highest-degree metabolites until modularity peaks); I use a closed-form degree test instead because the scorer must not depend on a
   clustering run.

What the two parts do on the panel (computed):

| model | species | hub cutoff c | hubs | share of all edges on hubs | hubs found (name, first 12 characters) |
|---|---|---|---|---|---|
| brenda_merged | 178 | 11 | 16 | 0.43 | H2O (65), ATP (43), ADP (42), phosphate (36), NH3 (24), NADH (21), NAD+ (21), H+ (20) ... |
| e_coli_core | 72 | 12 | 4 | 0.23 | H+ (35), H2O H2O (18), H+ (17), ATP C10H12N5 (13) |
| R-HSA-77289 | 101 | 7 | 9 | 0.42 | H+ [mitochon (20), CoA-SH [mito (13), FADH2 [mitoc (12), Trifunctiona (12), FAD [mitocho (12), H2O [mitocho (10), NAD+ [mitoch (10), NADH [mitoch (10) ... |
| skm_pss | 961 | 10 | 14 | 0.13 | VPg (53), SARD1 (31), CBP60G (30), CDPK (25), MYC2 (21), Ca2+ (20), SNRK2 (18), HC-Pro (18) ... |
| skm_ja | 117 | 7 | 3 | 0.17 | MYC2 (21), JAZ (9), JAM (8) |
| skm_sa | 193 | 8 | 4 | 0.17 | SARD1 (31), SA (16), ICS (13), NPR1 (10) |

With raw counts the hubs carry 13-43 % of all edges; with the weight each hub carries 1/|S| of the total.
On the optimiser's partition the raw cut is 0.13-0.41 and the weighted cut 0.03-0.19.

**Levels.** P_k(v) = the group of v at level k, or v's deepest group above k if v has none at level k
(a node placed directly in a parent that also has child groups forms a *residual group* with the parent's
id; `residual_node_share` reports how many nodes). "Leaf" = the deepest level used anywhere. Every
per-level metric is returned with suffixes `_L1`, `_L2`, `_L3`, `_leaf`; an unused level is NaN.

**Per group g at a level.** n_g = number of home members; in(g) = weight of edges with both ends in g;
KR_g, KS_g = summed strength of its reactions and of its species; vol(g) = KR_g + KS_g;
cut(g) = vol(g) - 2·in(g).

**Clones.** A copy (`clone_in`) changes only G1 and G7, and only for hub species: edge (r, s) is local at
level k if s ∈ H and a copy of s lies in P_k(r). Copies of other species are drawn but never make an edge
local. All other rules use the home position, so cloning cannot improve them.

**A level with one group separates nothing.** Cut, modularity, conductance are NaN there. If the leaf
level has one group, G1, G5 and G6 are **failed, not skipped**.

## 2. The rules

Verdicts: pass / warn / fail / na. "Judgement" marks a threshold that is mine, not from a source.

### G1 - Edges stay inside boxes (cut share at the leaf)
- **Curator check.** Pick any species: most of the reactions it takes part in are drawn in its own box.
- **Metric** `cut_w_leaf` = (1/W) Σ_(r,s) w_rs · [P_leaf(r) ≠ P_leaf(s) and no honoured copy of s in
  P_leaf(r)]. Equivalently: the mean, over species, of the share of its reactions drawn outside its box.
- **Direction, range.** Lower is better; [0, 1].
- **Thresholds.** pass ≤ 0.25; warn ≤ 0.50; fail > 0.50. The 0.50 bound follows from Radicchi et al.
  2004: a group is a community in the weak sense when its internal degree sum exceeds its external one;
  summed over groups that is Σ 2·in(g) > Σ cut(g), i.e. cut share < 0.5. Above it the grouping as a whole
  is not a set of communities. The 0.25 bound is judgement, placed between what random groupings (0.92-0.99)
  and module-based groupings (0.02-0.33) reach on the panel.
- **Levels.** Non-decreasing with depth (an edge cut at level k is cut at every deeper level):
  `cut_w_L1` ≤ `cut_w_L2` ≤ `cut_w_L3`. Judged at the leaf; the levels are combined in G7, not here.
- **Does not apply.** Model without edges (na). One group: fail.
- **Gaming.** (a) Few huge boxes cut nothing: one box per compartment reaches 0.00-0.17. Counter-rules G2,
  G8. (b) Copies of everything: blocked inside the metric (only hub copies count) and judged by G11.

### G2 - No dominant box (largest group's share) - gate
- **Curator check.** No single box holds more than a third of the map.
- **Metric** `largest_share_leaf` = max_g n_g / N over leaf groups.
- **Direction, range.** Lower is better; (0, 1].
- **Thresholds.** pass ≤ 1/3; warn ≤ 1/2; fail > 1/2. Judgement. Anchor from the brief: a 42-node model
  drawn as one box (share 1.0) is the stated example of an unclear result.
- **Levels.** Non-increasing with depth. Judged at the leaf only: at L1 a localisation group can
  legitimately hold most of the model (cytosol), and localisation is given, not chosen.
- **Does not apply.** Never. For a model too small to split (N < 9) read the verdict with care.
- **Gaming.** Splitting into crumbs passes: counter-rules G3, G8, G9.

### G3 - Few single-node boxes (singleton share) - gate
- **Curator check.** Hardly any box holds a single node.
- **Metric** `singleton_group_share_leaf` = #{leaf groups with n_g = 1} / #{leaf groups}. A node left
  alone in a parent beside child groups counts as a singleton. `singleton_node_share_leaf` is the same
  per node.
- **Direction, range.** Lower is better; [0, 1].
- **Thresholds.** pass ≤ 0.05; warn ≤ 0.20; fail > 0.20. Judgement. Anchor from the brief: the current
  pipeline's 167 singletons of 346 groups (0.48) is the stated example of an unclear result.
- **Levels.** Reported per level; judged at the leaf.
- **Does not apply.** Never.
- **Gaming.** Sweeping singletons into one "rest" box passes: that box is disconnected and leaky, caught
  by G4 and G6.

### G4 - Every box is one connected piece
- **Curator check.** Inside a box you can walk from any node to any other along edges without leaving it.
- **Metric** `connected_node_share_leaf` = Σ_g |largest connected component of the subgraph induced by
  the home members of g| / Σ_g n_g, over leaf groups with n_g ≥ 2. `disconnected_group_share_leaf` =
  share of those groups in more than one piece.
- **Direction, range.** Higher is better; (0, 1].
- **Thresholds.** pass ≥ 0.95; warn ≥ 0.80; fail < 0.80. The requirement itself (a community must be
  internally connected) is from Traag, Waltman & van Eck 2019, who treat disconnected communities as a
  defect of the Louvain algorithm and report up to 16 % of them in practice; the tolerances 0.95 and 0.80
  are judgement.
- **Levels.** Usually rises toward L1 (unions of linked groups), not guaranteed. Judged at the leaf.
- **Does not apply.** No group with two nodes (na; G3 then fails).
- **Gaming.** (a) Singletons are trivially connected: excluded from the metric, caught by G3. (b) One box
  is connected whenever the model is: caught by G2. (c) A bag of unrelated reactions held together only
  by a currency species that lives there: see the diagnostic `connected_nohub_node_share_leaf` (hub
  species removed before the test). It is not judged, because in signalling models the hub is often the
  legitimate centre of the group (a transcription factor with its targets). Copies never connect a group.

### G5 - Boxes follow the network's own clusters (modularity)
- **Curator check.** The boxes agree with the clusters a dedicated clustering algorithm finds.
- **Metric.** Barber's bipartite modularity with the species-normalised weight, per level k:
  Q_k = Σ_g [ in(g)/W - KR_g · KS_g / W² ]. `Q_w_best` = max_k Q_k over levels with ≥ 2 groups.
  `Q_ref` = Q of the partition found by a deterministic Louvain-type optimiser of the same function
  (in `metrics_graph.py`). **`Q_eff` = `Q_w_best` / `Q_ref`** is judged.
  *Null model:* the bipartite configuration model - the expected weight between reaction r and species s
  is k_r·k_s/W, and nothing is expected between two species or two reactions (Barber 2007). *Hubs:*
  handled by the weight, every species has strength 1. *Clones:* ignored.
- **Direction, range.** Higher is better. Q < 1; 0 for a random grouping and for groups holding one kind
  of node only. Q_eff ≈ 0 (random) to ≈ 1; it can pass 1 slightly if a strategy beats the greedy optimiser.
- **Thresholds.** pass ≥ 0.70; warn ≥ 0.40; fail < 0.40. Judgement, checked on the panel (random
  -0.02-0.02; one box per compartment 0.25-0.60; compartment × modules 0.83-1.00). A relative threshold is
  used because absolute modularity is not comparable across networks: random graphs already show sizeable
  modularity from fluctuations (Guimerà, Sales-Pardo & Amaral 2004). The customary absolute reading
  (0.3-0.7 for networks with clear community structure, attributed to Newman & Girvan 2004) applies to
  unweighted modularity; `Q_raw_best` is returned for that comparison.
- **Levels.** Not monotone. Modularity has an intrinsic scale and cannot resolve modules below it
  (Fortunato & Barthélemy 2007), so a fine leaf level scores lower although it may be good; hence the
  best level is taken, not the leaf.
- **Does not apply.** Model without edges (na). One group: fail.
- **Gaming.** High modularity with unreadable sizes or disconnected groups: G2, G8, G4. Adding levels in
  the hope that one scores well is limited to three levels and costs G10 when a level is redundant.

### G6 - Each box is a community (conductance)
- **Curator check.** Each box has more edge weight inside it than leaving it.
- **Metric.** For a leaf group, conductance φ(g) = cut(g) / min(vol(g), 2W - vol(g)); g is a *weak
  community* when 2·in(g) > cut(g) (for vol(g) ≤ W this is φ(g) < 1/2).
  `weak_node_share_leaf` = Σ_{g weak} n_g / N. `conductance_mean_leaf` (node-weighted mean φ) is returned.
- **Direction, range.** Higher is better; [0, 1].
- **Thresholds.** pass ≥ 0.80; warn ≥ 0.50; fail < 0.50. The per-group criterion is the weak-community
  definition of Radicchi et al. 2004; conductance is one of the two structural scores that best recover
  ground-truth communities in Yang & Leskovec 2012. The shares 0.80 and 0.50 are judgement.
- **Levels.** Rises toward L1 (larger groups keep more inside). Judged at the leaf.
- **Does not apply.** Model without edges (na). One group: fail.
- **Gaming.** Huge boxes are always weak communities (one box per compartment: 1.00): G2, G8. Differs
  from G1 in that G1 is an average and hides individual leaky boxes.

### G7 - Crossing edges stay close; boxes are small against the map (edge span over the hierarchy)
- **Curator check.** An edge that leaves its box ends in a sibling box under the same parent, not on the
  far side of the map.
- **Metric.** Dasgupta's cost of the group tree, normalised. box(e) = the smallest group holding both
  ends of edge e (an honoured hub copy counts as an end); the whole model if none.
  `edge_span_w` = Σ_e w_e · |box(e)| / (N · W). `edge_span_null` = the same with every (reaction, species)
  pair as an edge, i.e. nodes dealt at random into the same boxes. **`edge_span_ratio`** = their quotient.
- **Direction, range.** Lower is better; (0, ≈1]. 1 = no better than a random filling.
- **Thresholds.** pass ≤ 1/3; warn ≤ 2/3; fail > 2/3. Judgement, checked on the panel (module-based
  0.04-0.42; one box per compartment 0.51-1.00; random 0.99-1.02). The cost function is from Dasgupta 2016.
- **Levels.** This is how the levels are combined: one number for the whole tree; a crossing costs the
  size of the lowest box that still holds both ends. A dummy level gains nothing.
- **Does not apply.** Model without edges.
- **Gaming.** None from the two degenerate ends: one box and all singletons both give exactly 1.00 without
  any counter-rule. Many tiny tight boxes (one per reaction: 0.27-0.68) can pass: G3, G8, G9. The ratio is
  easier to lower in large models, where any box is small against N (see section 6).

### G8 - Boxes have a readable size
- **Curator check.** A box holds between 3 and 50 nodes.
- **Metric** `readable_node_share_leaf` = Σ_{g: 3 ≤ n_g ≤ 50} n_g / N over leaf groups.
- **Direction, range.** Higher is better; [0, 1].
- **Thresholds.** pass ≥ 0.80; warn ≥ 0.50; fail < 0.50. Upper size 50: Yoghourdjian et al. 2021 found
  that readers struggle with node-link diagrams above 50 nodes at high density and above 100 at low
  density; Ghoniem, Fekete & Castagliola 2004 found node-link diagrams lose to matrices on most tasks
  above 20 nodes; 50 is the lenient end of that range. Lower size 3 (a reaction with a substrate and a
  product) and the shares 0.80 / 0.50 are judgement.
- **Levels.** Judged at the leaf; reported per level.
- **Does not apply.** N < 3.
- **Gaming.** Random groups of ten pass (0.81-1.00): G1, G4-G7 fail them.

### G9 - A handful of boxes side by side at every level
- **Curator check.** At every level there are at least 2 and at most 9 boxes side by side.
- **Metric** `max_siblings` = the largest sibling set: the top-level groups, and for every group its child
  groups (plus one if it also holds nodes directly). Also `n_groups_L1`.
- **Direction, range.** Lower is better; [1, N]. Exactly one top-level group fails.
- **Thresholds.** pass 2-9; warn 10-50; fail > 50 or a single top-level group. 9 is the upper end of
  Miller's 7 ± 2 (Miller 1956); 50 treats the overview as a node-link diagram of boxes (Yoghourdjian et
  al. 2021). Applying a memory-span figure to boxes on a map is judgement.
- **Levels.** Evaluated over all levels at once.
- **Does not apply.** Never.
- **Gaming.** Wrapping boxes in extra boxes to shorten sibling lists: limited by the three levels and
  judged by G10. Intended effect: a large model cannot pass with one flat level.

### G10 - No box that holds just one box
- **Curator check.** No box contains exactly one box and nothing else.
- **Metric** `single_child_share` = share of parent groups with exactly one child group and no node of
  their own; a nested grouping with a single top-level group adds the undrawn root as such a parent.
- **Direction, range.** Lower is better; [0, 1].
- **Thresholds.** pass = 0; warn ≤ 0.20; fail > 0.20. Judgement. Always avoidable: L2 and L3 are optional
  per node, so a group that would have one child simply has none.
- **Does not apply.** Groupings with one level (na).
- **Gaming.** None; this is the counter-rule for dummy levels (G5, G7, G9).

### G11 - Copies are spent on hubs and every copy is used
- **Curator check.** Only species used all over the map are drawn more than once, and every copy has one
  of its reactions in the same box.
- **Metric** `clone_misuse_share` = share of `clone_in` entries (species, group) that are idle (no reaction
  of the species lies in that group) or belong to a species outside H. Also `clone_idle_share`,
  `clone_nonhub_share`, `clones_per_species`, `clone_absorbed_w` (how much the copies lowered the cut).
- **Direction, range.** Lower is better; [0, 1].
- **Thresholds.** pass ≤ 0.10; warn ≤ 0.30; fail > 0.30. Judgement. H is data-derived (section 1). The
  reason to keep non-hub species single: in Guimerà & Amaral 2005 the low-degree metabolites that link
  modules turned out better conserved across species than hubs linked mostly inside one module, so they
  carry information; drawing such a species twice hides the link between two boxes.
- **Does not apply.** `clone_in` not used (na). Not using copies is not penalised: under the weight, all
  hub edges together cost at most |H|/|S| of the cut.
- **Gaming.** None known. The hub test has a grey zone just below the cutoff (see section 6).

### G12 - Agreement with the curated reference (optional, not part of the graph score)
- **Curator check.** Nodes the curators put in one pathway share a box, and a box does not mix pathways.
- **Metric.** BCubed precision, recall and F1 over the labelled nodes, with "share at least one label" as
  the reference relation: precision_i = share of the labelled nodes in i's group that share a label with
  i; recall_i = share of the nodes sharing a label with i that are in i's group; P and R are the means;
  **`ref_bcubed_f1_best`** = the better F1 of L1 and leaf. Also `ref_cut_w`: the G1 cut of the reference
  itself (an unlabelled reaction takes the most frequent label of its labelled participants) - use it to
  see whether G1's thresholds are fair for that model.
- **Direction, range.** Higher is better; [0, 1].
- **Thresholds.** pass ≥ 0.60; warn ≥ 0.40; fail < 0.40. Judgement and **not calibrated**: I had no
  access to the reference. Prefer ranking strategies by this value to thresholding it. BCubed is the
  family that satisfies all formal constraints in Amigó et al. 2009; my variant differs from their
  extension for overlapping labels, which would penalise any partition for not duplicating a node that
  has two labels.
- **Does not apply.** No reference for the model (na).
- **Gaming.** One box gives recall 1, singletons give precision 1; F1 handles both. It is a held-out
  check, so it cannot be tuned to without breaking the task rules.

## 3. Behaviour by level, in one table

| metric | toward the leaf it ... | judged at | why |
|---|---|---|---|
| cut share (G1) | never decreases | leaf | the finest boxes are what is drawn |
| largest share (G2) | never increases | leaf | L1 may be a compartment, which is given |
| singleton share (G3) | usually increases | leaf | |
| connectedness (G4) | usually decreases | leaf | |
| modularity (G5) | has a peak | best level | resolution limit |
| weak-community share (G6) | usually decreases | leaf | |
| edge span (G7) | - | whole tree | integrates all levels |
| readable size (G8) | changes either way | leaf | |
| siblings (G9), single child (G10) | - | whole tree | |

`cut_w_mean` (mean of the per-level cuts) is returned as a diagnostic only: a dummy or lopsided top level
lowers it, which G7 does not allow.

## 4. Combining the rules into one score

Implemented as `judge(metrics)` and `score(metrics)` in `metrics_graph.py`.

1. Each rule gives pass = 1, warn = 0.5, fail = 0; na is left out.
2. Weighted mean over the applicable rules G1-G11 with weights G1 3; G2, G3, G4, G5, G7 2; G6, G8, G9,
   G10, G11 1. G12 is reported beside the score, not in it, so models with and without a reference stay
   on one scale.
3. Gates: if G2 or G3 fails, the score is capped at 0.40. A map with one dominant box or mostly
   single-node boxes is not a grouping, whatever its edges do.
4. One group at the leaf: G1, G5, G6 are fails, not na.
5. Over the six models: report the per-model scores, their mean and their minimum, and the full
   rule × model verdict table. Rank by the mean; break ties by the minimum, then by `edge_span_ratio`.
6. Do not average raw metric values across rules or models; their scales differ.

With these rules neither degenerate grouping comes out best on any of the six models (tested):

| synthetic grouping | score | G1 | G2 | G3 | G4 | G5 | G6 | G7 | G8 | G9 |
|---|---|---|---|---|---|---|---|---|---|---|
| `optimiser_nested` | 0.94-1.00 | 0.03-0.19 | 0.04-0.18 | 0.00-0.00 | 0.99-1.00 | 1.00-1.00 | 1.00-1.00 | 0.04-0.27 | 0.81-1.00 | 9-56 |
| `optimiser_clone_hubs` | 0.94-0.97 | 0.02-0.15 | 0.04-0.18 | 0.00-0.00 | 0.99-1.00 | 1.00-1.00 | 1.00-1.00 | 0.04-0.24 | 0.81-1.00 | 17-106 |
| `optimiser` | 0.94-0.97 | 0.03-0.19 | 0.04-0.18 | 0.00-0.00 | 0.99-1.00 | 1.00-1.00 | 1.00-1.00 | 0.05-0.28 | 0.81-1.00 | 17-106 |
| `dummy_root_optimiser` | 0.88-0.88 | 0.03-0.19 | 0.04-0.18 | 0.00-0.00 | 0.99-1.00 | 1.00-1.00 | 1.00-1.00 | 0.05-0.28 | 0.81-1.00 | 17-106 |
| `compartment_optimiser` | 0.74-0.91 | 0.12-0.33 | 0.03-0.18 | 0.00-0.18 | 0.83-1.00 | 0.83-1.00 | 0.88-1.00 | 0.14-0.42 | 0.79-1.00 | 14-106 |
| `compartment` | 0.25-0.66 | 0.00-0.17 | 0.49-1.00 | 0.00-0.17 | 0.49-1.00 | 0.25-0.60 | 1.00-1.00 | 0.51-1.00 | 0.00-0.51 | 1-10 |
| `per_reaction` | 0.31-0.40 | 0.25-0.66 | 0.00-0.07 | 0.21-0.56 | 1.00-1.00 | 0.44-0.83 | 0.28-0.86 | 0.27-0.68 | 0.36-0.80 | 47-856 |
| `random_sized` | 0.31-0.34 | 0.92-0.99 | 0.04-0.18 | 0.00-0.00 | 0.07-0.23 | -0.02-0.02 | 0.00-0.00 | 0.99-1.02 | 0.81-1.00 | 17-106 |
| `random_clone_all` | 0.29-0.32 | 0.84-0.97 | 0.04-0.18 | 0.00-0.00 | 0.07-0.23 | -0.02-0.02 | 0.00-0.00 | 0.90-0.99 | 0.81-1.00 | 17-106 |
| `species_apart_from_reactions` | 0.16-0.28 | 1.00-1.00 | 0.31-0.59 | 0.00-0.20 | 0.01-0.06 | 0.00-0.00 | 0.00-0.00 | 1.00-1.00 | 0.00-0.68 | 2-19 |
| `one_group` | 0.19-0.25 | n/a | 1.00-1.00 | 0.00-0.00 | 0.89-1.00 | n/a | n/a | 1.00-1.00 | 0.00-0.00 | 1-1 |
| `singletons` | 0.14-0.14 | 1.00-1.00 | 0.00-0.01 | 1.00-1.00 | n/a | 0.00-0.00 | 0.00-0.00 | 1.00-1.00 | 0.00-0.00 | 148-1817 |

Ranges are minimum-maximum over the six panel models; "n/a" = NaN. G1 `cut_w_leaf`, G2
`largest_share_leaf`, G3 `singleton_group_share_leaf`, G4 `connected_node_share_leaf`, G5 `Q_eff`, G6
`weak_node_share_leaf`, G7 `edge_span_ratio`, G8 `readable_node_share_leaf`, G9 `max_siblings`.

Rules that reward a degenerate grouping, and their counter-rules: G1 rewards few huge boxes (G2, G8);
G2 rewards crumbs (G3, G8, G9); G3 rewards one box (G2); G4 rewards one box and singletons (G2, G3);
G6 rewards huge boxes (G2, G8); G8 rewards random groups of ten (G1, G4-G7). G5 and G7 need none.

## 5. Synthetic groupings used for checking

Built by `metrics_graph.synthetic()`; none is a strategy output.

- `one_group`, `singletons`, `compartment` (one group per `loc_space`), `random_sized` (the optimiser's
  group sizes, members shuffled, seed 0).
- `optimiser` (the partition behind `Q_ref`), `optimiser_nested` (the same leaves under at most 9 merged
  top-level groups), `compartment_optimiser` (compartments at L1, modules cut by compartment at L2),
  `optimiser_clone_hubs` (hub species copied into every group that uses them).
- Adversarial: `dummy_root_optimiser` (one top-level group above the optimiser's partition),
  `random_clone_all` (random groups, every species copied wherever it is used), `per_reaction` (one group
  per reaction with the species that meet it first), `species_apart_from_reactions` (species by
  compartment, reactions by SBO term).

`optimiser*` are positive controls of the graph viewpoint only. They score `Q_eff` = 1 by construction and
say nothing about naming, localisation or biological sense.

## 6. What I am unsure about

1. **Strength of the weight.** 1/d_s lets leaf species dominate where the median degree is 1 (the three
   SKM models). G1 is then easier there (optimiser 0.03-0.05) than in metabolic models (0.12-0.19), with
   one threshold for both. A milder weight (for instance 1/√d_s) would sit between raw and normalised; I
   did not find a source that prefers one.
2. **Hubs in signalling models are regulators** (MYC2, SARD1, VPg), not currency. The test treats them
   like water and ATP. For G1 that is harmless (one vote each); for G11 it means copies of a regulator
   are accepted, which a curator may not want.
3. **Modifier edges** weigh the same as reactant and product edges.
4. **The viewpoint favours topology-based strategies.** Groupings by function or reaction type will warn
   or fail on G4-G6 by design. How much that should count is for the combined RULES.md, not for me.
5. **Calibration is on my own synthetic groupings**, and `Q_ref` comes from a greedy optimiser, not a
   proven optimum.
6. **G1 at the leaf** is harsh on an honest, very fine third level; G7 compensates only partly.
7. **G7 depends on model size**: the same box size is a smaller share of a larger map. Thresholds were set
   on models of 148-1817 nodes.
8. **G9**: whether 7 ± 2 applies to boxes on a map is not established; 9 may be too strict at L2.
9. **G8**: the 50-node bound comes from path-finding tasks on general graphs, not from pathway maps.
10. **Residual groups** (nodes left directly in a parent) are scored as groups. A strategy that parks
    hub species there on purpose loses a little on G3, G4, G6.
11. **G12** is uncalibrated, and I could test it on constructed labels only.
12. **Hub cutoff grey zone**: a species one reaction below the cutoff (for example NADH in e_coli_core, 12
    reactions, cutoff 12) is not a hub; its copies count as misuse in G11.

## 7. Sources (retrieved by web search; abstract or record level unless stated)

- Amigó E, Gonzalo J, Artiles J, Verdejo F (2009). A comparison of extrinsic clustering evaluation metrics
  based on formal constraints. Information Retrieval 12:461-486. doi:10.1007/s10791-008-9066-8
- Barber MJ (2007). Modularity and community detection in bipartite networks. Phys Rev E 76:066102.
  https://arxiv.org/abs/0707.1616
- Croes D, Couche F, Wodak SJ, van Helden J (2006). Inferring meaningful pathways in weighted metabolic
  networks. J Mol Biol 356(1):222-236. https://www.sciencedirect.com/science/article/abs/pii/S0022283605011599
- Dasgupta S (2016). A cost function for similarity-based hierarchical clustering. STOC '16:118-127.
  doi:10.1145/2897518.2897527
- Fortunato S, Barthélemy M (2007). Resolution limit in community detection. PNAS 104(1):36-41.
  doi:10.1073/pnas.0605965104
- Pathway identification by network pruning in the metabolic network of Escherichia coli. Bioinformatics.
  doi:10.1093/bioinformatics/btp575 (the retrieved record showed title and DOI, not the authors; used
  only for its description of Huss & Holme's procedure)
- Ghoniem M, Fekete J-D, Castagliola P (2004). A comparison of the readability of graphs using node-link
  and matrix-based representations. IEEE InfoVis 2004:17-24. doi:10.1109/INFVIS.2004.1
- Guimerà R, Amaral LAN (2005). Functional cartography of complex metabolic networks. Nature 433:895-900.
  doi:10.1038/nature03288
- Guimerà R, Sales-Pardo M, Amaral LAN (2004). Modularity from fluctuations in random graphs and complex
  networks. Phys Rev E 70:025101. doi:10.1103/PhysRevE.70.025101
- Huss M, Holme P (2007). Currency and commodity metabolites: their identification and relation to the
  modularity of metabolic networks. IET Syst Biol 1(5):280-285. doi:10.1049/iet-syb:20060077
- Miller GA (1956). The magical number seven, plus or minus two. Psychol Rev 63(2):81-97.
  doi:10.1037/h0043158
- Newman MEJ (2001). Scientific collaboration networks. II. Shortest paths, weighted networks, and
  centrality. Phys Rev E 64:016132. doi:10.1103/PhysRevE.64.016132
- Newman MEJ, Girvan M (2004). Finding and evaluating community structure in networks. Phys Rev E
  69:026113. The 0.3-0.7 reading was retrieved from pages citing this paper, not from its text.
  https://www.researchgate.net/publication/7659762_Finding_and_Evaluating_Community_Structure_in_Networks
- Radicchi F, Castellano C, Cecconi F, Loreto V, Parisi D (2004). Defining and identifying communities in
  networks. PNAS 101(9):2658-2663. doi:10.1073/pnas.0400054101
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing well-connected communities.
  Sci Rep 9:5233. doi:10.1038/s41598-019-41695-z
- Yang J, Leskovec J (2012). Defining and evaluating network communities based on ground-truth.
  https://arxiv.org/abs/1205.6233
- Yoghourdjian V, Yang Y, Dwyer T, Lawrence L, Wybrow M, Marriott K (2021). Scalability of network
  visualisation from a cognitive load perspective. IEEE TVCG 27(2):1677-1687. doi:10.1109/TVCG.2020.3030459
