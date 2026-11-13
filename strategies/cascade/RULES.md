# cascade - rules

Flow-anchored modules. A top-level group is either the cascade of one anchor (the steps that produce a
signal and the steps downstream of it), or a metabolic route (steps linked through shared main metabolites),
or one of two fallbacks (currency steps, unconnected steps). One program, no model-specific branch; which
rules fire depends only on SBO terms, participant roles and graph structure.

Run:

```
python strategies/cascade/strategy.py --features panel/features/<model> --sbml panel/models/<model>.xml --out strategies/cascade/out/<model>
```

Vocabulary used below

- *step* = reaction. *Small molecule* = species whose `sbo_path` contains SBO:0000247 (simple chemical).
- *Reactant, product, modifier* = `role` in `edges.tsv`.
- *Flow graph*: reactant -> reaction -> product and modifier -> reaction; for signalling reactions also
  reaction -> reactant (a controlled binding, state change or degradation changes the amount of its
  substrate, so the signal passes on to the substrate). *Supply graph*: the same without that last edge.
  *Mass graph*: reactant -> reaction -> product over metabolic and transport reactions only.
  Currency species (R2) are left out of all three.
- *Flow share*: leaving a node with k onward edges costs log k, so the cost of a path is -log of the share of
  the source's flow that arrives if flow splits evenly at every branching. Ties: fewer steps, then the family
  with more sensing reactions, then id.

## Rules, in the order applied

The `rule` column of `groups.tsv` holds the id of the rule that placed the node: R5-R8, R10-R14 for
reactions, R17a-d (or R14) for species.

**R1 - reaction class.**
- If `sbo_path` contains SBO:0000655 (transport reaction), SBO:0000185, SBO:0000627, SBO:0000628 or
  SBO:0000632 (translocation, exchange, demand, sink): class *transport*.
- Else if it contains SBO:0000589 or SBO:0000205 (genetic production, composite biochemical process): *genetic*.
- Else if it contains SBO:0000177 (non-covalent binding): *binding*.
- Else if every reactant and product is a small molecule: *metabolic*.
- Else: *state* (a macromolecule is changed: activation, inactivation, degradation, dissociation, ...).
- binding, state and genetic are the *signalling* classes.
- Reads: `nodes.tsv` `sbo_path`, `loc_kind`; `edges.tsv` `role`; species `sbo_path`.
- Fallback: a reaction without SBO term skips the first three tests, except that one reactant and one
  product in different spaces (`loc_kind` not `inside`) give *transport*. A species without SBO term counts
  as a small molecule.

**R2 - currency species.**
- R2a carrier pair: if two species stand on opposite sides of at least PAIR_MIN reactions, and these are more
  than half of the reactions of each of the two, then - taking candidate pairs in descending order of that
  count - the pair is a carrier pair if more than half of those reactions still have a reactant and a
  product once the pair and the carriers accepted before it are left out. (ATP/ADP, NAD+/NADH are found
  this way; a main substrate/product pair that occurs with three different cofactors is not, because
  removing it empties the reactions.)
- R2b hub: a small molecule that is reactant or product of more than CAP reactions is currency (it could not
  lie inside any route, see R10).
- R2c transported copy: in a transport reaction where all reactants and products but one are currency, the
  remaining one is currency too; repeated until nothing changes (extracellular H2O, phosphate, ...).
- Action: currency species never link steps, are never anchors, and are placed by R17.
- Reads: `edges.tsv` (reactant, product), R1 class.
- Fallback: a model in which no pair and no hub qualifies has no currency species (all three SKM models).

**R3 - anchors.**
- If a small molecule that is not currency is a modifier of a signalling reaction, or a reactant of a binding
  reaction, it is an anchor (a sensed molecule: ligand, second messenger).
- If a species has SBO:0000405 (perturbing agent) in its `sbo_path` and takes part in a reaction, it is an
  anchor (a stimulus).
- The reactions named in the condition are the anchor's *sensing reactions*.
- Reads: species `sbo_path`; `edges.tsv` `role`; R1 class.
- Fallback: a model without such species has no anchors; R4-R9 and R15 do nothing and every step goes
  through R10-R14 (brenda_merged, e_coli_core, R-HSA-77289).

**R4 - anchor families.** If the sensing reactions of two anchors share reactants or products (not counting
anchors) and the shared set is at least FAMILY_OVERLAP of the smaller of the two sets, the anchors form one
family (several ligands of one receptor). The family is led by the anchor with the most sensing reactions
(ties: id). Reads: `edges.tsv`. Fallback: an anchor sharing nothing is a family of one.

**R5 - production core.** If a metabolic or transport reaction lies on a mass-graph path into a
small-molecule anchor (it makes the anchor or a precursor) or out of it (turnover, conjugation, export),
it is placed in that anchor's family, stage *production*. Several anchors: the one fewest steps away, then
the family with more sensing reactions, then id. Reads: `edges.tsv`, R1, R2. Fallback: a stimulus anchor
and a ligand that no modelled reaction makes have no production core.

**R6 - catalyst supply.** If a reaction not yet placed makes, moves or activates a modifier of a production
core - it lies on a reactant -> product path into that modifier, a path that stops at anchors - it is placed
in that family, stage *production* (genetic production of the enzyme, its activation, its transport).
Several families: fewest steps, then family order. Reads: `edges.tsv`. Fallback: none needed.

**R7 - downstream.** If a reaction not yet placed can be reached from an anchor in the flow graph, it is
placed in the family that sends it the largest flow share. Stage: *genetic production* if the best path
contains a genetic reaction (the reaction itself included), otherwise *molecular interaction*.
This is the rule for nodes reachable from several anchors: exactly one family, the one with the largest
share. Reads: `edges.tsv`, R1. Fallback: none needed.

**R8 - side inputs.** If a reaction not yet placed reaches a placed reaction in the supply graph (it feeds
it), it is placed in the family of the nearest placed reaction by flow share (computed on the reversed
supply graph). Stage by what the step does: genetic reaction -> *genetic production*; metabolic reaction or
transport of small molecules -> *production*; otherwise *molecular interaction*.

**R9 - weak anchors.** If a family holds fewer than MIN_RX reactions after R5-R8, its anchors are removed
and R4-R8 are repeated, until every family holds at least MIN_RX reactions.

**R10 - routes.** Applied to the metabolic and transport reactions that no family took. Every non-currency
species that is reactant or product of two or more of them is a possible link. Links are taken in ascending
order of the number of reactions they join (weakest branch points first; ties: id). If the routes a link
would fuse hold together at most CAP reactions, they are fused; otherwise the species stays a branch point
between routes. So linear chains form first, minor branch points join them, and hubs separate the modules.
Reads: `edges.tsv`, R1, R2. Fallback: after R12 the rule runs once more on reactions of any class that are
still unplaced (for a signalling reaction the modifiers count as link species too); this is what groups
signalling content that has no anchor.

**R11 - small routes.** If a route holds fewer than MIN_RX reactions, it joins the route with which it shares
a non-currency species; the species joining the fewest reactions is tried first, and among the routes on it
the one with the most reactions on that species wins (ties: id). Smallest routes first, repeated until
nothing moves. Fallback: a small route that shares nothing is left to R12-R14.

**R12 - collateral steps.** If a reaction is still unplaced and is connected, in the flow graph taken without
direction, to a placed reaction, it joins the top-level group of the nearest one by flow share (outputs of a
regulator that sits in a cascade, a side reaction on a precursor). In a family its stage is set as in R8.
Applied after R11 and again after the second pass of R10.

**R13 - currency module.** If the model has currency species, they form one group together with the unplaced
reactions whose reactants and products are all currency species (ATP synthase, transhydrogenase, exchange
of H+ and water). The group is made only if it would hold at least MIN_LOC nodes. Fallback: otherwise there
is no currency module; such reactions fall to R14 and currency species are placed by R17c-d.

**R14 - unconnected steps.** Every reaction still unplaced, and every species that takes part in no reaction,
is placed in a localisation group of its own `loc_space` (one group per space, basis `localisation`).
A species homed with such a reaction by R17 is likewise placed in the group of its own `loc_space`.
Reads: `nodes.tsv` `loc_space`. This is the last fallback; no node is left out.

**R15 - stages (level 2 of a family).** The three stages set by R5-R8 and R12 are the level-2 groups of a
family. If a stage holds fewer than MIN_RX reactions it is merged into the neighbouring stage (production ->
molecular interaction; genetic production -> molecular interaction; molecular interaction -> genetic
production, else production). If only one stage remains the family has no level 2.

**R16 - route elements (level 2 of a route).** Within a route, the reactions lying on a directed cycle through
at least three non-currency species form the element *cycle* (direction reactant -> product as written; a
reaction with `reversible` = 1 may be passed either way). The other reactions fall into pieces connected
through species off the cycle; a piece of at least MIN_RX reactions is a *linear branch* element. If the
cycle or the rest of the route holds fewer than MIN_RX reactions, the route has no level 2. Reactions in
pieces below MIN_RX stay directly in the route. Reads: `edges.tsv`, `nodes.tsv` `reversible`.
Fallback: `reversible` empty counts as not reversible.

**R17 - species.** A species has one home, the group (levels 1-2) of
- R17a: for an anchor, the reaction of its own family that makes it (else that uses it);
- R17b: for a currency species, the currency module (when the model has one, R13);
- R17c: otherwise the reactions that make it (it is a product) - the group holding most of them, ties to the
  smaller group id;
- R17d: if no reaction makes it, the reactions that use it (reactant or modifier), chosen the same way.
- Fallback: a species in no reaction is placed by R14.

**R18 - clones.** If a species takes part in reactions of top-level groups other than its home, `clone_in`
lists one group per such top-level group: the deepest group there holding most of those reactions (ties:
id). Currency species and regulators shared between cascades therefore have one home and one copy per
other module that uses them. Species used only inside their home top-level group are never cloned.

**R19 - catch-all.** If a top-level group holds more than CATCH_ALL of all nodes and has at least two level-2
groups, it is replaced by its level-2 groups, which become top-level groups labelled
"<label of the group>: <own label>". Nodes that sat in it without a level-2 group keep the old group.
(A single-pathway model has one anchor; its three stages become the top level.) Applied once.

**R20 - localisation level.** Under every deepest flow group, if at least two `loc_space` values hold
MIN_LOC or more of its nodes each, each of those spaces becomes a child group (basis `localisation`,
`compartment_id` = the space) holding exactly the nodes placed directly in the group that have that
`loc_space`; boundary and between nodes go with their `loc_space` as nodes.tsv gives it. Nodes in spaces below MIN_LOC stay directly
in the flow group. Nesting of compartments is never encoded: all spaces are siblings and the renderer nests
them by compartments.tsv. Fallback: a group lying in one space (every group of a one-compartment model)
has no localisation level.

**R21 - names.** Every label comes from an identifier of the model or of SBO.
- Family: names of its anchors (at most three, then "+n"); `name_source` = member, `name_xref` = id of the
  leading anchor.
- Stage: the SBO term of the stage - production (SBO:0000393), molecular interaction (SBO:0000344),
  genetic production (SBO:0000589); `name_source` = sbo.
- Route and linear branch: names of its two most characteristic member species, "A / B"; cycle: "cycle
  through A". A species characterises a group by n^3/u^2, n = reactions of the group it is reactant or
  product of, u = all reactions it takes part in (currency species are skipped). `name_source` = member,
  `name_xref` = id of A.
- Currency module: names of the three currency species that are reactant or product of the most reactions
  in the model (same name listed once); `name_source` = member, `name_xref` = id of the first.
- Localisation group: the compartment name; `name_source` = compartment, `name_xref` = compartment id.
- Fallback: a route with no member species is unnamed (`name_source` = none). Did not occur on the panel.
- `basis`: flow for families, stages, routes, linear branches and the currency module; topology for a cycle;
  localisation for R14 and R20 groups.

## Parameters

| parameter | value | used in | how the value was chosen |
|---|---|---|---|
| CAP | ceil(sqrt(number of reactions)) | R2b, R10 | derived from the data: with about sqrt(N) routes of at most sqrt(N) reactions the page and the box are equally crowded. Values: brenda_merged 15, e_coli_core 10, R-HSA-77289 7, skm_pss 30, skm_ja 10, skm_sa 13 |
| MIN_RX | 3 reactions | R9, R11, R15, R16 | justified: three steps are the shortest chain with an inner step; a one- or two-step group shows no route |
| MIN_LOC | 5 nodes | R13, R20 | justified: the smallest complete two-step chain has five nodes (three species, two reactions) |
| PAIR_MIN | 3 reactions | R2a | justified: two reactions with the same opposite pair are typically one conversion written twice (forward/backward, isoenzymes); three is the smallest count that is not |
| majority tests | more than 1/2 | R2a | justified: "most of its reactions"; no free value |
| FAMILY_OVERLAP | 1/2 | R4 | justified: at least half of the smaller target set is shared |
| CATCH_ALL | 1/2 | R19 | justified: "a box holding most of the model" = more than half of the nodes |
| flow-share cost | log(k) per node with k onward edges | R7, R8, R12 | no free value: even split of flow |
| characteristic score | n^3/u^2 | R21 | justified, used for labels only: n reactions in the group times the squared share n/u of the species' reactions that lie in the group |
| anchors / species shown in a label | 3 / 2 | R21 | display only |

No value was tuned per model, and no value was fitted to reference labels (none were read). CAP with factor
1.5 and 2 was looked at once on the panel (top-level groups at factor 1 / 1.5 / 2: brenda_merged 22/18/17, e_coli_core 12/11/9, R-HSA-77289 11/10/9;
the three SKM models unchanged); factor 1 was kept because it is the only one with a derivation.

## How species are placed, and when clones are used

Reactions are placed first (R5-R14); a species follows its reactions (R17): it sits where it is made, and if
nothing makes it, where it is used. An anchor always sits in its own family, a currency species in the
currency module. Each species has exactly one home. A clone (R18) is declared only
when a species takes part in a reaction of another top-level group - one clone per such group, in the
deepest group that uses it most. On the panel 7 to 132 species per model carry clones; the species with the
most clones is a currency species (H2O in brenda_merged, 20 groups).

## Declared fallbacks (contract rule 7)

| missing | what happens | rule id on the node |
|---|---|---|
| reaction without SBO term | class from participants: metabolic if all small molecules, transport if one species crosses spaces, else state | as placed (R10-R14) |
| species without SBO term | treated as a small molecule | as placed |
| no anchor in the model | routes only | R10-R14 |
| reaction linked to nothing but currency | currency module if all its reactants and products are currency (R13), else unconnected | R13 / R14 |
| reaction linked to nothing | group of its `loc_space` | R14 |
| species in no reaction | group of its `loc_space` | R14 |
| one compartment only | no localisation level | - |
| `reversible` not set | not reversible | - |
| group without member species | unnamed | - |
