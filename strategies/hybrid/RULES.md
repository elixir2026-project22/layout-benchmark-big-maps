# RULES - strategy `hybrid`

A size-bounded hierarchical grouping, computed from the feature tables of one SBML model
(`nodes.tsv`, `edges.tsv`, `compartments.tsv`). One program, no model-specific branch.

- **Level 1 - functional module.** Annotation seeds, connectivity propagates: function terms written in the
  SBML notes start the modules where the network confirms them (R3); connectivity carries a term to the
  unannotated nodes (R4) and builds the modules alone where there is no annotation (R5).
- **Below level 1** - only where a group is above the leaf bound or spans several compartments:
  finer function terms (R9), reaction type from the SBO hierarchy (R10), localisation as given (R11),
  connectivity parts (R12, R13).
- **Split and merge rules** (R6, R7, R12, R13) keep the bounds below on every model.

## Bounds

| bound | value | enforced by | holds on the panel |
|---|---|---|---|
| number of level-1 groups | at most 12; aim `ceil(n / 40)` kept inside 4..12 | R7 | yes: 10, 5, 4, 12, 5, 9 |
| largest level-1 group | at most `max(40, floor(n / 3))` nodes | R6, and R7 never merges beyond it | yes: largest share 0.331 |
| every leaf group | at most 40 nodes | R9-R13 | yes: largest leaf 40 |
| every group | at least 3 nodes, **except localisation groups** | R7, R10, R12 | yes, with that exception (see NOTES.md) |

`n` = number of species plus reactions. Clones are not counted in sizes. The upper bounds are guaranteed by
construction; the aim of at least 4 level-1 groups and the 3-node minimum are not (NOTES.md, failure cases 5 and 12).

## Rules, in the order applied

Each rule: **if** (condition) **then** (action); **reads** (SBML attribute or feature column); **fallback**.

### R1 - function terms
- **If** the notes of a node hold the key `mapman`,
- **then** every entry `<code>_<name path>` in it (entries separated by `;`) is a function term of the node.
  A code `10.7.1` has the ancestors `10` and `10.7`. Term names are taken from the name path when it has as
  many dot-separated parts as the code (otherwise only the top-level name is taken). HTML entities are
  unescaped first; an entry starts only where `;` is followed by a numeric code and `_`, so a `;` inside a
  name does not split it.
- **Reads:** `nodes.tsv` column `notes` (JSON), key `mapman`.
- **Fallback:** a node without the key is unannotated; it is placed by R4 or R5. A model without the key has
  no terms; all its nodes are placed by R5 and R7.

### R2 - weighted graph and connectivity communities
- **If** always,
- **then** species and reactions form a bipartite graph, one edge per participation (reactant, product or
  modifier). An edge weighs `1 / (number of reactions of its species)`, so a species in 60 reactions ties
  them 30 times more weakly than a species in 2. Communities: Leiden, modularity, resolution 1, seed 0,
  nodes in id order.
- **Reads:** `edges.tsv` columns `reaction_id`, `species_id`.
- **Fallback:** a node without edges is a community of its own and is merged by R7.

### R3 - seeds: annotation confirmed by partners
- **If** an annotated node carries term `t` and at least one of its partners (nodes that share a reaction
  with it) carries `t` too,
- **then** the node becomes a seed of the term with the largest `observed - expected`, where `observed` is the
  number of partners carrying `t` and `expected` is the number of annotated partners times the share of all
  annotated nodes that carry `t`. Ties: the term carried by fewer nodes, then the smaller code. Terms are
  compared at the top of the term tree (`10`, `14`, ...), or one level below a term that R6 has split.
- **Reads:** `notes.mapman` (R1), `edges.tsv`.
- **Fallback:** a node none of whose terms is shared by a partner is treated as unannotated (R4, R5).

### R4 - propagation inside a community
- **If** a node has no term yet and has a neighbour with a term in the same community (R2),
- **then** it takes the term with the largest summed edge weight among those neighbours (ties: smaller code).
  Repeated in rounds, all nodes of a round at once, until nothing changes. A term never crosses a community
  border.
- **Reads:** `edges.tsv`, communities of R2.
- **Fallback:** a node in a community without any seed is placed by R5.

### R5 - connectivity module
- **If** a community holds no seed,
- **then** the community is a group of its own.
- **Reads:** communities of R2.
- **Fallback:** none needed; every node is in exactly one community.

### R6 - split a group above the share cap
- **If** a group holds more than `max(40, floor(n / 3))` nodes,
- **then** (a) if it is a term group and the term has child terms in the model, the comparison of R3 moves one
  level down for that term and R3-R5 are run again (a node annotated only at the split level loses that
  term); (b) otherwise the group is cut into connectivity parts within the cap, as in R12.
  Repeated until no group above the cap can be split.
- **Reads:** group sizes; the term tree of R1.
- **Fallback:** (b).

### R7 - merge down to the aim
- **If** there are more groups than the aim `K = ceil(n / 40)` kept inside 4..12, or a group has fewer than 3
  nodes,
- **then** the smallest group (ties: smaller key) joins another group, chosen among the groups it shares an
  edge with and whose union with it stays within the cap of R6: first the one with the most closely related
  term (number of leading code parts in common; applies to term groups only), then the one with the highest
  modularity gain (shared edge weight minus the product of the two groups' total edge weights divided by the
  total of the node set). If it shares no edge with any group that fits: the highest gain among all groups
  that fit. If no group fits the cap, it stays as it is. The larger group keeps its identity. Repeated.
- **Reads:** `edges.tsv`, group sizes, terms.
- **Fallback:** as stated (no connection: highest gain; nothing fits: unchanged).

### R8 - name a module after a function term
- **If** at least 3 members of a group carry term `t` (any depth of the tree) and they are more than half of
  the group's annotated members,
- **then** `t` may name the group. Candidates are ranked by F1 between the carriers of `t` in the model and
  the annotated members of the group (ties: the deeper term). Each term names one group only, the one it fits
  best. Label = the name path of the term as written in the notes (`Phytohormone action: jasmonic acid`);
  `name_source` = `notes`, `name_xref` = `mapman:<code>`, `basis` = `function`.
- **Reads:** `notes.mapman`.
- **Fallback:** a group no term qualifies for, or whose terms are taken, is named by R14.

### R9 - finer function terms inside a module (level 2)
- **If** a level-1 group has more than 40 nodes and is named after term `t` (R8),
- **then** R2-R8 are run inside it with the terms below `t`: aim `ceil(size / 40)` kept inside 2..12, cap
  `max(40, floor(size / 3))`. The result is kept if at least two children are named after a term; the
  children are level 2 (a child's label drops the part its parent shows) and each child is then handled by
  R11, R12 and R13 at level 3.
- **Reads:** as R1-R8.
- **Fallback:** fewer than two term-named children: the group goes to R10 and R11.

### R10 - reaction type (SBO)
- **If** a group above 40 nodes was not split by R9 and its reactions carry at least two different SBO terms,
- **then** its reactions are grouped by SBO term. A type group with fewer than 3 nodes is rolled up: it is
  united with the type that shares the deepest SBO ancestor and the union takes that ancestor's term (for
  example activation + inactivation -> conversion). The level is skipped when fewer than two types remain or
  one of two would stay below 3 nodes. Label = SBO term name, `name_source` = `sbo`, `basis` = `reaction_type`.
  **Order with R11:** both orders are built (type then localisation; localisation then type) and the order
  with fewer leaf groups below 3 nodes is used; equal: type first, as the track brief states. In the order
  localisation-first, a localisation group above 40 nodes is split by type only if every type group then
  fits 40 nodes; otherwise R12.
- **Reads:** `nodes.tsv` columns `sbo_path`, `sbo_path_names` of reactions.
- **Fallback:** a reaction without SBO term takes the type most frequent among the reactions that share a
  species with it (ties: smaller SBO id). No SBO term on any reaction, or one type only: the level is skipped.

### R11 - localisation
- **If** the nodes of a group lie in at least two places,
- **then** the group is split by the given place of each node: nodes with `loc_kind` `inside` or `between`
  go to the group of their `loc_space`; nodes with `loc_kind` `boundary` go to `<loc_space>@boundary` when
  there are at least 3 of them, otherwise to the group of their `loc_space` (CONTRACT.md rule 5 allows a
  space group to hold nodes whose `loc_space` it is). The level is skipped when it would give one child or
  when more than half of the children would be singletons. Nesting of compartments is never encoded.
  Label = compartment name (` (boundary)` appended), `name_source` = `compartment`, `name_xref` and
  `compartment_id` as CONTRACT.md requires, `basis` = `localisation`.
- **Reads:** `nodes.tsv` columns `loc_space`, `loc_kind`; `compartments.tsv` column `name`.
- **Fallback:** the level is skipped (one place only, or mostly singletons); the group stays a leaf, or goes
  to R12 if it has more than 40 nodes.

### R12 - connectivity parts keep the leaf bound
- **If** a group has more than 40 nodes, no rule above splits it and a level is free,
- **then** it is cut into connectivity parts of at most 40 nodes one level down: Leiden inside the group; a
  community above 40 is cut again (if Leiden returns one block the resolution is doubled, up to 64; last
  resort: the id-sorted list is cut into equal pieces); the parts are then merged by R7 down to the fewest
  the bound allows, `ceil(size / 40)`. Nodes without a neighbour inside the group are packed in id order.
  Parts are named by R14, `basis` = `topology`.
- **Reads:** `edges.tsv`.
- **Fallback:** the id-order cut.

### R13 - parts of one compartment (no level left)
- **If** a localisation group at level 3 has more than 40 nodes,
- **then** it is cut as in R12, and every part is written as a localisation group of the same compartment
  (same `compartment_id`), labelled `<compartment name> [<member name>]` with the member of R14.
- **Reads:** `edges.tsv`, `nodes.tsv` `loc_space`.
- **Fallback:** as R12.

### R14 - name a connectivity group after one member
- **If** a group has no term, SBO or compartment name (R5 groups without a qualifying term, R12 parts),
- **then** it is named after the species with the highest `k * (k / d)^2`, `k` = its reactions inside the
  group, `d` = all its reactions (ties: smaller id). Label = the species name (its id if it has no name),
  `name_source` = `member`, `name_xref` = the node id, `basis` = `topology`.
- **Reads:** `nodes.tsv` column `name`, `edges.tsv`.
- **Fallback:** a group without species is named after its first reaction in id order.

### R15 - clones of hub species
- **If** a species has more reactions than `Q3 + 1.5 * IQR` of the species degrees of its model,
- **then** it keeps its one home and gets a copy (`clone_in`) for each of its reactions outside its own groups:
  in the deepest group of that reaction that can hold the species - a group that is not a localisation
  group, or a localisation group of the species' own place. When the reaction sits in a localisation group
  of another place, the copy goes to the sibling group of the species' own place if there is one, else to
  the parent.
- **Reads:** `nodes.tsv` columns `degree` (recomputed from `edges.tsv`), `loc_space`, `loc_kind`.
- **Fallback:** no species beyond the fence: no clones. No other species is ever cloned.

## Parameters

Every value is the same for all models. None was tuned per model. "Fixed in advance" means the value was
set from the argument given and not varied; "chosen on the panel" means alternatives were looked at on all
six models together.

| parameter | value | used in | how the value was chosen |
|---|---|---|---|
| leaf bound `LEAF_MAX` | 40 nodes | bounds, R6, R7, R9-R13 | Fixed in advance. A textbook pathway drawn as one box has about this size (glycolysis: 10 reactions, 11 intermediates, cofactors: 30-40 nodes). Not varied. |
| minimum group size `GROUP_MIN` | 3 nodes | R7, R8, R10, R11, R12 | Fixed in advance. The smallest process is substrate, reaction, product. Also the support a term needs to name a group. |
| level-1 aim `K` | `ceil(n / LEAF_MAX)` inside 4..12 | R7 | Derived from the data (model size). Upper end 12: a page as a 3 x 4 grid of boxes; fixed in advance. Lower end 4: with the share cap of 1/3 three groups would have to be exactly equal, four leave room; fixed in advance. |
| share cap `SHARE_CAP` | 1/3 of the parent | R6, R7, R9 | Fixed in advance: no box holds as much as the two next ones could. Chosen on the panel in one respect: an extra cap of `2 / K` was tried and removed, because on skm_pss it split `Phytohormone action` at level 1 and R7 then re-merged the hormones in arbitrary pairs. |
| edge weight | `1 / degree of the species` | R2, R4, R7, R12 | Chosen on the panel. The alternative, removing hub species above the Tukey fence before clustering, broke the SKM networks into many fragments (187 communities on skm_pss against 110). The weight has no free value. |
| Leiden resolution, seed | 1.0, 0 | R2, R9, R12 | Defaults (plain modularity); the seed only makes runs repeatable. Not varied. |
| seed score | observed - expected partners | R3 | No free value. Chosen on the panel over a plain count of partners, which let the most frequent bin win every tie. |
| naming majority | more than 1/2 of annotated members, at least `GROUP_MIN` carriers | R8 | Fixed in advance (majority rule). |
| "mostly singletons" | more than half of the children | R11 | The wording of the track brief. |
| type/localisation order | fewer leaf groups below `GROUP_MIN`; ties: type first | R10 | No free value. |
| member-name score | `k * (k / d)^2` | R14 | Chosen on the panel. With exponent 1, five of the nine level-1 groups of e_coli_core and R-HSA-77289 were named after ATP, NAD, CoA-SH, FAD or the trifunctional protein (5 to 8 of their 12 to 13 reactions inside). |
| hub fence `HUB_FENCE` | `Q3 + 1.5 * IQR` of species degrees | R15 | Tukey's standard outlier fence, computed from each model's own degree distribution; the factor 1.5 was not varied. |
| function note keys | `mapman` | R1 | The only notes key in the panel that holds a hierarchical vocabulary together with its names. Other keys can be added to `FUNCTION_NOTE_KEYS` if they use the form `<dotted code>_<dotted name path>`. |

## How species are placed

A species has exactly one home, decided level by level:

1. **Level 1.** By its own annotation if a partner confirms it (R3). Otherwise by its reactions: it takes the
   term that most of its edge weight points to inside its community (R4) - since all edges of one species
   weigh the same, this is the majority of its labelled reactions there. Otherwise by its community (R5).
   Then it moves with its group through R6 and R7. The rule that decided this is written to the column `rule`
   of `groups.tsv` (R3, R4, R5, R6 or R7) for species and reactions alike.
2. **Function level below (R9).** As level 1, inside the module.
3. **Reaction-type level (R10).** A species goes to the type group that holds most of its reactions within
   the parent group; ties: the group where it is a product, then the smaller SBO id. A species with no
   reaction inside the parent goes to the largest type group.
4. **Localisation level (R11, R13).** By its own `loc_space` and `loc_kind`, never by its reactions.
5. **Connectivity parts (R12, R13).** By the Leiden community it falls in.

Fallbacks (CONTRACT.md rule 7): no function term -> R4/R5; no SBO term on a reaction -> type of the
neighbouring reactions, or the level is skipped; one compartment only -> the localisation level is skipped;
no edges -> a community of its own, merged by R7.

## When clones are used

Only for hub species (R15): degree above `Q3 + 1.5 * IQR` of the model's species degrees. The fence is
9.5 (brenda_merged), 12 (e_coli_core), 3.5 (R-HSA-77289), 6 (skm_pss), 3.5 (skm_ja), 3.5 (skm_sa) reactions.
The home of a hub is decided like that of any other species. `clone_in` lists the further groups in which a
copy is drawn; a copy never sits in a localisation group of a place other than the species' own.

## What the output columns record

- `rule` (groups.tsv): the rule that decided the node's level-1 group - R3, R4, R5, R6 or R7.
- `basis` (groups_meta.tsv): `function` = named after a function term that the majority of the group's
  annotated members carry (R8, R9); `topology` = a connectivity group named after a member (R5, R12, R14);
  `reaction_type` = R10; `localisation` = R11, R13.
- Group ids: `M01`, `M02`, ... are the level-1 groups, largest first. Children append `.f<i>` (finer
  function), `.sbo<id>` (reaction type), `.<compartment_id>` or `.<compartment_id>_b` (localisation,
  boundary), `.p<i>` (connectivity part).
