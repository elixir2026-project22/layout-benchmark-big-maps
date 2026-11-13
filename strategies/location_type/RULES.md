# RULES - strategy `location_type`

Cell localisation on top, reaction type inside each compartment. This is the inverse of the
type-then-location hierarchy and serves as its control.

- Level 1 (L1): where a node is drawn in the cell - one group per space, one `@boundary` group per space
  whose boundary carries transport or membrane steps. The drawing nests these groups as one cell from
  `compartments.tsv`; the levels of this strategy never encode which compartment lies inside which.
- Level 2 (L2): the type of the reactions inside each L1 group - the SBO process term, or the EC class
  where the SBO term does not divide the group.
- Level 3 (L3): only inside an L2 box that is too large to read - the next finer type.

The rules are applied in the order R1 to R12. Rules R8 to R10 are applied every time a type level is
made (after R4 to R7). The `rule` column of `groups.tsv` holds the id of the rule that decided the
node's innermost group (see "The rule column" below).

Inputs are the three feature tables only. The SBML file is not read.

## Rules

### R1 - localisation group (level 1)
- **If** a node has `loc_kind` = `inside` or `between`,
- **then** its L1 group is the group of its `loc_space`: id `loc_<compartment_id>`, basis `localisation`,
  `compartment_id` = the space. Label = the compartment's `name`; `name_source` = `compartment`,
  `name_xref` = the compartment id. A `between` node (it connects sibling spaces) belongs to their common
  parent, which is what `loc_space` already holds.
- **Reads:** `nodes.tsv` `loc_space`, `loc_kind`; `compartments.tsv` `compartment_id`, `name`.
- **Fallback:** a compartment without a name is labelled with its id. Every node has a `loc_space` (the
  feature tables give the root space to a node without compartment), so no node is left out.

### R2 - boundary group (level 1)
- **If** a node has `loc_kind` = `boundary` (a reaction whose participants lie on both sides of the
  boundary of a space, or a node in a membrane compartment),
- **then** its L1 group is `loc_<loc_space>@boundary`, basis `localisation`, `compartment_id` =
  `<loc_space>@boundary`. Label: the name of the membrane compartment that encloses the space, if the
  model has one (then `name_xref` = that membrane's compartment id); otherwise "`<space name>` boundary"
  with `name_xref` = the space's compartment id. `name_source` = `compartment` in both cases.
- **Reads:** `nodes.tsv` `loc_space`, `loc_kind`; `compartments.tsv` `kind`, `encloses_id`, `name`.
- **Fallback:** none needed; a boundary group exists only where a node has `loc_kind` = `boundary`.
  A boundary group is kept however small it is: localisation is given, not chosen.

### R3 - one compartment only (the localisation level is skipped)
- **If** R1 and R2 give exactly one group for the whole model (a one-compartment model with no boundary node),
- **then** localisation says nothing and is not used as a level. The first type level (R4, else R5)
  becomes L1, and R11 may add L2 and L3. The nodes that no type box takes (R8, R10) are placed in the
  compartment's own group (as in R1), which then stands at level 1 beside the type boxes.
- **Reads:** the set of L1 groups from R1 and R2.
- **Fallback:** if no type key divides the model either, the single compartment group holds every node
  (L1 only, rule id R1).
- **Variant:** `--keep-single-compartment` switches R3 off. The single compartment is then L1 and holds
  100% of the nodes; type is L2 and L3 as in every other model.

### R4 - type by SBO process term (first type level)
- **If** the reactions of an L1 group (under R3: of the model) fall into at least two SBO boxes,
- **then** each box becomes a group: id `<parent>/sbo_<number>`, basis `reaction_type`, label = the SBO
  term's name, `name_source` = `sbo`, `name_xref` = the SBO id.
  - A reaction's term is its own `sbo`.
  - Roll-up: a term carried by fewer than MIN_RX reactions of the group is replaced by its SBO parent
    (deepest terms first, repeated), so that rare terms join their parent's box. The roll-up never
    reaches "process" (SBO:0000375) or a direct child of it (such as "biochemical or transport
    reaction"): these say too little to name a box.
  - A box is an SBO term that, after roll-up, holds at least MIN_RX reactions.
- **Reads:** `nodes.tsv` `sbo`, `sbo_name`, `sbo_path`, `sbo_path_names` (reactions).
- **Fallback:** a reaction without SBO term, or whose term stays rare, is handled by R8. If fewer than two
  boxes result (no SBO term at all, or one term for the whole group), R4 does not apply to this group
  and R5 is tried.

### R5 - type by EC class (fallback of R4; also the next key of the ladder)
- **If** R4 does not divide the group, and its reactions fall into at least two EC classes with at least
  MIN_RX reactions each,
- **then** each class becomes a group: id `<parent>/ec_<class>`, basis `reaction_type`, label
  "EC `<class>` `<class name>`" (e.g. "EC 1 oxidoreductases"), `name_source` = `ec`, `name_xref` =
  `ec-code:<class>`. A reaction with several EC numbers goes to the class most of them fall in (ties: the
  lowest class number).
- **Reads:** `nodes.tsv` `ec` (reactions).
- **Fallback:** a reaction without EC number, or of a rare class, is handled by R8. If fewer than two
  boxes result, the group is not divided, unless it holds more than MAX_BOX nodes: then R7 is tried.

### R6 - type by EC sub-class and sub-sub-class (ladder keys, used through R11 only)
- **If** R11 asks for one more level inside an EC box,
- **then** its reactions are divided by the next field of the EC number: sub-class (`1.1`, `2.7`) inside a
  class box, sub-sub-class (`1.1.1`) inside a sub-class box. Id `<parent>/ec_<prefix>`, label
  "EC `<prefix>`", `name_source` = `ec`, `name_xref` = `ec-code:<prefix>`. Only EC numbers under the
  parent box's prefix count; most frequent prefix wins, ties go to the numerically lowest.
- **Reads:** `nodes.tsv` `ec` (reactions).
- **Fallback:** a reaction whose EC number is not that deep (e.g. `1.2.-.-`), or with a rare prefix: R8.
  Fewer than two boxes: the next key (R7) is tried. The fourth EC field is a serial number, not a class,
  and is never used.

### R7 - type by the SBO role of the modifiers (last ladder key)
- **If** R11 asks for one more level and the earlier keys did not divide the box (or, at the first type
  level, an L1 group of more than MAX_BOX nodes that R4 and R5 did not divide),
- **then** its reactions are divided by the SBO term of their modifier references ("stimulator",
  "inhibitor", "catalyst"). Id `<parent>/role_<number>`, basis `reaction_type`, label "with `<SBO name>`",
  `name_source` = `sbo`, `name_xref` = the SBO id. A reaction with modifiers of several roles takes the
  most frequent one (ties: lowest SBO id). The generic term "modifier" (SBO:0000019) is ignored: it is the
  root of the role terms and says too little to name a box.
- **Reads:** `edges.tsv` `role`, `sbo`, `sbo_name`.
- **Fallback:** a reaction without modifier, or whose modifiers carry no specific role term: R8. Fewer
  than two boxes: the box is not divided.

### R8 - reactions that get no box
- **If**, at a level being made, a reaction (a) has no value for the key in use, or (b) has a value that
  fewer than MIN_RX reactions of the group share, or (c) belongs to a box that would hold fewer than
  MIN_BOX nodes once species are placed (R9),
- **then** it gets no box at this level and stays directly in the enclosing group, drawn beside the
  boxes. In case (c) the box is withdrawn and the key is applied again without that value; if fewer than
  two boxes remain, the key does not divide the group and the next key is tried.
- **Reads:** the same columns as the key in use; group sizes after R9.
- **Fallback:** this rule is the fallback of R4 to R7. Under R3 the enclosing group is the compartment's
  own group.

### R9 - home of a species
- **If** a species takes part in at least one reaction of its own L1 group (under R3: in any reaction),
- **then**, at each type level, its home is the box holding most of its reactions, counted over the
  reactions inside the group being divided (one count per row of `edges.tsv`). "Directly in the
  enclosing group" (reactions placed by R8) competes as one more candidate. Ties are broken by: most
  rows as product, then most rows as reactant, then a box before "directly in the enclosing group",
  then the lowest group id. A species never leaves its L1 localisation group: R1 and R2 place it, R9
  only chooses the type box inside.
- **Reads:** `edges.tsv` `reaction_id`, `species_id`, `role`; `nodes.tsv` `loc_space`, `loc_kind`.
- **Fallback:** R10.

### R10 - species with no reaction in its own level-1 group
- **If** a species takes part in no reaction of its own L1 group (it is only transported across a
  boundary, only acts in another compartment, or has no reaction at all),
- **then** it stays directly in its L1 group, outside the type boxes. Under R3 this is the compartment's
  own group.
- **Reads:** `edges.tsv`; `nodes.tsv` `loc_space`, `loc_kind`.
- **Fallback:** none needed.

### R11 - one more level inside a box that is too large
- **If** a type box holds more than MAX_BOX nodes (reactions plus the species at home in it) and a level
  is free below it (at most three levels in all),
- **then** its reactions are divided by the next key of the type ladder that divides them, and its
  species are placed again by R9 among the new boxes. The ladder is: SBO term (R4) > EC class (R5) > EC
  sub-class > EC sub-sub-class (R6) > modifier role (R7). A key already used above the box is not used
  again; an EC sub-level is used only inside an EC box. R8 applies to the reactions the key does not place.
- **Reads:** group sizes; the columns of the key in use.
- **Fallback:** if no remaining key gives at least two boxes, the box stays as it is, larger than MAX_BOX
  (see NOTES.md, failure case 2).

### R12 - clones of a species
- **If** a species takes part in reactions placed in groups other than its home, inside its own L1
  localisation group,
- **then** a copy is drawn (`clone_in`) in such a group when (a) the species is only a modifier there, or
  (b) the species takes part in reactions of HUB_BOXES or more different groups (its home included).
  No copy is made in a group that encloses the home, and none outside the species' own L1 localisation
  group: a clone never crosses a compartment boundary.
- **Reads:** `edges.tsv` `role`; the groups from R1 to R11.
- **Fallback:** in every other case the species is drawn once, in its home, and its edges cross the box borders.

## The rule column of `groups.tsv`

| value | node | meaning |
|---|---|---|
| R1, R2 | reaction or species | its L1 localisation group was not divided; the node sits directly in it |
| R4, R5, R6, R7 | reaction | the key that defines its innermost type box |
| R8 | reaction | it got no box at the first type level and sits directly in the L1 group (under R3: in the compartment's own group). A reaction left over at a deeper level keeps the rule id of the box it stays in |
| R9 | species | placed by majority of its reactions (in a type box, or directly in the enclosing group if most of its reactions are there) |
| R10 | species | no reaction in its own L1 group; sits directly in it |

R3, R11 and R12 do not appear in the column: R3 and R11 decide which levels exist, R12 fills `clone_in`.

## Parameters

No parameter is set per model; the same values are used for all six.

| parameter | value | how it was chosen |
|---|---|---|
| MAX_BOX | 50 nodes | Justified, not derived: the most nodes one box can hold and still be read as one diagram (about 20 reactions with their participants). It only decides where a further level is tried. Changing it to 30 or 80 does not change any Stress Knowledge Map group at L2 and moves the panel total of L3 groups between 10 and 2 (table below). |
| MIN_BOX | 5 nodes | Justified: one tenth of MAX_BOX, so that the boxes of one level differ in size by at most one order of magnitude. Five nodes is also one reaction with two substrates and two products; a smaller box frames less than one such reaction. |
| MIN_RX | 3 reactions | Justified: a type box stands for a class of reactions. One reaction is an instance and two are a pair; three is the smallest number that shows a repeated type. Also the threshold of the SBO roll-up in R4. |
| HUB_BOXES | 3 groups | Justified by geometry: a species shared by two boxes is their link and can be drawn once, on their common border. One shared by three or more boxes cannot border them all, so it is copied. |
| SBO roll-up ceiling | "process" (SBO:0000375) and its direct children | A property of the ontology, not of a model: these are the terms shared by nearly all reaction types, so a box named after one would not say what it holds. |
| ignored role term | "modifier" (SBO:0000019) | The root of the modifier role terms; same reason as above. |
| EC depth | 3 fields | The fourth EC field numbers the individual enzyme; it is not a class. |
| number of levels | 3 | Fixed by the contract. |
| minimum number of boxes for a key to divide a group | 2 | A level with one box repeats its parent. |

Effect of each parameter on the panel totals (sum over the six models, one parameter changed at a time;
produced by running `strategy.py` with `--max-box`, `--min-box`, `--min-rx`, `--hub-boxes`):

| setting (others at default) | L1 groups | L2 groups | L3 groups | type boxes < 5 nodes | innermost boxes > 50 nodes | nodes outside sub-boxes (reactions) | clone copies |
|---|---|---|---|---|---|---|---|
| default | 50 | 66 | 6 | 0 | 9 | 138 (19) | 414 |
| MAX_BOX=30 | 50 | 71 | 10 | 0 | 9 | 170 (36) | 457 |
| MAX_BOX=80 | 50 | 62 | 2 | 0 | 11 | 135 (18) | 386 |
| MIN_RX=2 | 50 | 67 | 6 | 0 | 9 | 156 (26) | 414 |
| MIN_RX=5 | 50 | 50 | 6 | 0 | 9 | 226 (57) | 375 |
| MIN_BOX=3 | 50 | 73 | 6 | 5 | 9 | 139 (20) | 414 |
| MIN_BOX=8 | 50 | 53 | 6 | 0 | 9 | 208 (60) | 370 |
| HUB_BOXES=2 | 50 | 66 | 6 | 0 | 9 | 138 (19) | 640 |
| HUB_BOXES=4 | 50 | 66 | 6 | 0 | 9 | 138 (19) | 300 |

The defaults are the row "default". The delivered outputs were produced with the defaults only.

## How species are placed, and when they are cloned

1. A species is first placed by localisation: its L1 group is the group of its `loc_space` (R1), or the
   `@boundary` group of that space if it lies in a membrane compartment (R2). This never changes.
2. Inside that L1 group, the species goes to the type box that holds most of its reactions (R9). Only
   reactions of the same L1 group count, so a transport reaction (which sits in a boundary group) never
   pulls a species out of its compartment.
3. A species all of whose reactions lie outside its own L1 group stays directly in that group (R10).
4. If the home box is divided once more (R11), the species is placed again among the new boxes by the same count.
5. A species has exactly one home. Copies are listed in `clone_in` (R12) in two cases: the species is only
   a modifier in the other group (an enzyme or regulator is drawn next to the reaction it acts on), or
   the species takes part in three or more groups (a hub such as ATP or water in a metabolic model).
   A species shared by exactly two groups as reactant or product is not cloned: its edge between the
   two boxes shows how they connect.
6. Clones stay inside the species' own L1 localisation group. A cytosolic regulator of a nuclear reaction
   is drawn once, in the cytosol, and its edge crosses the nuclear boundary.

## What the strategy does with a one-compartment model

Localisation gives one group, which would hold the whole model and say nothing. R3 therefore skips the
level: the type hierarchy starts at level 1. For `brenda_merged` (one compartment, every reaction
"biochemical reaction") R4 does not divide the model either, so level 1 is the EC class (R5) and level
2 the EC sub-class (R6, through R11). The eight nodes no EC box takes sit in the compartment's own
group. The variant `--keep-single-compartment` keeps the compartment as level 1 instead.

## Where names come from

| group | label | `name_source` | `name_xref` |
|---|---|---|---|
| space (R1) | compartment name from the SBML | `compartment` | compartment id |
| boundary (R2) | name of the enclosing membrane compartment, else "`<space name>` boundary" | `compartment` | membrane compartment id, else the space's id |
| SBO type (R4) | SBO term name | `sbo` | SBO id |
| EC class (R5) | "EC `<n>` `<class name>`" | `ec` | `ec-code:<n>` |
| EC sub-class (R6) | "EC `<n.m>`" | `ec` | `ec-code:<n.m>` |
| modifier role (R7) | "with `<SBO role name>`" | `sbo` | SBO id |

No group is unnamed and no name is guessed from words. The names of the seven top-level EC classes
(oxidoreductases, transferases, hydrolases, lyases, isomerases, ligases, translocases) are the one
vocabulary written into `strategy.py`; SBO names come from `nodes.tsv` and `edges.tsv`.
