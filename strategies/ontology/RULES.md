# Strategy `ontology` - rules

Boxes come from the controlled identifiers the SBML carries. A reaction is the unit that is grouped; a
species follows its own term or its reactions. Rules are listed in the order `strategy.py` applies them.
The `rule` column of `groups.tsv` holds the rule that gave the node its level-1 box
(reactions: R1, R2, R3, R8, R9 or R10; species: R11, R12 or R10).

Words used below
- **term path**: a term with its ancestors, top of the tree first (`mapman:10 > mapman:10.7 > mapman:10.7.1`).
- **box**: a group. **Candidate box** of a path: the first term on the path, from the top, that is still allowed
  to be a box (R5-R7 decide which terms are allowed).
- **process ontologies**: bins, pathway cross-references, GO. **EC** describes the chemistry of one reaction.

## Rules

| id | condition (if) | action (then) | reads | fallback |
|---|---|---|---|---|
| **R1** own term | The node carries an identifier of one of four ontologies: (a) **bins** - a notes key whose whole value is a `;`-list of `<dotted code>_<dotted names>` (MapMan: `10.7.1_Phytohormone action.jasmonic acid.biosynthesis`); (b) **EC** - an `ec-code` cross-reference, the notes key `ec_numbers`, or a bin whose name segments are `EC_<n>[-<n>...] <text>` (MapMan writes the EC tree as the bin "Enzyme classification"; such a path is read as EC, and `<text>` becomes the name of that EC term); (c) **pathway cross-reference** - database `kegg.pathway`, `wikipathways`, `smpdb`, `panther.pathway`, `unipathway`, `pw`, `pathwaycommons`; a `kegg` id of the form `<2-4 letters><5 digits>` (`ath00906`, `map00908`, `ko04075` are all read as `kegg.pathway:map<5 digits>`); a `biocyc`/`metacyc` id or a `...?object=<id>` URL whose id contains the token `PWY`; a `reactome` id with qualifier `isPartOf`; (d) **GO** - a `go` cross-reference. | The node gets one term path per identifier. Bins and EC give trees (EC: `1 > 1.3 > 1.3.99 > 1.3.99.3`; a trailing `-` is cut). Pathway and GO terms are flat (one level): the bundle has no GO graph and no pathway hierarchy. A reaction also gets its SBO path (used only by R9). | `nodes.tsv`: `notes`, `ec`, `xrefs`, `go`, `sbo_path`, `sbo_path_names` | No identifier of an ontology: reaction -> R2 for that ontology; species -> no own term (R12). |
| **R2** enzyme inheritance | A reaction has no own term of an ontology, and at least one of its modifier species has. | The reaction takes the terms of those modifiers for that ontology (each modifier is one *source*). | `edges.tsv`: `role` = `modifier`; R1 terms of the species | No modifier with a term: R3 (process ontologies), none (EC). |
| **R3** substrate/product inheritance | As R2 but no modifier has a term, the ontology is a process ontology, and a reactant or product has a term. | The reaction takes the terms of those reactants and products. Not applied to EC: the EC number of a substrate protein says nothing about this reaction. | `edges.tsv`: `role` = `reactant`, `product` | No term reaches the reaction: R8. |
| **R4** choice between ontologies and between terms | A reaction has candidate boxes. | **a.** Ontology order: bins, then EC, then pathway cross-reference, then GO; the first ontology that offers a candidate box is used, the others are ignored for this reaction. **b.** Inside that ontology: the box named by the most sources (R1: the reaction itself; R2: its modifiers; R3: its reactants and products). **c.** Still tied: the tied box in which most of the reaction's participant species already live (species homes by R11/R12 from the reactions decided so far); repeated in rounds, each round using the reactions decided before it. **d.** Still tied after the rounds stop: the reactions with the fewest tied boxes take the tied box that already holds more reactions (then the one more reactions name, then the smaller id in natural order), and c resumes. | sources from R1-R3; `edges.tsv` | - |
| **R5** depth: when a term is a box | A term is a candidate box. | It stays a box only if at least `M_MIN` = 3 reactions end up in it. Otherwise it is dropped and R4 is run again without it: its reactions take another candidate, or the next ontology, or go on to R8. Level 1 starts at the top terms of each tree (bin `10`, EC class `1`). **Single-child chain**: if every reaction that reached a box through a term path continues to the same child term, the box is that child instead (repeated) - so a box is named by the deepest term that still holds all of it (`EC 4.2.1.17`, not `EC 4`, when every lyase of the model is 4.2.1.17). Drop order when several boxes are too small: ontology order of R4a, deepest term first, boxes that fewer than `M_MIN` reactions name at all first, then the smallest; one box per pass. | reaction counts per box | Dropped box: its reactions fall through R4 -> R8 -> R9 -> R10. |
| **R6** page rule | More than `PAGE` = 12 level-1 boxes. | The box with the fewest reactions is dropped (as in R5), until 12 remain. | reaction counts per box | as R5 |
| **R7** catch-all rule | A level-1 box of a tree ontology (bins, EC) holds more than `S_MAX` = 0.5 of all nodes (species and reactions, after R8-R12) and its term has child terms. | The term is *opened*: its child terms become level-1 candidates in its place. Reactions with no child term, or whose child was dropped, stay in a box labelled `<term> (other)` (itself subject to R5). R5-R7 are repeated until nothing changes. This is how depth follows model size: a one-pathway model is opened down to the sub-pathway terms. | node counts per box | Flat boxes (GO, pathway, SBO) and `(other)` boxes cannot be opened; they stay large. |
| **R8** propagation to reactions without a term | A reaction has no box after R1-R7. | It takes the box in which most of its participant species live, counting only boxes of process ontologies (ties: more modifiers there, then id). Repeated in rounds until no reaction changes. EC boxes are not propagated: a transport or exchange reaction is not an oxidoreductase because its substrate is used by one. | `edges.tsv`; species homes (R11/R12) | No participant lives in a process box: R9. |
| **R9** fallback: reaction type | A reaction has no box after R8 and has an SBO term. | Box = the deepest term on its SBO path that gathers at least `M_MIN` such reactions (its own term first; a term with too few is dropped by R5 and its reactions move to the parent term). | `nodes.tsv`: `sbo_path`, `sbo_path_names` | No SBO term, or no term on the path gathers 3: R10. |
| **R10** last resort | A reaction has no box after R9. | It takes the box in which most of its participants live, whatever its ontology (rounds). If none of its participants has a box, it goes to one leftover box (id `unassigned`), which is named by the deepest SBO term all its reactions share, or unnamed (`name_source` = `none`, basis `mixed`) if they have none. A species none of whose reactions has a box, or with no reaction, goes there too. | `edges.tsv`; `sbo_path` | - (every node is placed) |
| **R11** species with an own term | A species has an own term (R1) whose candidate box exists (holds reactions). | It lives in that box. Several such boxes: the one holding the greatest weight of its reactions (weights as in R12), then the box with more reactions, then id. Only the first ontology of the R4a order that gives an existing box is used. | R1 terms; `edges.tsv` | Own box does not exist: R12. |
| **R12** species by its reactions | A species has no own term, or none whose box exists. | It lives in the box with the greatest weight of its reactions; a reaction weighs 1 / (number of its species), so the reactions a species shares with few others count most. If a reaction names its main metabolites (notes key `main_metabolites`, tokens matched to the one participant id ending in the token) and the species is a main metabolite of at least one reaction, only those reactions vote. | `edges.tsv`; `nodes.tsv`: `notes` | No reaction with a box: R10. |
| **R13** second level | A level-1 box is a term of a tree ontology (not an `(other)` box). | R4, R5 (with `M_MIN`, single-child chain) and R6 (at most `PAGE` sub-boxes) are repeated one tree level below the box's term. Reactions that came by R8/R10 take the sub-box of most of their participants; species follow R11/R12 inside the box. A single sub-box is kept only if it holds more than `S_MAX` of the box's reactions. Nodes without a sub-box stay directly in the level-1 box. R7 is not applied at level 2. | as R4-R6 | No sub-box: the box has no function level 2 (R14 may still divide it). |
| **R14** localisation as the innermost level | A function box with no function sub-box (a level-2 box, or a level-1 box without level 2) has nodes in at least two places. A place is `loc_space`, with `loc_kind` = `boundary` counted as a separate place `<loc_space>@boundary`. | Each place with at least `M_MIN` nodes becomes a localisation box one level down (`compartment_id` = `loc_space` or `loc_space@boundary`). A boundary part with fewer than `M_MIN` nodes joins the box of its space when that space is also present. Boxes are made only if at least two places qualify. `between` nodes belong to their `loc_space`. No compartment nesting is encoded: every localisation box sits directly in its function box. | `nodes.tsv`: `loc_space`, `loc_kind`; `compartments.tsv`: `name` | Fewer than two qualifying places: no localisation level; nodes of small places stay directly in the function box. |
| **R15** clones | A species has reactions in at least `CLONE_MIN_BOXES` = 3 different level-1 boxes. | `clone_in` lists, for each of its reactions outside its home level-1 box, the function box of that reaction - or that box's localisation part for the species' own place, if there is one. | `edges.tsv`; R14 boxes | Fewer than 3 level-1 boxes: no clone; the species is drawn once and its edges cross boxes. |

## Depth, per ontology

| ontology | tree in the bundle | level 1 | level 2 | how depth is chosen |
|---|---|---|---|---|
| bins (notes; MapMan) | yes, with names | top bin, or deeper after R5 chain / R7 | next bin level (R13) | R5 (>= 3 reactions), R7 (open a catch-all), R6 (<= 12) |
| EC | yes; names only for the 7 classes, plus whatever the notes give | EC class, or deeper after R5 chain / R7 | next EC level (R13) | same |
| pathway cross-reference | no | the pathway id | - | none possible: exact identifier |
| GO | no (no GO graph in `common/`) | the GO id | - | none possible: exact identifier |
| SBO (fallback) | yes (`sbo_path`) | deepest term with >= 3 unplaced reactions | - | R9: bottom-up roll-up to the parent |

## Names

| box | label | `name_source` | `name_xref` |
|---|---|---|---|
| bin | the bin's names as written in the notes, joined by `.`; at level 2 only the part below the level-1 term | `notes` | `<notes key>:<code>`, e.g. `mapman:10.7` |
| EC | `EC <number> (<name>)`; name = text the notes give for that EC term, else the class name for a class or for a level-1 box | `ec` | `ec-code:<number>` |
| pathway | `KEGG pathway map00908`, `<db> pathway <id>` - the identifier itself, the SBML has no name for it | `xref` | `kegg.pathway:map00908` |
| GO | the GO id | `go` | the GO id |
| SBO | the SBO term name | `sbo` | the SBO id |
| localisation | compartment name, plus ` boundary` | `compartment` | compartment id |
| `(other)` | label of the opened term + ` (other)` when a sibling box of that term exists | as the term | the opened term |
| leftover | SBO name, or empty | `sbo` / `none` | SBO id / empty |

Basis: bins, EC, pathway, GO -> `function`; SBO -> `reaction_type`; localisation -> `localisation`;
unnamed leftover -> `mixed`. No name is derived from words in node names or free text.

## Parameters

| parameter | value | how chosen |
|---|---|---|
| `M_MIN` | 3 | Justified, not fitted: two reactions are a pair, three is the smallest set that is a group. Counted in reactions for function boxes (so one enzyme with its substrates is never a box - the "four-digit EC" failure), in nodes for localisation boxes. Set before the first run and never changed. Sensitivity (2, 5) is in NOTES.md. |
| `PAGE` | 12 | Justified: a 3 x 4 grid of labelled boxes is what one landscape page holds legibly. Applies to level 1 and to the sub-boxes of any box. Set before the first run and never changed. Only `skm_pss` reaches it. Sensitivity (9, 16) in NOTES.md. |
| `S_MAX` | 0.5 | Justified by the definition of a catch-all: a box holding "most of the model" holds more than half. Also used for "most of the box" in R13. Set before the first run and never changed; 1/3 gives the same group counts and largest level-1 share on all six models (NOTES.md). |
| `CLONE_MIN_BOXES` | 3 | Justified: a species touching two boxes is one visible link between them; from three boxes on its edges cross the page. The value was never changed; the unit was changed once during development (from finest function boxes to level-1 boxes), which lowered the number of cloned species (e_coli_core 34 -> 25). |
| ontology order (R4a) | bins, EC, pathway, GO; SBO last | Justified: first the ontologies whose terms have a tree in the bundle (so depth can be chosen) and names (so a reader can read the box), then flat bare identifiers; SBO says how a reaction happens, not what it is for, so it is the fallback. |
| process ontologies (R3, R8) | bins, pathway, GO | Justified: a process term is shared by the reactions and molecules of the process, an EC number belongs to one reaction. |
| reaction weight (R12) | 1 / number of species in the reaction | Derived from the data (`edges.tsv`), no free value. |
| EC class names | 7 names | The only vocabulary embedded in `strategy.py` (oxidoreductases ... translocases). Everything else comes from the feature tables. |

No parameter, and no branch of the code, depends on which model is processed.

## Where species go, and when they are cloned

1. A species with an own term whose box exists lives in that box (R11) - an annotated protein is drawn where its
   annotation says, even if most of its reactions are elsewhere.
2. Every other species lives where the weighted majority of its reactions is (R12); main-metabolite lists are
   respected when the model has them.
3. Inside its level-1 box the same two rules give the level-2 box (R13); the localisation box is not chosen, it is
   the species' own `loc_space` / `loc_kind` (R14).
4. A species has exactly one home. It is cloned only when its reactions lie in three or more level-1 boxes
   (R15): then one copy is listed for each function box, outside the home level-1 box, that holds one of its
   reactions. Species touching two level-1 boxes are not cloned, so links between two boxes stay visible as edges.
