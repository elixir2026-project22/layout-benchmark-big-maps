# type_location - rules

Top level (L1): the type of reaction or signalling step, read from the reaction's SBO term at a depth chosen
by rule R9, and refined by EC class where SBO does not divide. Second level (L2): cell localisation, one
group per space inside each type and an `@boundary` group for transport and membrane steps. No third level.

Run:

```
python strategy.py --features panel/features/MODEL --sbml panel/models/MODEL.xml --out out/MODEL
python common/contract.py validate panel/features/MODEL out/MODEL
python common/contract.py write-sbml panel/models/MODEL.xml out/MODEL --features panel/features/MODEL --strategy type_location -o out/MODEL/MODEL_grouped.xml
```

One program, no model-specific branch, standard library only. The `rule` column of `groups.tsv` holds the
rule that decided the node's L1 group: R1-R8 or R10 for a reaction, R11 or R12 for a species.

## Rules, in the order applied

### A. The type path of every reaction

A type path is a list of steps from general to specific: first the SBO ancestors of the reaction's term, then
(R8) its EC class and EC subclass.

**R1 - SBO term present.**
- Condition: the reaction's `sbo` is not empty.
- Action: the path is `sbo_path` (root > ... > term); step names from `sbo_path_names`.
- Reads: `nodes.tsv` `sbo`, `sbo_path`, `sbo_path_names`.
- Fallback: no SBO term -> the first of R2-R7 whose condition holds. A term that SBO does not know gives a
  one-step path and so a type of its own.

**R2 - no SBO term, the reaction moves matter between spaces.**
- Condition: `sbo` empty; at least one reactant and one product; reactants and products lie in two or more
  spaces (a membrane compartment counts as the space it encloses).
- Action: path ends at `SBO:0000185` translocation reaction.
- Reads: `nodes.tsv` `compartments`; `compartments.tsv` `kind`, `encloses_id`; `edges.tsv` `role`.
- Fallback: R3.

**R3 - no SBO term, a source or a sink.**
- Condition: `sbo` empty; no reactant or no product; no modifier.
- Action: path ends at `SBO:0000631` pseudoreaction.
- Reads: `edges.tsv` `role`.
- Fallback: R4.

**R4 - no SBO term, complex formation.**
- Condition: `sbo` empty; two or more different reactants; exactly one product; the product species is a
  non-covalent complex (`SBO:0000253` is in its `sbo_path`; protein complex and macromolecular complex are
  below that term).
- Action: path ends at `SBO:0000177` non-covalent binding.
- Reads: `edges.tsv` `role`; species `sbo_path`.
- Fallback: R5.

**R5 - no SBO term, complex dissociation.**
- Condition: `sbo` empty; exactly one reactant, which is a non-covalent complex; two or more different products.
- Action: path ends at `SBO:0000180` dissociation.
- Reads: as R4.
- Fallback: R6.

**R6 - no SBO term, reactants become products in one space.**
- Condition: `sbo` empty; at least one reactant and one product; R2-R5 do not hold.
- Action: path ends at `SBO:0000176` biochemical reaction.
- Reads: `edges.tsv` `role`.
- Fallback: R7.

**R7 - no SBO term, nothing else known.**
- Condition: `sbo` empty; R2-R6 do not hold (for example a product and a modifier, no reactant).
- Action: path ends at `SBO:0000375` process.
- Reads: `edges.tsv` `role`.
- Fallback: none needed; every reaction now has a path.

**R8 - EC class as two further steps.**
- Condition: the reaction's `ec` is not empty.
- Action: append the EC class (first digit) and the EC subclass (first two digits) to the path. With several
  EC numbers: the class most of them share, then the subclass most of those share; a tie goes to the lower
  number. A reaction whose L1 group is named by an EC class gets `rule` = R8.
- Reads: `nodes.tsv` `ec`.
- Fallback: no EC number -> the path ends at the SBO term.

### B. The depth of the type tree: L1

**R9 - choose the depth.** Start with one type that holds every reaction. Repeat (a) and (b) until neither applies.
- (a) Descend. Condition: all reactions of a type continue to the same next step and none stops at the
  current one. Action: the type moves to that step (a type is always named by the most specific step its
  reactions share; this also carries a type into an EC class or subclass that all its reactions share).
- (b) Split. Condition: a type holds more than `theta` of all reactions and can be divided. Take the largest
  such type (tie: lower step id). Action: its SBO subtypes with at least `n_min` reactions become types of
  their own. If the term has no such SBO subtype - or the type is what is left of a term after its SBO
  subtypes were taken out - its EC classes with at least `n_min` reactions become types instead. Everything
  else stays together (R10). A split that would give more than `k_max` types is not carried out.
- Reads: the type paths of R1-R8.
- Fallback: a type with no subtype of at least `n_min` reactions stays whole, even above `theta`.

**R10 - subtypes too small to stand alone.**
- Condition: in a split, a subtype has fewer than `n_min` reactions, or a reaction stops at the divided term.
- Action: these reactions form one group. It is named by the most specific step all of them share; when that
  is the divided term itself the label is "other <term>" with the term's identifier. A reaction whose own
  subtype was pooled into an "other" group gets `rule` = R10.
- Reads: the type paths.
- Fallback: none.

**R11 - home type of a species.** A species has no type; it lives with its reactions.
- Condition: the species takes part in at least one reaction.
- Action: L1 = the type that holds most of its reactions (reactant, product or modifier, each reaction
  counted once). Ties, in this order: the type with more of its reactions located in or on the species' own
  space; the type with more reactions that produce it; the lower group id.
- Reads: `edges.tsv`; `nodes.tsv` `loc_space`.
- Fallback: R12.

**R12 - a species in no reaction.**
- Condition: degree 0.
- Action: L1 = `S_unconnected`, an unnamed group (basis `mixed`, `name_source` `none`); L2 by R13-R15.
- Reads: `edges.tsv`.
- Fallback: none. (No such species in the panel.)

### C. Localisation: L2

**R13 - a node inside a space.**
- Condition: `loc_kind` is `inside` or `between`, and the model has more than one localisation cell (R15).
- Action: L2 = the group (L1 type, `loc_space`), `compartment_id` = `loc_space`, basis `localisation`,
  label = the compartment's name. A `between` node goes to its `loc_space`, the common parent of the spaces
  it connects.
- Reads: `nodes.tsv` `loc_space`, `loc_kind`; `compartments.tsv` `name`.
- Fallback: R15.

**R14 - a node on a boundary (transport across it, or a participant in its membrane).**
- Condition: `loc_kind` is `boundary`.
- Action: L2 = the group (L1 type, `loc_space`@boundary), `compartment_id` = `<loc_space>@boundary`.
  Label: the name of the membrane compartment that encloses the space when the model has one, otherwise
  "<space name> boundary". This holds for every type, so a binding step across the nuclear envelope sits
  on the nucleus boundary of the "non-covalent binding" type, not in a transport type.
- Reads: `nodes.tsv` `loc_space`, `loc_kind`; `compartments.tsv` `name`, `kind`, `encloses_id`.
- Fallback: none.

**R15 - one compartment only.**
- Condition: every node of the model has the same (`loc_space`, boundary or not): localisation cannot divide
  anything.
- Action: L2 is one more step of the type tree instead: one split as in R9(b) (no `theta` test) inside each
  L1 type, groups with basis `reaction_type`; the rest becomes "other <L1 label>". An L1 type with no
  subtype of at least `n_min` reactions has no L2; its nodes sit directly in L1. A species goes to the L2
  group that holds most of its reactions within its home type (tie: lower group id).
- Reads: the type paths; `nodes.tsv` `loc_space`, `loc_kind`.
- Fallback: none.

**R16 - clones.**
- Condition: a species takes part in reactions of a type other than its home type.
- Action: one copy in each such type, written to `clone_in`: the L2 group of that type for the species' own
  localisation cell if that group exists, otherwise the L1 group of that type. (Under R15: the L2 group of
  that type holding most of the species' reactions.) No clone inside the home type.
- Reads: `edges.tsv`; the groups made by R9-R15.
- Fallback: a species with reactions of one type only has no clone.

## Parameters

| parameter | value | used in | how the value was chosen |
|---|---|---|---|
| `theta` | 1/3 | R9(b): a divisible type holding more than this share of the reactions is split | Justified, then checked on the whole panel. A bound of 1/k guarantees at least k top-level boxes wherever the annotation allows a division (with fewer, one box would hold more than 1/k and be split); k = 3 is the smallest top level that is more than a pair of boxes. The obvious alternative, 1/2 ("no box holds most of the model"), was run on all six models: it changes two of them (table below) - skm_sa stops at three types of which one holds 58% of the nodes, and R-HSA-77289 keeps all oxidoreductases in one box. The value is one number for all models and was not fitted to any score. |
| `n_min` | 3 reactions | R9(b), R10, R15: a subtype with fewer reactions is pooled | Justified: three is the smallest group that is neither a single reaction nor a pair. Sensitivity for 2 and 5 below. |
| `k_max` | 11 types | R9(b): guard on the number of L1 groups | Justified: fewer than a dozen top-level boxes. It never binds on the panel (the outputs are byte-identical with `k_max` = 99; the largest L1 count is 10). |

Sensitivity, same program, all six models (L1 groups / L2 groups / largest L1 share of nodes):

| model | default | theta = 1/2 | n_min = 2 | n_min = 5 |
|---|---|---|---|---|
| brenda_merged | 7 / 28 / 0.29 | 7 / 28 / 0.29 | 7 / 29 / 0.29 | 7 / 23 / 0.29 |
| e_coli_core | 10 / 12 / 0.22 | 10 / 12 / 0.22 | 10 / 12 / 0.22 | 8 / 10 / 0.22 |
| R-HSA-77289 | 8 / 12 / 0.25 | 6 / 10 / 0.39 | 10 / 14 / 0.25 | 6 / 10 / 0.35 |
| skm_pss | 6 / 52 / 0.28 | 6 / 52 / 0.28 | 6 / 52 / 0.28 | 6 / 52 / 0.28 |
| skm_ja | 5 / 24 / 0.31 | 5 / 24 / 0.31 | 5 / 24 / 0.31 | 5 / 24 / 0.31 |
| skm_sa | 5 / 25 / 0.34 | 3 / 20 / 0.58 | 6 / 26 / 0.34 | 5 / 25 / 0.34 |

Conventions that are not thresholds: a membrane counts as the space it encloses (R2, from `compartments.tsv`);
majority vote over several EC numbers with ties to the lower number (R8); the tie-break order of R11.

## How species are placed, and when they are cloned

- A species has exactly one home: L1 by R11 (the type of most of its reactions), L2 by R13/R14 (its own
  compartment inside that type; a species in a membrane compartment goes to the `@boundary` group of the
  space the membrane encloses). Its L2 group therefore always has its own `loc_space` as compartment.
- A species-only localisation group can result: for example the cytosolic and nuclear forms of a protein
  whose only reactions are transport steps sit in "transport reaction / cytoplasm" and "transport reaction /
  nucleus", with the transport reactions in "transport reaction / nucleus boundary" between them.
- A species is cloned once into every other L1 type in which it takes part (R16), never within a type:
  inside one type the localisation groups are drawn as one cell, so a reaction on a boundary reaches its
  participants on both sides without a copy. With the clones, every reaction finds all its participants
  inside its own L1 box (checked on all six models: 100% of the participant edges).

## Names

| group | `basis` | `name_source` | `name_xref` | label |
|---|---|---|---|---|
| L1 type named by an SBO term | `reaction_type` | `sbo` | the SBO id | the SBO term name; "other <name>" for the pooled rest (R10) |
| L1 or L2 type named by an EC class or subclass | `reaction_type` | `ec` | `ec-code:<class>` or `ec-code:<class>.<subclass>` | "EC 1 Oxidoreductases", "EC 2.7 Transferring phosphorus-containing groups" |
| L2 localisation, inside | `localisation` | `compartment` | the compartment id | the compartment name |
| L2 localisation, boundary | `localisation` | `compartment` | the membrane compartment id if the model has one, else the space id | the membrane's name, else "<space name> boundary" |
| `S_unconnected` (R12) | `mixed` | `none` | empty | empty |

SBO names are those in `nodes.tsv`; for the fallback rules R2-R7 the six SBO ids and names are constants in
`strategy.py`, checked against `common/SBO_OBO.obo`. EC class and subclass names are the Enzyme Commission
nomenclature, a constant table in `strategy.py` keyed by the EC identifier; a subclass missing from the table
is labelled "EC x.y". No label is derived from the name of a molecule, a reaction or a model.
