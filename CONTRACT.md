# Contract for a grouping strategy

A strategy turns one SBML model into a nested grouping of all its species and reactions. Every strategy
starts from the same tables and writes the same two files, so strategies differ in grouping logic only.

## Input (per model, in `panel/features/<model>/`)

Produced by `common/sbml_features.py` from `panel/models/<model>.xml`; column definitions are in that
script's docstring.

- `nodes.tsv` - one row per species and reaction: compartment, SBO term with its ancestors, EC, GO,
  all annotation cross-references, the SBML notes as key/value pairs, degree, and the node's place in
  the cell (`loc_space`, `loc_kind`, `loc_outer`).
- `edges.tsv` - one row per participant: reaction, species, role (reactant, product, modifier),
  stoichiometry, the SBO term of the role.
- `compartments.tsv` - the compartment containment tree: `parent_id` is the space a compartment lies
  inside; a membrane has `encloses_id`, the space it bounds; `virtual` marks spaces added to keep the
  biology (see `common/compartments.py`).

You may also read `panel/models/<model>.xml` directly, but anything you take from it must be something
any SBML file could carry.

## Output (per model, in `strategies/<name>/out/<model>/`)

`groups.tsv` - one row per node of `nodes.tsv`, every species and reaction exactly once.

| column | meaning |
|---|---|
| `node_id` | id from `nodes.tsv` |
| `L1`, `L2`, `L3` | group id at each level. `L1` required; `L2`, `L3` optional, no gaps. A group has exactly one parent. |
| `rule` | id of the rule in your `RULES.md` that placed this node, e.g. `R3` |
| `clone_in` | optional, species only: `\|`-joined ids of further groups where a copy of the species is drawn |

`groups_meta.tsv` - one row per group id.

| column | meaning |
|---|---|
| `group_id` | unique over all levels |
| `level` | 1, 2 or 3 |
| `parent` | group id one level up; empty at level 1 |
| `label` | the group's name; empty exactly when `name_source` is `none` |
| `basis` | `reaction_type`, `localisation`, `function`, `topology`, `flow` or `mixed` |
| `name_source` | `sbo`, `ec`, `go`, `xref`, `notes`, `compartment`, `member` or `none` |
| `name_xref` | the identifier behind the label: `SBO:0000176`, `ec-code:2.7.1`, `GO:0006096`, `mapman:10.5`, a compartment id, or (for `member`) a node id |
| `compartment_id` | only for basis `localisation`: a `compartment_id`, or `<compartment_id>@boundary` |

Both files: tab-separated, UTF-8, `\n` line ends, rows sorted by id.

## Rules every strategy must keep

1. **One program for all models.** `python strategy.py --features <dir> --sbml <file> --out <dir>` with no
   model-specific branch, no list of ids or names of particular molecules, pathways or models.
2. **Deterministic.** Two runs give byte-identical files. Fix every seed; sort before iterating.
3. **From the SBML only.** No network, no external database, no file outside `panel/` and `common/`.
4. **Names come from identifiers, never guessed from words.** A group is named only from an SBO term, EC
   class, GO term, another cross-reference, a notes key, a compartment, or one member node. Otherwise it
   stays unnamed (`name_source` = `none`, empty label).
5. **Localisation is given, not chosen.** `nodes.tsv` places every node. A localisation group holds only
   nodes whose `loc_space` is its compartment; an `@boundary` group only nodes with `loc_kind` =
   `boundary` (transport across that space's boundary, or participants in its membrane). Nodes with
   `loc_kind` = `between` belong to their `loc_space` (the common parent of the spaces they connect).
   Localisation groups are drawn nested by `compartments.tsv` - extracellular around cytosol, organelles
   inside cytosol - so never encode nesting of compartments in your own levels.
6. **Every node is placed**, including species. State the rule that places species (a species shared by
   several groups has one home; further copies go in `clone_in`).
7. **Declare the fallback.** State what happens to a node that lacks the attribute your rule needs
   (no SBO term, no EC, one compartment only), and which rule id it gets.

## Check and export

```bash
python common/contract.py validate panel/features/<model> strategies/<name>/out/<model>
python common/contract.py write-sbml panel/models/<model>.xml strategies/<name>/out/<model> \
    --features panel/features/<model> --strategy <name> -o strategies/<name>/out/<model>/<model>_grouped.xml
```

`validate` must print `valid`. `write-sbml` writes the grouping as nested SBML Level 3 `groups`
(`kind="partonomy"`, `SBO:0000633`), readable by layout-benchmark's grouped layouts.
