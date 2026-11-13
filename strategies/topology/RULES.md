# Topology strategy: rules

**What it is.** Groups are modules of the species-reaction graph, found by community detection after the
currency species have been cut out. Three algorithms are delivered as separate variants
(`topology_leiden`, `topology_infomap`, `topology_walktrap`). They differ in rule R4 only; every other
rule, parameter and line of code is shared (`strategy.py --algorithm leiden|infomap|walktrap`).

**What it is not.** No rule reads a compartment, so no group has basis `localisation` and
`compartment_id` is always empty. Group `basis` is `topology` for Leiden and Walktrap and `flow` for
Infomap.

Rules are applied in the order R1 to R9. The `rule` column of `groups.tsv` holds the rule that fixed the
node's deepest group: `R4`, `R5a`, `R5b`, `R5c`, `R6`, `R7` or `R8`.

## Rules

### R1 - the graph
- **If** a row of `edges.tsv` names a reaction and a species that are both in `nodes.tsv`,
  **then** the two nodes are linked. One link per (reaction, species) pair, weight 1, whatever the role and
  however many roles the species has in that reaction. Stoichiometry is not used.
- **Direction** (used by Infomap only; Leiden and Walktrap use the same links without direction):
  reactant = species -> reaction; product = reaction -> species; modifier = species -> reaction (the
  modifier acts on the reaction and is not produced by it). **If** the reaction is reversible, **then** its
  reactant and product links run both ways; modifier links stay one-way.
- **Reads:** `edges.tsv` `reaction_id`, `species_id`, `role`; `nodes.tsv` `kind`, `reversible`.
- **Fallback:** a participant row with an id missing from `nodes.tsv` is skipped. An empty `reversible`
  counts as irreversible. A node without any link is a part of its own and ends under R5c.

### R2 - currency candidates
- **If** a species has degree >= d*, **then** it is a currency candidate. d* is the smallest degree that
  fewer than one species of the model would reach if the model's links were spread over its species at
  random: the smallest d with n x P(Poisson(mean species degree) >= d) < 1, n = number of species.
- **Reads:** the species degree from R1 (equal to `nodes.tsv` `degree`).
- **Fallback:** none needed, every species has a degree. A model in which no species reaches d* has no
  currency species.

### R3 - glue test (candidate becomes currency, or stays as a hub)
- Candidates are taken by falling degree (ties: `node_id`). A reaction is **stranded** by a candidate
  when, once the candidate and the species already cut are removed, the reaction's connected component
  holds less than half of the nodes its component has in the full graph.
- **If** fewer than half of the candidate's reactions are stranded, **then** the species is a currency
  species: it is cut out of the graph for R4-R7 and placed by R8.
  **Else** it stays in the graph: it is the centre of a module, not glue between modules.
- **Reads:** the links of R1.
- **Fallback:** none needed.

### R4 - modules (the only rule that differs between the variants)
- **Then** the graph without currency species is partitioned:
  - `leiden`: Leiden optimisation of modularity (resolution 1) on the undirected links, iterated until no
    further improvement, seed 1;
  - `infomap`: Infomap on the directed links of R1, best of 10 trials, seed 1;
  - `walktrap`: Walktrap on the undirected links, random walks of 4 steps, dendrogram cut where modularity
    is highest.
- A node whose top-level part is the algorithm's own community, untouched by R5 and R6, gets rule `R4`.
- **Reads:** the links of R1 without the currency species of R3.
- **Fallback:** a node set without any link gives one part per node (handled by R5).

### R5 - merge parts that are too small
Applied to the parts of R4, smallest part first (ties: lowest `node_id`), until no part changes.
- **R5a If** a part has fewer than S_MIN nodes and has direct links to other parts, **then** it joins the
  part with the most direct links to it (ties: the smaller part, then the lowest `node_id`).
- **R5b If** it has no direct link to another part but shares currency species with one, **then** it
  joins the part with the largest weight through currency species. That weight is the projection of each
  currency species c onto its d reactions in these parts: each pair of its reactions weighs 1/(d - 1), so a
  reaction keeps the weight 1 that its link to c had (the weighting Newman introduced for collaboration
  networks, 2001).
- **R5c If** it has neither, **then** it stays as it is: a separate small component of the model.
- **Reads:** the links of R1; the currency species of R3.
- **Fallback:** R5c is the fallback.

### R6 - at most K_MAX groups under one parent
- **If** there are more than K_MAX parts, **then** the two connected parts with the largest modularity
  gain are merged, and this repeats until K_MAX parts remain or no two parts are connected.
  Connection weight w(i,j) = direct links + weight through currency species (as in R5b).
  Gain = w(i,j)/m - s(i) s(j)/(2 m^2), with m the total weight and s the weight attached to a part
  (ties: the smaller merged size, then the lowest `node_id`s). The nodes of the smaller part get rule `R6`.
- **Reads:** the links of R1; the currency species of R3.
- **Fallback:** when no two parts are connected the rule stops, and the parent keeps more than K_MAX
  groups. This needs more than K_MAX separate components; no panel model has that many.

### R7 - next level for groups that are too large
- **If** a group without children has more than S_MAX nodes (currency species not counted) and its level
  is below 3, **then** R4, R5 and R6 are applied to the subgraph of its nodes, with the same algorithm and
  seed, and the parts become its children. Nodes get rule `R7` (or `R5a`, `R5b`, `R5c`, `R6` when that
  rule moved them at the new level).
- **Reads:** the links of R1 inside the group.
- **Fallback:** when the subgraph gives fewer than two parts the group stays a leaf, larger than S_MAX.
  The same holds for a group at level 3.

### R8 - currency species: home and copies
- **Then** each currency species is placed in the leaf group that holds the most of its reactions (ties:
  lowest group id) and inherits that group's L1..L3. Every other leaf group that holds one of its reactions
  is listed in `clone_in`. Rule `R8`.
- **Reads:** the links of R1; the groups of R4-R7.
- **Fallback:** a currency species always has reactions; should it have none in any group it goes to the
  first L1 group.

### R9 - names
Groups are named level by level, larger groups first. A term or a label is used once per model.
- **R9a If** an annotation term is carried by at least half of the group's reactions or by at least half of
  its species, by at least two members, and the group holds more than half of all nodes of the model that
  carry the term, **then** the group is named after it. Among several such terms the one with the fewest
  carriers in the model wins, then the deepest one in its hierarchy. Terms and their sources:

  | `name_source` | term | read from | label | `name_xref` |
  |---|---|---|---|---|
  | `notes` | classification code: a notes value of the form `1.2.3_name.name.name`, several separated by `;`; every prefix of the code is a term | `nodes.tsv` `notes`, any key except free text | the name path up to that depth (key and code when the name path has a different depth than the code) | `<notes key>:<code>`, e.g. `mapman:10.7.2` |
  | `sbo` | the SBO term of a reaction and each of its ancestors | `sbo_path`, `sbo_path_names` of reactions | the SBO name | `SBO:0000589` |
  | `ec` | EC number and each of its prefixes | `ec`, any node | `EC 2.7.1` | `ec-code:2.7.1` |
  | `go` | GO id | `go`, any node | the GO id (the tables carry no GO names) | `GO:0006635` |

- **R9b Else** the group is named after its hub (`name_source` = `member`, `name_xref` = the node id,
  label = the node's `name`): the species that takes part in the most reactions of the group (ties: the
  larger share of its own links inside the group, then `node_id`). When no species takes part in two
  reactions of the group, the reaction with the most participants in the group is the hub. Currency
  species are never hubs. A candidate whose name is empty or already used is passed over.
- **R9c Else** the group is unnamed (`name_source` = `none`, empty label).
- **Reads:** `nodes.tsv` `notes`, `sbo_path`, `sbo_path_names`, `ec`, `go`, `name`, `kind`.
- **Fallback:** R9b when no term qualifies or the node has no annotation; R9c when no member has a usable name.

### Group ids
`M01`, `M02`, ... at level 1, larger groups first (ties: lowest `node_id`); children `M01.01`,
`M01.01.01`.

## Parameters

No parameter is set per model. The panel informed the choice between candidate definitions for R2
(see the first row); it did not set any number.

| Parameter | Value | Rule | How the value was chosen |
|---|---|---|---|
| expected species beyond the degree bound (`CHANCE`) | 1 | R2 | Derived from the data: the bound d* comes from each model's own number of species and mean degree. "Fewer than one species by chance" is the boundary itself, not a fitted number. Resulting d*: brenda_merged 13, e_coli_core 11, R-HSA-77289 8, skm_pss 9, skm_ja 7, skm_sa 8. Tukey's fences on the same degree distribution were tried first and dropped: at Q3 + 3 IQR, E. coli core lost only cytosolic H+ and kept water and ATP in the graph; at Q3 + 1.5 IQR the fence was 3.5 in three models (IQR of 1 on integer degrees), which makes a species with four reactions an outlier. |
| "most of" (`MAJORITY`) | 0.5 | R3, R9 | Justified: a simple majority, used in four places (a reaction is stranded when it keeps less than half of its network; a candidate is a hub when at least half of its reactions are stranded; a term must cover at least half of the group; the group must hold more than half of the term). |
| smallest module (`S_MIN`) | 5 nodes | R5 | Justified: two reactions in a row with their three species (A -> B -> C) is the smallest thing that is a path and not a single reaction. |
| groups per page (`K_MAX`) | 12 | R6 | Justified: the largest qualitative colour palettes have 12 classes (ColorBrewer Paired, Set3), so every group on a page can have its own colour; 12 boxes also fill a 4 x 3 grid on a landscape page. |
| largest leaf group (`S_MAX`) | 60 nodes | R7 | Derived from the two values above: K_MAX x S_MIN. A group that could hold a full page of smallest modules is itself a page and gets a level of its own. |
| levels (`MAX_LEVEL`) | 3 | R7 | The contract allows L1 to L3. |
| link weight | 1 for every role | R1 | Justified: the simplest choice with no free number; stoichiometry would let one biomass reaction dominate. |
| Leiden resolution | 1 | R4 | Standard modularity; any other value would be a free number. |
| Walktrap walk length (`WALKTRAP_STEPS`) | 4 | R4 | The value of the original paper and the igraph default. NOTES.md gives the effect of 3 and 5. |
| Infomap trials (`INFOMAP_TRIALS`) | 10 | R4 | igraph default. |
| seed | 1 | R4 | Fixed for determinism. NOTES.md gives the spread over seeds 1 to 10. |

## How species are placed, and when clones are used

- **Species that are not currency species** are nodes of the graph like reactions. The algorithm places
  them (R4), and R5, R6 and R7 move them together with their part. They have one home and **no clones**:
  a link from such a species to a reaction in another group is a real connection between the two groups
  and stays visible as one.
- **Currency species** (R2 and R3) are not in the graph while modules are found. R8 gives each one home,
  the leaf group with most of its reactions, and a clone in every other leaf group where it reacts.
  Clones are used for currency species only.
- Reactions are never cloned.

## Contract rule 7 - what happens to a node that lacks an attribute

| Missing | Effect | Rule id |
|---|---|---|
| no SBO term, no EC, no GO, no notes | none on grouping; the node offers no term for R9a | group named by R9b |
| no name (reactions of the SKM models) | the node cannot be the hub of R9b; the next candidate is taken | R9b, else R9c |
| no link at all | the node is a part of its own and stays so | R5c |
| part of a small separate component | the component stays one small group | R5c |
| one compartment only, or many | none: compartments are not read | - |
