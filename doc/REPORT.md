# Which grouping makes the clearest pathway map?

BioHackathon Europe 2026, project 22 (*Open Pathway Maps*). Eight grouping strategies, each computed from the
SBML alone, run on the same six models, scored by one rulebook of 27 rules that three other agents wrote
without seeing any strategy. Everything here is reproducible from this repository; `README.md` lists the
commands.

## 1. What was compared

| strategy | top level | inside |
|---|---|---|
| `type_location` | reaction / signalling type from the SBO hierarchy (EC class where SBO does not divide) | localisation per compartment, with boundary groups |
| `location_type` | localisation per compartment | reaction type, EC class as fallback |
| `ontology` | annotation trees: notes bins (MapMan), EC, pathway cross-references, GO; SBO fallback | next level of the same tree, then localisation |
| `topology_leiden` / `_infomap` / `_walktrap` | communities of the species-reaction graph | re-run of the same algorithm inside oversized groups |
| `cascade` | flow anchors (signalling) and routes/cycles (metabolism) | stage, then localisation |
| `hybrid` | function term where annotated, connectivity where not, size-bounded | finer function, reaction type, localisation |
| `baseline` | the project's current grouping (`../grouping_output`) | flat |

Panel: `brenda_merged` (your merged BRENDA model, 178 species / 218 reactions), `e_coli_core`, `R-HSA-77289`,
and the Stress Knowledge Map `skm_ja`, `skm_sa`, `skm_pss` (961 / 856). Held out for the scale check:
Human-GEM (`Human2.sbml`, 21,122 nodes, 146 curated subsystems).

All 48 strategy outputs (8 x 6) are contract-valid and byte-identical on a second run, verified here, and each
is also written as nested SBML Level 3 `groups` that layout-benchmark's grouped layouts can read.

## 2. Ranking

Mean over the six models, merged-rulebook score (`scoring/scores_summary.csv`; per model and per rule in
`scoring/scores.csv`):

| rank | strategy | mean | worst model | structure A | graph B | localisation C | naming D | reference E |
|---|---|---|---|---|---|---|---|---|
| 1 | `hybrid` | 0.622 | skm_ja (0.56) | 0.95 | 0.74 | 0.33 | 0.87 | 0.38 |
| 2 | `location_type` | 0.611 | skm_pss (0.35) | 0.91 | 0.52 | 1.00 | 0.94 | 0.46 |
| 3 | `topology_infomap` | 0.558 | skm_ja (0.37) | 0.94 | 0.92 | 0.04 | 0.79 | 0.28 |
| 4 | `topology_leiden` | 0.558 | e_coli_core (0.51) | 0.95 | 0.93 | 0.04 | 0.78 | 0.24 |
| 5 | `ontology` | 0.551 | skm_pss (0.36) | 0.92 | 0.57 | 0.26 | 0.89 | 0.52 |
| 6 | `type_location` | 0.532 | skm_pss (0.29) | 0.85 | 0.54 | 0.61 | 0.91 | 0.45 |
| 7 | `cascade` | 0.474 | skm_sa (0.31) | 0.88 | 0.69 | 0.25 | 0.81 | 0.35 |
| 8 | `topology_walktrap` | 0.467 | e_coli_core (0.33) | 0.88 | 0.91 | 0.04 | 0.80 | 0.28 |
| 9 | `baseline` | 0.307 | skm_pss (0.19) | 0.67 | 0.81 | 0.08 | 0.43 | 0.56 |

The current grouping is last on every model, and by a wide margin: 0.307 against 0.47-0.62 for the eight
strategies. Its weaknesses are the ones you saw in the report page - unnamed groups (naming block 0.43 against
0.78-0.94) and, on `skm_pss`, 167 singleton groups out of 346.

**The top two are not separated by this evidence.** `hybrid` leads by 0.011; on `score_core` (the same score
without the reference block) `location_type` leads, 0.646 to 0.630. Treat ranks 1-2 as a tie, and ranks 3-4
(`topology_infomap` 0.5585, `topology_leiden` 0.5580) as indistinguishable.

Rule verdicts behind the ranking (counted over the 6 models x 26 scored rules, `scores.csv`): `location_type`
fails fewest, 17, then `hybrid` 27, the three topology variants 24 each (with more rules not applying),
`type_location` 33, `ontology` 36, `cascade` 35, and the baseline 25 of a smaller applicable set.

## 3. How much the ranking depends on the weights

`scoring/sensitivity.py` re-scores every row with the scorer's own arithmetic under 11 weightings (equal
blocks, each block heavy, each block removed) and with each of the 26 scored rules dropped in turn:

| strategy | rank as published | best | worst | first place in |
|---|---|---|---|---|
| `hybrid` | 1 | 1 | 2 | 33 of 38 |
| `location_type` | 2 | 1 | 5 | 5 of 38 |
| `topology_infomap` | 3 | 3 | 6 | 0 of 38 |
| `topology_leiden` | 4 | 2 | 6 | 0 of 38 |
| `ontology` | 5 | 3 | 6 | 0 of 38 |
| `type_location` | 6 | 3 | 7 | 0 of 38 |
| `cascade` | 7 | 7 | 8 | 0 of 38 |
| `topology_walktrap` | 8 | 6 | 8 | 0 of 38 |
| `baseline` | 9 | 9 | 9 | 0 of 38 |

`hybrid` is first in 33 of 38 variants and `location_type` in the other 5 - the five where localisation or the
curated reference carries the most weight. No other strategy ever takes first place, and the baseline is last
in all 38. So the choice *between the top two* is a choice about what you weight; the gap to the rest, and to
the current grouping, is not.

## 4. The finding that matters most

**No strategy reproduces curated pathway labels better than drawing one box per compartment.** Rule M25 scores
agreement with the held-out labels as skill against the better of two trivial groupings (one box per
compartment, everything in one box). It is negative in 19 of the 34 rows where a reference exists. On the
held-out Human-GEM, run unchanged on 21,122 nodes:

| | `hybrid` | `location_type` | one box per compartment |
|---|---|---|---|
| F1 against the 146 subsystems | 0.340 | 0.495 | 0.489 |
| skill over the trivial grouping | -0.293 | +0.010 | - |
| rulebook score | 0.617 | 0.339 | - |

`location_type` matches the compartment baseline because at that scale it very nearly *is* the compartment
baseline (largest top-level group 6,617 of 21,122 nodes). `hybrid` keeps its size bounds (12 groups,
1,175-2,391 nodes each, largest 11%) and scores higher on the rulebook, but agrees with curation less.

Reading: curated subsystems are not recoverable from SBML annotation and topology alone on these models. What
the strategies can do reliably is make a map *readable* - bounded boxes, named from identifiers, localisation
drawn as a cell. Claiming they recover curation would not be supported by these numbers.

## 5. Where each strategy fails

Each track stated its own failures in `strategies/<name>/NOTES.md`; the ones confirmed by the scores:

- **`type_location`** (your proposal). Readable types, good names (naming 0.91), but on the SKM models the
  type boxes cut pathways: on `skm_pss` it is the lowest of the eight strategies (0.287, with only the baseline below it at 0.191). A type is a property of one reaction, so
  the top level never groups a cascade.
- **`location_type`**. Perfect localisation block (1.00, by construction) and the best names, but one
  compartment swallows the map: largest top-level group 0.86 (R-HSA-77289), 0.73 (`skm_pss`), 0.31 on
  Human-GEM. Nine innermost boxes stay above 50 nodes because no type attribute divides SKM reactions further.
- **`ontology`**. The best reference/biology block of the eight (0.52) and recognisable names where MapMan
  bins exist ("Phytohormone action.jasmonic acid.biosynthesis"). Below the top EC classes, labels are bare
  identifiers; its page rule dropped 19 top-level boxes on `skm_pss` (118 reactions re-placed).
- **`topology_*`**. Best graph block (0.91-0.93) - these are the groups the graph actually has - but
  localisation 0.04 (they declare none) and 89% of groups named after a hub member, which locates a module
  without saying what it does. Infomap and Leiden are indistinguishable here; Walktrap is clearly worse (0.467).
- **`cascade`**. Rank 7 of 9. Hubs cut textbook pathways (glycolysis lands in three top-level
  groups of `e_coli_core`), and route labels are member species.
- **`hybrid`**. Best structure block (0.95) and the only strategy whose bounds hold from 148 to 21,122 nodes,
  but its reaction-type level exists in only 3 of 6 models, 31 leaf groups fall below its own minimum size,
  and its reference agreement is the weakest of the two leaders.

## 6. Recommendation

For the SKM protocol, **function first, localisation inside, with size bounds** - the `hybrid` shape - with two
corrections taken from the other tracks:

1. Seed the function level from *every* controlled identifier, not just MapMan notes: EC, GO and pathway
   cross-references too (`ontology` does this; `hybrid` does not, which is why its level 1 falls back to
   connectivity on the three metabolic models).
2. Keep `location_type`'s localisation level as the innermost level everywhere, drawn as a cell from the
   containment tree. It costs nothing in readability and is the only level that a curator can check by eye.

Do not adopt a pure topology strategy for a published map: it scores best on the graph but names its boxes
after molecules. Do not adopt localisation-first alone: one compartment takes the whole map.

## 7. Limits

- **The revision round was not completed.** Scores are first-round only. Six revisions were dispatched and
  stopped before they delivered, at your request to wrap up.
- Determinism was verified by rerunning every strategy (48 of 48 byte-identical), but is recorded in the
  scorer as verified input rather than measured inside it.
- The scorer's test suite was first reported as 49 passing by the agent that wrote it, in its own flat layout;
  in this repository it initially ran 42 passed / 7 skipped because the tests looked for `panel/features`
  beside the script. With the paths corrected it runs 49 passed here. The three rule-proposal
  suites pass in full: 42, 48 and 67 tests.
- Thresholds: of the 27 rules, only the box-size limits, the 3-12 items-per-level range and the modularity
  range have a literature source; the rest are the rule authors' judgement, marked as such in `RULES.md`, and
  calibrated on synthetic groupings only. Several proved uninformative on this panel (M04, M05, M20, M27
  separate nothing between strategies) - stated, not retuned.
- The reference is thin: curated labels exist for 4 of 6 panel models, cover species only in the SKM models,
  and the baseline's own labels for `brenda_merged` and `skm_pss` came from that same curation, so its
  reference score there is not a prediction (flagged `ref_contaminated` in `scores.csv`).
- `topology` ran three algorithms as asked; `cascade`, `ontology` and `hybrid` delivered one variant each, so
  their parameter sensitivity is reported in their own notes rather than scored.
- Two tracks (`location_type`, `ontology`) fixed parameters after inspecting results on this panel, so their
  thresholds are not blind to it. No strategy ever saw the curated labels.

## 8. Files

- `CONTRACT.md` - what a strategy reads and writes; `common/` - feature extractor, compartment containment
  tree, contract validator, SBML groups writer.
- `strategies/<name>/` - code, `RULES.md` (the if-then rules), `NOTES.md` (per-model statistics and failures),
  `out/<model>/` (groups.tsv, groups_meta.tsv, `<model>_grouped.xml`).
- `rules/` - the three rule proposals with their own metrics and tests; `scoring/RULES.md` - the merged
  rulebook; `scoring/score_grouping.py`, `scores.csv`, `scores_summary.csv`, `sensitivity/`.
- `examples/index.html` - every model under every strategy, one renderer.
- `heldout/` - the curated labels (scoring only) and the Human-GEM feature tables.
