# SCORING_NOTES.md - what was merged, what was dropped, what the numbers say

BioHackathon Europe 2026, project 22. Companion to `RULES.md` (the merged rulebook), `score_grouping.py`
(the one scorer), `test_score_grouping.py` (49 tests, all passing), `scores.csv` (52 rows x 184 columns)
and `scores_summary.csv`.

```
python -m pytest test_score_grouping.py -q
python score_grouping.py --all --panel panel/features --strategies strategies \
    --baseline baseline --heldout heldout --out scores.csv --summary scores_summary.csv
python score_grouping.py --features panel/features/skm_ja --grouping strategies/hybrid/skm_ja \
    --reference heldout/skm_ja.reference.tsv --determinism 1 --json skm_ja_hybrid.json
```

## 1. What was merged

35 proposed rules -> **27 kept** (M01-M27), **8 dropped** as duplicates of a kept rule, **8 conflicts**
stated and decided (`RULES.md` section 2). The rulebook was finished and the weights fixed before any
strategy was scored; no threshold or weight was changed afterwards.

Dropped: G3 (singletons, covered by RD03 which exempts one-node compartments), G8 (box size 3-50, covered
by RD02 above and RD03 below), G9 (siblings at every level, covered by RD01 at the top and RD08 inside),
G12 (reference BCubed F1, covered by BIO-1 which is chance-corrected and handles multi-label references),
RD05 (long arcs, covered by G1 and G7, and defeated by cloning everything), RD06 (share of cloned species,
covered by G11's data-derived hub test), RD09 (backed names per group, covered by FAIR-3 per node),
RD10 (compartment redraw inferred from box purity, covered by LOC-2 counting declared boxes).
All eight are still computed and written to `scores.csv` as `diag_*` columns, so nothing is lost but the
double counting.

The eight conflicts and their decisions: sibling upper bound 9 vs 12 (-> 12), singleton threshold
0.05/0.20 vs 0.10/0.25 (-> 0.10/0.25), relative vs absolute box size (-> keep both, M04 and M01),
reference agreement in or out of the score (-> in, block E weight 0.20, with `score_core` reported
beside it), name lift vs name distinctness (-> keep both, M22 down-weighted for the declared overlap),
three different metric-to-number mappings (-> the biology four-knot mapping for all), gate as a cap or as
a factor (-> both, factor then cap), determinism as a gate (-> kept, set from verified input).

Reuse: 26 of the 27 rules read a metric straight from `metrics_graph.py`, `metrics_readability.py` or
`metrics_biology.py`; the `module` column of `rule_table()` and the "From" line of every rule in
`RULES.md` name which. The single exception is **M05** (`fan_max12_share`), computed in
`score_grouping.py` because the published `fan_ok_share` folds in the lower bound that is M06 here;
it reuses `metrics_readability._Model` for the box convention.

## 2. Thresholds that proved uninformative on this panel

Findings, not corrections - no threshold was retuned (hard rule of the task). Counts are over the 48
strategy rows (8 strategies x 6 models); the baseline's 4 rows are listed separately where they differ.

| rule | metric | pass / warn / fail over 48 strategy rows | finding |
|---|---|---|---|
| **M04** | `largest_share_leaf` | 48 / 0 / 0 | **Uninformative.** No strategy puts more than a third of a model in one leaf box; the baseline passes too (52 of 52 rows). The rule still earns its place as a gate - it is the one that stops "everything in one box" - but on this panel it separates nothing and its weight of 3 in block A is spent. |
| **M20** | `named_share` | 48 / 0 / 0 | **Uninformative between strategies**, informative against the baseline (2 of 4 baseline rows fail). Every strategy names every group; the project's motivating failure (34 named groups of 346) is exactly what it catches. |
| **M05** | `fan_max12_share` | 35 / 0 / 0, 13 na | **Uninformative where it applies.** No grouping with a second level has a parent with more than twelve child boxes. |
| **M27** | `fair_deterministic` | 48 pass (verified input) | Inert by construction in this run - see section 5. |
| **M02** | `singleton_share` | 46 / 2 / 0 | Almost uninformative: two warns, no fail. Again it is the baseline (`skm_pss`, 0.48) that the rule was written for. |
| **M11, M13, M24, M08** | `Q_eff`, `edge_span_ratio`, `stranded_poor`, `depth_excess` | no fail in any row | Not uninformative - the utilities spread (`edge_span_ratio` 0.04-0.65, `Q_eff` pass 38 / warn 10) - but the *fail* zone is empty, so read them as ordinal values, not as verdicts. |
| **M26** | `hom` | 9 / 36 / 3 | The pass limit of 0.21 is nearly out of reach, as the biology proposal predicted: 36 of 48 rows warn. On the SKM models reactions carry only an SBO type. Treat `hom` and `hom_gain` as a ranking value there. |
| **M25** | `ref_skill` | 2 / 6 / 24 (16 na) | Harsh: only two strategy rows (cascade and ontology on `skm_pss`, 0.322) beat the best trivial grouping by the pass margin, and 24 rows score below the compartment baseline that `ref_skill` subtracts. This says the curated pathways of the panel are not recoverable blind, not that the strategies are broken; the proposal itself called the threshold uncalibrated. |

Rules that discriminate well on this panel: M15-M19 (localisation), M14 (clone misuse, 27 fails), M10
(connected boxes, 22 fails), M09, M16, M22, M23, M01.

## 3. The ranking

Mean over the six panel models; `score_core` is the same score without M25, `refblind` scores the two
contaminated baseline rows without M25 as well (`RULES.md` 4.4).

| rank | strategy | models | mean | min | core | refblind | best model | worst model |
|---|---|---|---|---|---|---|---|---|
| 1 | hybrid | 6 | **0.622** | 0.558 | 0.630 | 0.622 | brenda_merged 0.733 | skm_ja 0.558 |
| 2 | location_type | 6 | **0.611** | 0.349 | 0.646 | 0.611 | R-HSA-77289 0.812 | skm_pss 0.349 |
| 3 | topology_infomap | 6 | 0.558 | 0.371 | 0.570 | 0.558 | brenda_merged 0.717 | skm_ja 0.371 |
| 4 | topology_leiden | 6 | 0.558 | 0.509 | 0.573 | 0.558 | R-HSA-77289 0.612 | e_coli_core 0.509 |
| 5 | ontology | 6 | 0.551 | 0.358 | 0.588 | 0.551 | R-HSA-77289 0.694 | skm_pss 0.358 |
| 6 | type_location | 6 | 0.532 | 0.287 | 0.573 | 0.532 | R-HSA-77289 0.708 | skm_pss 0.287 |
| 7 | cascade | 6 | 0.474 | 0.311 | 0.482 | 0.474 | R-HSA-77289 0.654 | skm_sa 0.311 |
| 8 | topology_walktrap | 6 | 0.467 | 0.325 | 0.477 | 0.467 | R-HSA-77289 0.618 | e_coli_core 0.325 |
| 9 | baseline | 4 | 0.307 | 0.191 | 0.293 | 0.293 | R-HSA-77289 0.390 | skm_pss 0.191 |

Ranks 3 and 4 are separated by 0.0005 on the mean (infomap 0.5585, leiden 0.5580) - indistinguishable;
had the means been equal the rulebook's tie-break on the minimum would have put `topology_leiden` first
(0.509 against 0.371). The three topology algorithms are separate rows throughout and do **not** come out
together: infomap and leiden at 0.558, walktrap at 0.467.

Block means (`scores_summary.csv`) show where each strategy loses:

| strategy | A structure | B graph | C localisation | D naming | E biology |
|---|---|---|---|---|---|
| hybrid | 0.955 | 0.741 | 0.332 | 0.868 | 0.383 |
| location_type | 0.908 | 0.525 | **1.000** | 0.936 | 0.455 |
| topology_infomap | 0.940 | 0.923 | **0.044** | 0.790 | 0.280 |
| topology_leiden | 0.951 | 0.927 | **0.043** | 0.784 | 0.239 |
| ontology | 0.919 | 0.574 | 0.255 | 0.891 | 0.524 |
| type_location | 0.849 | 0.538 | 0.607 | 0.915 | 0.454 |
| cascade | 0.877 | 0.691 | 0.255 | 0.811 | 0.349 |
| topology_walktrap | 0.879 | 0.915 | 0.044 | 0.798 | 0.277 |
| baseline | 0.674 | 0.808 | 0.085 | 0.428 | 0.557 |

The pattern is the one both the graph and the biology proposals predicted for their own viewpoints: the
three pure topology strategies declare no localisation group at all (`loc_expressed` = 0 on every model,
`loc_repeat` therefore `na`), so block C collapses to about 0.04 for them, while `location_type` is
perfect there (1.000) and pays on block B (0.525). `hybrid` ranks first because it is the only strategy
that is strong in A, respectable in B (0.741) and not empty in C (0.332); it wins by not failing a block,
not by leading one. The weights that produce this were fixed before scoring, and two of the five block
weights (B and C, both 0.20) decide the top of the table: at 0.30/0.10 the order of the first three would
change. That sensitivity belongs to the reader, not to the data.

Gates bit 13 of the 52 rows (9 strategy rows, all 4 baseline rows): capped at 0.40 by a failed structural
gate. The failures are M01 (a box over 50 nodes: `cascade`, `location_type`, `ontology`, `type_location`
on the two largest SKM models, `topology_walktrap` on `e_coli_core`, the baseline on `brenda_merged`) and
M03 (top-box count outside 3-20: the baseline on three models, `cascade` on two).

## 4. Caveats that belong to this ranking

1. **The baseline is not a peer.** It exists for 4 of 6 models, is a converted one-level artefact of an
   earlier pipeline, and for `brenda_merged` and `skm_pss` it was built from the very labels held out
   here. Those two rows are flagged `ref_contaminated = 1`; its M25 values there (0.252 and 0.332) are the
   two highest of the whole table and are **not predictions**. Its mean is reported twice, 0.307 with M25
   and 0.293 without; it is last either way.
2. **Block E is thin.** M25 applies to 32 of 48 strategy rows and labels only species on the three SKM
   models; M26's pass zone is nearly unreachable there. Compare `score` and `score_core` before reading
   any gap under 0.03 as real - the top two differ by 0.011 on the mean and swap on `score_core`
   (location_type 0.646 vs hybrid 0.630).
3. **The localisation block is a compliance test, not a quality test.** A strategy that simply declares
   one localisation group per compartment gets block C for free; one that groups by topology alone gets
   nothing, however good its boxes are. Whether that is the right trade depends on the renderer question
   in `RULES.md` 8.3, which the examples gallery should settle.
4. **M10 (connected boxes) fails 22 of 48 strategy rows**, and all 22 belong to the five
   annotation-driven or localisation-driven strategies (location_type, ontology and type_location 5 rows
   each, hybrid 4, cascade 3); no topology row fails it. That is by design - a box of one enzyme class is
   not a connected subgraph - but it means block B partly measures "is this a topology strategy".
5. **No threshold comes from a user study of nested pathway maps.** The ranking should be re-read with the
   pass limits moved by a third before any of it is called settled.
6. **One run per strategy.** Determinism is verified input (section 5), and no revision round was scored
   (section 6).
7. **The two large SKM models hold the scores down.** Per-model score ranges over all nine rows:
   R-HSA-77289 0.390-0.812, e_coli_core 0.292-0.737, brenda_merged 0.356-0.733, skm_ja 0.355-0.694,
   skm_pss 0.191-0.614, skm_sa 0.311-0.587. The two lowest ceilings are `skm_sa` (0.587) and `skm_pss`
   (0.614); `skm_pss` has the lowest floor (0.191, the baseline; 0.287 among strategies) and is where four
   of the nine rows (location_type, ontology, type_location and the baseline) have their minimum, and it
   is the only model on which a structural gate fails for five different groupings (those four plus
   cascade). Ranking by the mean, as prescribed, keeps each model to one sixth of the weight.

## 5. Determinism (M27): verified input

The bundle contains one run per strategy, so this scorer could not observe a second run. Byte-identical
reruns for all 8 strategies x 6 models were **verified by the project lead** before the bundle was
assembled; M27 is recorded from that fact, not measured here. Every strategy row carries
`fair_deterministic = 1` with `determinism_source = verified_input`; the baseline rows carry NaN and
`determinism_source = not_verified`, and M27 is `na` for them, so the gate neither helps nor hurts the
baseline. `--rerun <dir>` compares two output directories byte for byte and records
`determinism_source = measured`; no row of `scores.csv` was produced that way.

Contract validation was run separately (`python common/contract.py validate panel/features/<model>
<grouping>`): all 48 strategy outputs and all 4 baseline outputs print `valid`.

## 6. Deviations

- **No revision round.** The plan foresaw returning each strategy its scorecard for one revision and
  re-scoring. This consolidation agent cannot rerun the six strategy agents, and the bundle holds a single
  run per strategy, so only first-round scores exist. `scores.csv` is therefore one block of first-round
  numbers with no "revised" column; the effect of tuning to the rules is not measured.
- **Determinism not measured**, recorded as verified input (section 5) as the task instructed.
- **FAIR-4's annotation-ablation test is not implemented** (it needs the strategies to be rerun); carried
  from the proposal as a stated limit under M24.
- **`heldout/human2`** (21,122 nodes, seventh model with its own reference) was **not scored**: it is not
  a panel model, no strategy produced a grouping for it, and there is nothing to score. Its reference is
  untouched.
- Metric values are those of the three reference implementations as shipped; none of their code was
  modified.
