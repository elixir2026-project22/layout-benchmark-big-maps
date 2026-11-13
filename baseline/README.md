# Baseline: the current grouping

The grouping the project's `pathway-grouping` step produced (`../../grouping_output/`), converted to the
contract of this repository (one flat level) so it is scored with the same script as the new strategies.

| model | source file | groups |
|---|---|---|
| brenda_merged | `grouping_output/brenda.groups.tsv` | 11 |
| e_coli_core | `grouping_output/e_coli_core_grouped.xml` | 28 |
| R-HSA-77289 | `grouping_output/R-HSA-77289_grouped.xml` | 23 |
| skm_pss | `grouping_output/skm_pss_grouped.xml` | 346 |

There is no baseline for skm_ja and skm_sa (the step was not run on them).

**Read with care:** for brenda_merged and skm_pss the baseline was built *from* the curated pathway labels
that this comparison holds out as the reference (`source` = curated / propagated in the pipeline's own
tables). Its agreement with the reference is therefore not a prediction and is not comparable with the
strategies, which never see those labels.
