# Documentation

| document | what it settles |
|---|---|
| [`../README.md`](../README.md) | what the repository is, the two stages, the options, and the commands that reproduce everything |
| [`../CONTRACT.md`](../CONTRACT.md) | what a grouping strategy reads and writes, and the seven rules every strategy keeps |
| [`../layout/README.md`](../layout/README.md) | the layout stage: the protocol in five steps, the eight options, the files a run writes |
| [`../layout/PROTOCOL.md`](../layout/PROTOCOL.md) | the placement rule by rule, each checkable against `positions.tsv` and `boxes.tsv`, and the determinism contract |
| [`REPORT.md`](REPORT.md) | the result of the grouping comparison: what was compared, the ranking, how much it depends on the weights |
| [`../scoring/RULES.md`](../scoring/RULES.md) | the merged rulebook: the 27 rules, their metrics and their verdicts |
| [`../scoring/SCORING_NOTES.md`](../scoring/SCORING_NOTES.md) | how the rulebook was merged from the three proposals, and what the scorer does with each rule |
| `../rules/rules_graph.md`, `rules_biology.md`, `rules_readability.md` | the three independent rule proposals the rulebook was merged from, each with its metrics module and tests |
| `../strategies/<name>/RULES.md` | the numbered rules one strategy's `rule` column refers to; `NOTES.md` next to it is its working record |
| [`../baseline/README.md`](../baseline/README.md) | where the baseline grouping comes from, and why its agreement with the reference is not a prediction |
| [`../FILESTRUCTURE.md`](../FILESTRUCTURE.md) | what every file and folder of the repository is |
| [`../PUSH.md`](../PUSH.md) | the commands that push this repository to its GitHub remote |

`REPORT.md` was written while the work was in one repository together with the models; it refers to
numbers in `scoring/scores.csv` and `scoring/scores_summary.csv`, both of which are kept here.
