"""How much of the ranking depends on one choice: alternative block weights, and each rule dropped in turn.

The merged rulebook combines 27 rules in five blocks with weights fixed before scoring (A structure 0.25,
B graph 0.20, C localisation 0.20, D naming 0.15, E reference 0.20). This script rebuilds each row's rule
table from the per-rule utilities and verdicts in scores.csv and re-scores it with score_grouping.score -
the same arithmetic as the scorer, including the core-structure factor and the gates - under

  * alternative block weights, including equal weights and each block removed;
  * each rule removed in turn (leave-one-out).

It reports every strategy's best, worst and published rank, and how often each takes first place.

Run:
    python scoring/sensitivity.py --scores scoring/scores.csv --out scoring/sensitivity
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "rules"), str(HERE.parent / "common")]
import score_grouping as SG  # noqa: E402

BLOCKS = list(SG.BLOCK_WEIGHTS)
PUBLISHED = "as published (A.25 B.20 C.20 D.15 E.20)"
SCHEMES = {
    PUBLISHED: dict(SG.BLOCK_WEIGHTS),
    "equal blocks": {b: 0.2 for b in BLOCKS},
    "structure heavy (A.40)": {"A": 0.40, "B": 0.15, "C": 0.15, "D": 0.15, "E": 0.15},
    "graph heavy (B.40)": {"A": 0.15, "B": 0.40, "C": 0.15, "D": 0.15, "E": 0.15},
    "localisation heavy (C.40)": {"A": 0.15, "B": 0.15, "C": 0.40, "D": 0.15, "E": 0.15},
    "naming heavy (D.40)": {"A": 0.15, "B": 0.15, "C": 0.15, "D": 0.40, "E": 0.15},
    "reference heavy (E.40)": {"A": 0.15, "B": 0.15, "C": 0.15, "D": 0.15, "E": 0.40},
    "without block A (structure)": {"A": 0.0, "B": 0.25, "C": 0.25, "D": 0.20, "E": 0.30},
    "without block B (graph)": {"A": 0.30, "B": 0.0, "C": 0.25, "D": 0.20, "E": 0.25},
    "without block C (localisation)": {"A": 0.30, "B": 0.25, "C": 0.0, "D": 0.20, "E": 0.25},
    "without block D (naming)": {"A": 0.30, "B": 0.25, "C": 0.25, "D": 0.0, "E": 0.20},
    "without block E (reference)": {"A": 0.30, "B": 0.25, "C": 0.25, "D": 0.20, "E": 0.0},
}


def tables(df: pd.DataFrame) -> dict:
    """(strategy, model) -> the rule table the scorer would have built, from scores.csv."""
    out = {}
    for _, row in df.iterrows():
        rows = [dict(rule=r.id, block=r.block, weight=r.weight,
                     utility=row.get(f"{r.id}_utility"), verdict=row.get(f"{r.id}_verdict"))
                for r in SG.RULES]
        out[(row["strategy"], row["model"])] = pd.DataFrame(rows)
    return out


def means(tabs: dict, weights: dict, drop_rules=()) -> pd.Series:
    """Mean score per strategy under these weights, using the scorer's own score()."""
    saved = dict(SG.BLOCK_WEIGHTS)
    SG.BLOCK_WEIGHTS.update({b: weights.get(b, 0.0) for b in saved})
    try:
        per = {}
        for (strategy, _model), table in tabs.items():
            per.setdefault(strategy, []).append(SG.score(table, drop_rules=drop_rules)["score"])
    finally:
        SG.BLOCK_WEIGHTS.update(saved)
    return pd.Series({k: sum(v) / len(v) for k, v in per.items()}).sort_values(ascending=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scores", default="scoring/scores.csv")
    ap.add_argument("--out", default="scoring/sensitivity")
    a = ap.parse_args(argv)
    df = pd.read_csv(a.scores)
    missing = [r.id for r in SG.RULES if f"{r.id}_utility" not in df.columns]
    if missing:
        raise SystemExit(f"scores.csv has no utility column for {missing}")
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tabs = tables(df)

    variants = {name: (w, ()) for name, w in SCHEMES.items()}
    for r in SG.RULES:
        if r.block != "gate":
            variants[f"without {r.id}"] = (dict(SG.BLOCK_WEIGHTS), (r.id,))
    rows = []
    for name, (w, drop) in variants.items():
        s = means(tabs, w, drop)
        kind = "weights" if name in SCHEMES else "rule dropped"
        for rank, (strategy, value) in enumerate(s.items(), start=1):
            rows.append(dict(variant=name, kind=kind, strategy=strategy, score=round(value, 4), rank=rank))
    ranks = pd.DataFrame(rows)
    ranks.to_csv(out / "sensitivity_ranks.csv", index=False, lineterminator="\n")

    piv = ranks.pivot_table(index="strategy", columns="variant", values="rank")
    summary = pd.DataFrame({
        "rank_as_published": piv[PUBLISHED].astype(int),
        "rank_best": piv.min(axis=1).astype(int),
        "rank_worst": piv.max(axis=1).astype(int),
        "times_first": (piv == 1).sum(axis=1).astype(int),
        "variants": piv.shape[1],
    }).sort_values("rank_as_published")
    summary.to_csv(out / "sensitivity_summary.csv", lineterminator="\n")
    check = means(tabs, SCHEMES[PUBLISHED])
    print(summary.to_string())
    print("\nfirst place over all variants:")
    print(ranks[ranks["rank"] == 1].groupby("strategy").size().sort_values(ascending=False).to_string())
    print("\npublished weights, recomputed mean score:")
    print(check.round(4).to_string())
    print(f"\nwrote {out}/sensitivity_ranks.csv and sensitivity_summary.csv")


if __name__ == "__main__":
    main()
