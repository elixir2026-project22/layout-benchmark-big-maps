#!/usr/bin/env python
"""score_grouping.py - one scorer for the merged rulebook (RULES.md), BioHackathon Europe 2026, project 22.

27 rules M01..M27 merged from the three blind proposals in rules/. Nothing is reimplemented: every metric
is read from metrics_graph.metrics(), metrics_readability.metrics() or metrics_biology.metrics(), except
M05 (fan_max12_share), which the published fan_ok_share mixes with the lower bound that is M06 here. The
"module" field of every rule says where its number comes from.

Every metric is NaN - never 0 - when its rule does not apply; an NaN rule is dropped and the weights of its
block are rescaled. Deterministic: no randomness, no model-specific branch, rows sorted before writing.

    python score_grouping.py --features panel/features/skm_ja --grouping strategies/hybrid/skm_ja
    python score_grouping.py --features panel/features/skm_ja --grouping strategies/hybrid/skm_ja \
        --reference heldout/skm_ja.reference.tsv --json skm_ja_hybrid.json
    python score_grouping.py --features panel/features/skm_ja --grouping strategies/hybrid/skm_ja \
        --rerun strategies/hybrid/skm_ja_run2
    python score_grouping.py --all --panel panel/features --strategies strategies --out scores.csv
    python -m pytest test_score_grouping.py -q
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_UP = os.path.dirname(HERE)
for _p in (os.path.join(HERE, "rules"), os.path.join(HERE, "common"), HERE,
           os.path.join(_UP, "rules"), os.path.join(_UP, "common")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import metrics_biology as MB          # noqa: E402
import metrics_graph as MG            # noqa: E402
import metrics_readability as MR      # noqa: E402

NAN = float("nan")
UTILITY_AT_WARN = 0.25
UTILITY_AT_PASS = 0.75
GATE_RULES = ("M01", "M02", "M03", "M04")      # structural gates (union of the two proposals' gate sets)
GATE_CAP = 0.40
FAIL_IF_ONE_GROUP = ("M09", "M11", "M12")      # a leaf level with one group fails these, it is not "na"

# --------------------------------------------------------------------------------- the rulebook
# id, source rule, module, metric key, direction, worst, warn, pass, best, block, weight inside block
R = collections.namedtuple("R", "id src module metric direction worst warn passes best block weight")

RULES = [
    # block A - structure and granularity
    R("M01", "RD02_BOXSIZE", "read", "oversize_node_share",        "lower", 0.70, 0.30, 0.10, 0.0, "A", 3.0),
    R("M02", "RD03_SINGLETON", "read", "singleton_share",          "lower", 0.50, 0.25, 0.10, 0.0, "A", 3.0),
    R("M03", "RD01_TOP", "read", "top_groups",                     "range", (1, 40), (2, 20), (3, 12), None, "A", 3.0),
    R("M04", "G2", "graph", "largest_share_leaf",                  "lower", 1.00, 0.50, 1.0 / 3.0, 0.0, "A", 3.0),
    R("M05", "RD08_FAN", "own", "fan_max12_share",                 "higher", 0.40, 0.70, 0.90, 1.0, "A", 2.0),
    R("M06", "G10", "graph", "single_child_share",                 "lower", 1.00, 0.20, 0.00, 0.0, "A", 1.0),
    R("M07", "RD11_BALANCE", "read", "size_balance",               "higher", 0.10, 0.25, 0.50, 1.0, "A", 1.0),
    R("M08", "RD12_DEPTH", "read", "depth_excess",                 "lower", 2.00, 1.00, 0.00, 0.0, "A", 1.0),
    # block B - graph structure
    R("M09", "G1", "graph", "cut_w_leaf",                          "lower", 1.00, 0.50, 0.25, 0.0, "B", 3.0),
    R("M10", "G4", "graph", "connected_node_share_leaf",           "higher", 0.00, 0.80, 0.95, 1.0, "B", 2.0),
    R("M11", "G5", "graph", "Q_eff",                               "higher", 0.00, 0.40, 0.70, 1.0, "B", 2.0),
    R("M12", "G6", "graph", "weak_node_share_leaf",                "higher", 0.00, 0.50, 0.80, 1.0, "B", 1.0),
    R("M13", "G7", "graph", "edge_span_ratio",                     "lower", 1.00, 2.0 / 3.0, 1.0 / 3.0, 0.0, "B", 2.0),
    R("M14", "G11", "graph", "clone_misuse_share",                 "lower", 1.00, 0.30, 0.10, 0.0, "B", 1.0),
    # block C - localisation fidelity
    R("M15", "LOC-1", "bio", "loc_expressed",                      "higher", 0.00, 0.90, 1.00, 1.0, "C", 0.25),
    R("M16", "LOC-2", "bio", "loc_repeat",                         "inverse", 0.00, 2.00, 1.00, 1.0, "C", 0.25),
    R("M17", "LOC-3", "bio", "loc_boundary",                       "higher", 0.00, 0.80, 1.00, 1.0, "C", 0.20),
    R("M18", "LOC-4", "bio", "loc_nesting",                        "higher", 0.00, 0.90, 1.00, 1.0, "C", 0.15),
    R("M19", "LOC-5", "bio", "loc_span_ok",                        "higher", 0.00, 0.90, 1.00, 1.0, "C", 0.15),
    # block D - naming and reuse
    R("M20", "RD04_NAMED", "read", "named_share",                  "higher", 0.20, 0.60, 0.90, 1.0, "D", 2.0),
    R("M21", "FAIR-2", "bio", "name_lift",                         "higher", 0.00, 0.21, 0.61, 1.0, "D", 2.0),
    R("M22", "RD07_DISTINCT", "read", "distinct_share",            "higher", 0.50, 0.80, 0.95, 1.0, "D", 1.0),
    R("M23", "FAIR-3", "bio", "name_grounded",                     "higher", 0.00, 0.50, 0.80, 1.0, "D", 2.0),
    R("M24", "FAIR-4", "bio", "stranded_poor",                     "lower", 1.00, 0.30, 0.10, 0.0, "D", 1.5),
    # block E - biology
    R("M25", "BIO-1", "bio", "ref_skill",                          "higher", 0.00, 0.05, 0.25, 1.0, "E", 0.6),
    R("M26", "BIO-2", "bio", "hom",                                "higher", 0.00, 0.05, 0.21, 1.0, "E", 0.4),
    # gate
    R("M27", "FAIR-1", "bio", "fair_deterministic",                "higher", 0.00, 1.00, 1.00, 1.0, "gate", 0.0),
]
RULE_BY_ID = {r.id: r for r in RULES}
BLOCK_WEIGHTS = {"A": 0.25, "B": 0.20, "C": 0.20, "D": 0.15, "E": 0.20}
BLOCK_NAMES = {"A": "structure", "B": "graph", "C": "localisation", "D": "naming", "E": "biology"}
REFERENCE_RULE = "M25"     # excluded from score_core so models with and without a reference compare

# Diagnostics kept in scores.csv: the eight dropped rules plus the per-rule diagnostics of RULES.md.
DIAGNOSTICS = [
    ("graph", "singleton_group_share_leaf"), ("graph", "readable_node_share_leaf"),
    ("graph", "max_siblings"), ("graph", "ref_bcubed_f1_best"), ("graph", "ref_cut_w"),
    ("read", "long_edge_share"), ("read", "cloned_species_share"), ("read", "backed_share"),
    ("read", "comp_redraw"),
    ("graph", "n_nodes"), ("graph", "n_species"), ("graph", "n_reactions"), ("graph", "n_edges"),
    ("graph", "depth"), ("graph", "n_groups_L1"), ("graph", "n_groups_leaf"),
    ("graph", "largest_share_L1"), ("graph", "cut_raw_leaf"), ("graph", "cut_w_L1"),
    ("graph", "Q_w_best"), ("graph", "Q_ref"), ("graph", "edge_span_w"), ("graph", "edge_span_null"),
    ("graph", "residual_node_share"), ("graph", "hub_degree_cutoff"), ("graph", "n_hubs"),
    ("graph", "clone_idle_share"), ("graph", "clone_nonhub_share"), ("graph", "n_clones"),
    ("graph", "disconnected_group_share_leaf"), ("graph", "connected_nohub_node_share_leaf"),
    ("graph", "conductance_mean_leaf"),
    ("read", "n_groups"), ("read", "box_size_max"), ("read", "box_over100_node_share"),
    ("read", "singleton_node_share"), ("read", "named_node_share"), ("read", "fan_ok_share"),
    ("read", "fan_max"), ("read", "only_child_share"), ("read", "l1_cross_share"),
    ("read", "cut_edge_share"), ("read", "modularity_L1"), ("read", "clone_demand"),
    ("read", "clone_load"), ("read", "comp_redraw_max"), ("read", "loc_box_repeat"),
    ("read", "mixed_box_node_share"), ("read", "depth_mean"), ("read", "depth_needed"),
    ("read", "top_all_localisation"), ("read", "sibling_dup_share"), ("read", "same_as_ancestor_share"),
    ("read", "name_recall_mean"), ("read", "ref_scatter"), ("read", "ref_mix"),
    ("bio", "ref_f1"), ("bio", "ref_precision"), ("bio", "ref_recall"), ("bio", "ref_best_level"),
    ("bio", "ref_f1_L1"), ("bio", "ref_f1_L2"), ("bio", "ref_f1_L3"), ("bio", "ref_baseline_f1"),
    ("bio", "ref_f1_one_group"), ("bio", "ref_f1_singletons"), ("bio", "ref_f1_compartments"),
    ("bio", "ref_labelled_share"), ("bio", "hom_L1"), ("bio", "hom_sbo"), ("bio", "hom_ec"),
    ("bio", "hom_go"), ("bio", "hom_class"), ("bio", "hom_compartments"), ("bio", "hom_gain"),
    ("bio", "annotated_share"), ("bio", "loc_pure_leaf"), ("bio", "loc_boundary_in_interior"),
    ("bio", "loc_boxes_mean"), ("bio", "loc_boxes_max"), ("bio", "loc_span_share"),
    ("bio", "name_support"), ("bio", "name_vacuous_share"), ("bio", "name_resolvable"),
    ("bio", "named_group_share"), ("bio", "poor_share"), ("bio", "stranded_annotated"),
    ("bio", "fair_sorted"), ("bio", "fair_rule_share"), ("bio", "coverage"),
]


# --------------------------------------------------------------------------------- helpers
def _isnan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def read_tsv(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def fan_max12_share(nodes, edges, comps, groups, meta, fan_max: int = 12) -> float:
    """M05: share of parent groups with at most `fan_max` child boxes (upper bound of RD08 only).

    Uses metrics_readability's own box convention (a localisation group and its sibling @boundary group
    are one box; a parent that also holds loose nodes counts one extra box). NaN for a flat grouping.
    """
    m = MR._Model(nodes, edges, comps, groups, meta)
    parents = [g for g in m.groups if m.children.get(g)]
    if not parents:
        return NAN
    return sum(m.n_boxes(g) <= fan_max for g in parents) / len(parents)


def unit(rule: R, v: float) -> float:
    """Metric -> utility in [0, 1]: 0 at the rule's worst value, 0.25 at the warn threshold, 0.75 at the
    pass threshold, 1 at the best value, linear in between (RULES.md 4.1). NaN stays NaN."""
    if _isnan(v):
        return NAN
    if rule.direction == "range":
        (zlo, zhi), (wlo, whi), (plo, phi) = rule.worst, rule.warn, rule.passes
        if plo <= v <= phi:
            return 1.0
        if v < plo:
            return float(np.interp(max(v, zlo), [zlo, wlo, plo], [0.0, UTILITY_AT_WARN, 1.0]))
        return float(np.interp(-min(v, zhi), [-zhi, -whi, -phi], [0.0, UTILITY_AT_WARN, 1.0]))
    worst, warn, ok, best = rule.worst, rule.warn, rule.passes, rule.best
    if rule.direction == "inverse":                       # loc_repeat: 1..inf, lower is better
        v, worst, warn, ok, best = 1.0 / max(v, 1.0), 0.0, 1.0 / warn, 1.0 / ok, 1.0
    elif rule.direction == "lower":
        v, worst, warn, ok, best = -v, -worst, -warn, -ok, -best
    xs, ys = [worst], [0.0]
    for x, y in ((warn, UTILITY_AT_WARN), (ok, UTILITY_AT_PASS), (best, 1.0)):
        if x > xs[-1] + 1e-12:
            xs.append(x)
            ys.append(y)
        else:
            ys[-1] = max(ys[-1], y)
    return float(np.interp(min(max(v, worst), best), xs, ys))


def verdict(rule: R, v: float) -> str:
    """'pass' | 'warn' | 'fail' | 'na', read from the metric and the rule's own thresholds."""
    if _isnan(v):
        return "na"
    if rule.direction == "range":
        (zlo, zhi), (wlo, whi), (plo, phi) = rule.worst, rule.warn, rule.passes
        if plo <= v <= phi:
            return "pass"
        if wlo <= v <= whi:
            return "warn"
        return "fail"
    if rule.direction in ("lower", "inverse"):
        return "pass" if v <= rule.passes + 1e-12 else "warn" if v <= rule.warn + 1e-12 else "fail"
    return "pass" if v >= rule.passes - 1e-12 else "warn" if v >= rule.warn - 1e-12 else "fail"


# --------------------------------------------------------------------------------- the scorer
def metrics(nodes, edges, comps, groups, meta, reference=None, rerun=None,
            determinism: float = NAN) -> dict:
    """All metrics of the three modules plus M05, in one flat dict with the module prefix stripped.

    Keys that two modules define under the same name (n_nodes, ref_f1 ...) keep the value of the module
    the merged rulebook reads them from, so the rule table is unambiguous: graph < readability < biology
    is the fill order, and each rule later reads its own module's dict.
    """
    mg = MG.metrics(nodes, edges, comps, groups, meta, reference=reference)
    mr = MR.metrics(nodes, edges, comps, groups, meta, reference=reference)
    mb = MB.metrics(nodes, edges, comps, groups, meta, reference=reference, rerun=rerun)
    if rerun is None and not _isnan(determinism):
        mb["fair_deterministic"] = float(determinism)      # verified input, see RULES.md section 5
    own = {"fan_max12_share": fan_max12_share(nodes, edges, comps, groups, meta)}
    return {"graph": mg, "read": mr, "bio": mb, "own": own}


def rule_table(mods: dict) -> pd.DataFrame:
    """One row per rule: metric value, verdict, utility, block and weight."""
    one_group = mods["graph"].get("n_groups_leaf", NAN) < 2
    rows = []
    for r in RULES:
        v = mods[r.module].get(r.metric, NAN)
        if r.metric == "top_groups" and mods["read"].get("top_all_localisation", 0.0) == 1.0 \
                and not _isnan(v) and v < r.passes[0]:
            # RD01: a top level made only of compartment boxes has no lower bound (the model gives it)
            vd, ut = "pass", 1.0
        else:
            vd, ut = verdict(r, v), unit(r, v)
        if _isnan(v) and one_group and r.id in FAIL_IF_ONE_GROUP:
            vd, ut = "fail", 0.0                            # a leaf level with one group separates nothing
        rows.append({"rule": r.id, "source_rule": r.src, "metric": r.metric, "module": r.module,
                     "value": float(v) if not _isnan(v) else NAN, "verdict": vd, "utility": ut,
                     "block": r.block, "weight": r.weight})
    return pd.DataFrame(rows)


def block_scores(table: pd.DataFrame, drop_rules=()) -> dict:
    """Weighted mean utility inside each block, over the rules that apply."""
    out = {}
    for b in BLOCK_WEIGHTS:
        sub = table[(table.block == b) & (~table.rule.isin(drop_rules)) & table.utility.notna()]
        out[b] = float((sub.utility * sub.weight).sum() / sub.weight.sum()) if len(sub) else NAN
    return out


def score(table: pd.DataFrame, drop_rules=()) -> dict:
    """One number per grouping and model (RULES.md section 4.3): block mean, core-structure factor,
    gate cap, determinism gate."""
    blocks = block_scores(table, drop_rules)
    num = sum(BLOCK_WEIGHTS[b] * s for b, s in blocks.items() if not _isnan(s))
    den = sum(BLOCK_WEIGHTS[b] for b, s in blocks.items() if not _isnan(s))
    mean = num / den if den else NAN
    util = dict(zip(table.rule, table.utility))
    vd = dict(zip(table.rule, table.verdict))
    core = [util[r] for r in GATE_RULES if not _isnan(util[r])]
    factor = 0.5 + 0.5 * min(core) if core else 1.0
    total = mean * factor
    if any(vd[r] == "fail" for r in GATE_RULES):
        total = min(total, GATE_CAP)
    if vd["M27"] == "fail":
        total = 0.0
    return {"blocks": blocks, "mean_unweighted": mean, "gate_factor": factor, "score": total}


def score_grouping(features_dir, grouping_dir, reference=None, rerun_dir=None,
                   determinism: float = NAN) -> dict:
    """Score one (model, grouping) pair. Returns metrics, the rule table and the scores."""
    nodes = read_tsv(os.path.join(features_dir, "nodes.tsv"))
    edges = read_tsv(os.path.join(features_dir, "edges.tsv"))
    comps = read_tsv(os.path.join(features_dir, "compartments.tsv"))
    groups = read_tsv(os.path.join(grouping_dir, "groups.tsv"))
    meta = read_tsv(os.path.join(grouping_dir, "groups_meta.tsv"))
    ref = read_tsv(reference) if isinstance(reference, str) else reference
    rerun = None
    if rerun_dir:
        rerun = (read_tsv(os.path.join(rerun_dir, "groups.tsv")),
                 read_tsv(os.path.join(rerun_dir, "groups_meta.tsv")))
    mods = metrics(nodes, edges, comps, groups, meta, reference=ref, rerun=rerun,
                   determinism=determinism)
    table = rule_table(mods)
    full = score(table)
    core = score(table, drop_rules=(REFERENCE_RULE,))
    return {"metrics": mods, "table": table, "score": full["score"], "score_core": core["score"],
            "blocks": full["blocks"], "gate_factor": full["gate_factor"],
            "mean_unweighted": full["mean_unweighted"],
            "verdicts": dict(zip(table.rule, table.verdict)),
            "utilities": dict(zip(table.rule, table.utility))}


def score_from_tables(nodes, edges, comps, groups, meta, reference=None, rerun=None,
                      determinism: float = NAN) -> dict:
    """Same as score_grouping but from data frames (used by the tests)."""
    mods = metrics(nodes, edges, comps, groups, meta, reference=reference, rerun=rerun,
                   determinism=determinism)
    table = rule_table(mods)
    full = score(table)
    core = score(table, drop_rules=(REFERENCE_RULE,))
    return {"metrics": mods, "table": table, "score": full["score"], "score_core": core["score"],
            "blocks": full["blocks"], "gate_factor": full["gate_factor"],
            "verdicts": dict(zip(table.rule, table.verdict)),
            "utilities": dict(zip(table.rule, table.utility))}


# --------------------------------------------------------------------------------- panel driver
# Provenance of the scored rows, read from the bundle's own README files, not from any metric:
#   * determinism: byte-identical reruns were verified outside this scorer for every strategy
#     (RULES.md section 5). The baseline is a converted artefact of an earlier pipeline with no
#     verified rerun, so its M27 stays NaN -> "na" and gates nothing.
#   * contaminated reference: baseline/README.md states that for brenda_merged and skm_pss the baseline
#     was built from the curated labels held out here, so its M25 is not a prediction. The flag only
#     annotates the output row; it changes no metric and no threshold.
DETERMINISM_VERIFIED_INPUT = 1.0
CONTAMINATED_REFERENCE = {("baseline", "brenda_merged"), ("baseline", "skm_pss")}


def _row(strategy, model, res, contaminated, determinism_source) -> dict:
    row = collections.OrderedDict()
    row["strategy"] = strategy
    row["model"] = model
    row["score"] = res["score"]
    row["score_core"] = res["score_core"]
    for b, name in BLOCK_NAMES.items():
        row["block_" + name] = res["blocks"][b]
    row["gate_factor"] = res["gate_factor"]
    row["verdict_string"] = "".join({"pass": "P", "warn": "w", "fail": "F", "na": "-"}[res["verdicts"][r.id]]
                                    for r in RULES)
    row["ref_contaminated"] = int(contaminated)
    row["determinism_source"] = determinism_source
    tbl = res["table"].set_index("rule")
    for r in RULES:
        row[r.id + "_" + r.metric] = tbl.at[r.id, "value"]
    for r in RULES:
        row[r.id + "_verdict"] = tbl.at[r.id, "verdict"]
    for r in RULES:
        row[r.id + "_utility"] = tbl.at[r.id, "utility"]
    for mod, key in DIAGNOSTICS:
        row["diag_" + key] = res["metrics"][mod].get(key, NAN)
    return row


def run_all(panel_dir, strategies_dir, baseline_dir=None, heldout_dir=None, out="scores.csv",
            summary=None) -> pd.DataFrame:
    """Score every strategy (and the baseline, if given) on every panel model. Deterministic order."""
    models = sorted(d for d in os.listdir(panel_dir) if os.path.isdir(os.path.join(panel_dir, d)))
    sources = []
    if strategies_dir and os.path.isdir(strategies_dir):
        sources += [(s, os.path.join(strategies_dir, s)) for s in
                    sorted(d for d in os.listdir(strategies_dir)
                           if os.path.isdir(os.path.join(strategies_dir, d)))]
    if baseline_dir and os.path.isdir(baseline_dir):
        sources.append(("baseline", baseline_dir))
    rows = []
    for name, root in sources:
        for model in models:
            gdir = os.path.join(root, model)
            if not (os.path.isfile(os.path.join(gdir, "groups.tsv"))
                    and os.path.isfile(os.path.join(gdir, "groups_meta.tsv"))):
                continue
            ref = None
            if heldout_dir:
                cand = os.path.join(heldout_dir, model + ".reference.tsv")
                ref = cand if os.path.isfile(cand) else None
            verified = name != "baseline"
            res = score_grouping(os.path.join(panel_dir, model), gdir, reference=ref,
                                 determinism=DETERMINISM_VERIFIED_INPUT if verified else NAN)
            rows.append(_row(name, model, res, (name, model) in CONTAMINATED_REFERENCE,
                             "verified_input" if verified else "not_verified"))
            print(f"{name:>18s}  {model:<14s} score {res['score']:.3f}  core {res['score_core']:.3f}  "
                  f"{rows[-1]['verdict_string']}")
    df = pd.DataFrame(rows).sort_values(["strategy", "model"]).reset_index(drop=True)
    if out:
        df.to_csv(out, index=False, lineterminator="\n")
        print(f"\nwrote {out}  ({len(df)} rows x {len(df.columns)} columns)")
    if summary:
        summarise(df).to_csv(summary, index=False, lineterminator="\n")
        print(f"wrote {summary}")
    return df


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """One row per strategy: mean and minimum score, block means, and the rules it never passes."""
    rows = []
    for name, sub in df.groupby("strategy", sort=True):
        row = collections.OrderedDict()
        row["strategy"] = name
        row["n_models"] = len(sub)
        row["mean_score"] = sub["score"].mean()
        row["min_score"] = sub["score"].min()
        row["mean_score_core"] = sub["score_core"].mean()
        blind = sub["score_core"].where(sub["ref_contaminated"] == 1, sub["score"])
        row["mean_score_refblind"] = blind.mean()      # contaminated rows scored without M25
        row["worst_model"] = sub.loc[sub["score"].idxmin(), "model"]
        row["best_model"] = sub.loc[sub["score"].idxmax(), "model"]
        for b, nm in BLOCK_NAMES.items():
            row["block_" + nm] = sub["block_" + nm].mean()
        row["mean_edge_span_ratio"] = sub["M13_edge_span_ratio"].mean()
        fails, never = [], []
        for r in RULES:
            v = sub[r.id + "_verdict"]
            if (v == "fail").any():
                fails.append(r.id)
            if not (v == "pass").any() and not (v == "na").all():
                never.append(r.id)
        row["rules_failed_somewhere"] = "|".join(fails)
        row["rules_never_passed"] = "|".join(never)
        row["ref_contaminated_models"] = "|".join(sorted(sub.loc[sub["ref_contaminated"] == 1, "model"]))
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["mean_score", "min_score"], ascending=False)
    out.insert(0, "rank", range(1, len(out) + 1))
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------------- command line
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Score a grouping against the merged rulebook (RULES.md).")
    ap.add_argument("--features", help="panel/features/<model>")
    ap.add_argument("--grouping", help="directory with groups.tsv and groups_meta.tsv")
    ap.add_argument("--reference", help="heldout/<model>.reference.tsv")
    ap.add_argument("--rerun", help="second run of the same strategy; sets M27 by byte comparison")
    ap.add_argument("--json", help="write the full metric and rule table here")
    ap.add_argument("--determinism", type=float, default=NAN,
                    help="set M27 from verified input (1 = verified byte-identical reruns)")
    ap.add_argument("--all", action="store_true", help="score every strategy on every panel model")
    ap.add_argument("--panel", default="panel/features")
    ap.add_argument("--strategies", default="strategies")
    ap.add_argument("--baseline", default="baseline")
    ap.add_argument("--heldout", default="heldout")
    ap.add_argument("--out", default="scores.csv")
    ap.add_argument("--summary", default="scores_summary.csv")
    a = ap.parse_args(argv)

    if a.all:
        run_all(a.panel, a.strategies, a.baseline, a.heldout, a.out, a.summary)
        return 0
    if not (a.features and a.grouping):
        ap.error("--features and --grouping are required unless --all is given")
    det = a.determinism
    if a.rerun:
        same = MB.same_bytes(a.grouping, a.rerun)
        print(f"rerun byte-identical: {same}")
    res = score_grouping(a.features, a.grouping, reference=a.reference, rerun_dir=a.rerun,
                         determinism=det)
    tbl = res["table"]
    print(f"{'rule':<5s} {'source':<14s} {'metric':<28s} {'value':>10s}  {'verdict':<8s} {'utility':>7s}")
    for _, r in tbl.iterrows():
        v = "na" if _isnan(r.value) else f"{r.value:.4f}"
        u = "na" if _isnan(r.utility) else f"{r.utility:.3f}"
        print(f"{r.rule:<5s} {r.source_rule:<14s} {r.metric:<28s} {v:>10s}  {r.verdict:<8s} {u:>7s}")
    print("\nblocks: " + "  ".join(
        f"{BLOCK_NAMES[b]} {('na' if _isnan(s) else f'{s:.3f}')}" for b, s in res["blocks"].items()))
    print(f"gate factor {res['gate_factor']:.3f}")
    print(f"score {res['score']:.3f}   score_core (without M25) {res['score_core']:.3f}")
    if a.json:
        payload = {"score": res["score"], "score_core": res["score_core"], "blocks": res["blocks"],
                   "gate_factor": res["gate_factor"], "verdicts": res["verdicts"],
                   "utilities": res["utilities"],
                   "rules": json.loads(tbl.to_json(orient="records")),
                   "metrics": {m: {k: (None if _isnan(v) else v) for k, v in d.items()}
                               for m, d in res["metrics"].items()}}
        with open(a.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        print(f"wrote {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
