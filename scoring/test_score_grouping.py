"""pytest for score_grouping.py (merged rulebook RULES.md).

Part 1  the utility mapping and the verdicts, at and between the thresholds of every rule.
Part 2  a hand-built toy model whose metric values are known by construction, including the
        not-applicable cases (every such metric must be NaN, never 0).
Part 3  the degenerate groupings on panel models: everything in one group, every node its own group,
        one group per compartment, random. None of them may rank first.

    python -m pytest test_score_grouping.py -q
"""
from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd
import pytest

import score_grouping as SG      # first: it puts rules/ and common/ on sys.path
import metrics_graph as MG

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = HERE if os.path.isdir(os.path.join(HERE, "panel", "features")) else os.path.dirname(HERE)
PANEL = os.path.join(_ROOT, "panel", "features")


def isnan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


# ============================================================ part 1: mapping and verdicts
def test_rulebook_has_27_rules():
    assert len(SG.RULES) == 27
    assert [r.id for r in SG.RULES] == [f"M{i:02d}" for i in range(1, 28)]


def test_every_rule_has_a_source_rule_and_a_module():
    for r in SG.RULES:
        assert r.src and r.module in ("graph", "read", "bio", "own")


def test_block_weights_sum_to_one():
    assert sum(SG.BLOCK_WEIGHTS.values()) == pytest.approx(1.0)


@pytest.mark.parametrize("rule_id", ["M01", "M02", "M09", "M13", "M24"])
def test_utility_knots_lower_is_better(rule_id):
    r = SG.RULE_BY_ID[rule_id]
    assert SG.unit(r, r.worst) == pytest.approx(0.0)
    assert SG.unit(r, r.warn) == pytest.approx(SG.UTILITY_AT_WARN)
    assert SG.unit(r, r.passes) == pytest.approx(SG.UTILITY_AT_PASS)
    assert SG.unit(r, r.best) == pytest.approx(1.0)


@pytest.mark.parametrize("rule_id", ["M07", "M10", "M11", "M20", "M21", "M23", "M25", "M26"])
def test_utility_knots_higher_is_better(rule_id):
    r = SG.RULE_BY_ID[rule_id]
    assert SG.unit(r, r.worst) == pytest.approx(0.0)
    assert SG.unit(r, r.warn) == pytest.approx(SG.UTILITY_AT_WARN)
    assert SG.unit(r, r.passes) == pytest.approx(SG.UTILITY_AT_PASS)
    assert SG.unit(r, r.best) == pytest.approx(1.0)


@pytest.mark.parametrize("rule_id", ["M06", "M08", "M15", "M17", "M18", "M19"])
def test_pass_equal_to_best_earns_full_utility(rule_id):
    r = SG.RULE_BY_ID[rule_id]
    assert SG.unit(r, r.passes) == pytest.approx(1.0)
    assert SG.unit(r, r.warn) == pytest.approx(SG.UTILITY_AT_WARN)


def test_utility_is_monotone_and_clipped():
    r = SG.RULE_BY_ID["M09"]                       # cut share, lower is better
    xs = [0.0, 0.1, 0.25, 0.4, 0.5, 0.8, 1.0]
    us = [SG.unit(r, x) for x in xs]
    assert us == sorted(us, reverse=True)
    assert SG.unit(r, 2.0) == pytest.approx(0.0)   # beyond the worst value
    assert SG.unit(r, -1.0) == pytest.approx(1.0)  # beyond the best value


def test_utility_of_nan_is_nan():
    for r in SG.RULES:
        assert isnan(SG.unit(r, float("nan")))


def test_loc_repeat_is_inverted():
    r = SG.RULE_BY_ID["M16"]
    assert SG.unit(r, 1.0) == pytest.approx(1.0)                  # one box per compartment
    assert SG.unit(r, 2.0) == pytest.approx(SG.UTILITY_AT_WARN)   # two equal halves
    assert SG.unit(r, 100.0) < 0.05
    assert SG.verdict(r, 1.0) == "pass" and SG.verdict(r, 2.0) == "warn" and SG.verdict(r, 2.5) == "fail"


def test_range_rule_top_groups():
    r = SG.RULE_BY_ID["M03"]
    for v in (3, 7, 12):
        assert SG.verdict(r, v) == "pass" and SG.unit(r, v) == pytest.approx(1.0)
    for v in (2, 13, 20):
        assert SG.verdict(r, v) == "warn"
    for v in (1, 21, 60):
        assert SG.verdict(r, v) == "fail"
    assert SG.unit(r, 1) == pytest.approx(0.0)
    assert SG.unit(r, 2) == pytest.approx(SG.UTILITY_AT_WARN)
    assert SG.unit(r, 20) == pytest.approx(SG.UTILITY_AT_WARN)
    assert SG.unit(r, 40) == pytest.approx(0.0)


@pytest.mark.parametrize("rule_id,below,at,above", [("M01", 0.05, 0.10, 0.31), ("M02", 0.05, 0.10, 0.26),
                                                    ("M09", 0.2, 0.25, 0.51), ("M24", 0.05, 0.10, 0.31)])
def test_verdicts_lower_is_better(rule_id, below, at, above):
    r = SG.RULE_BY_ID[rule_id]
    assert SG.verdict(r, below) == "pass"
    assert SG.verdict(r, at) == "pass"
    assert SG.verdict(r, above) == "fail"


def test_verdict_of_nan_is_na():
    for r in SG.RULES:
        assert SG.verdict(r, float("nan")) == "na"


# ============================================================ part 2: a toy model, known by construction
def _frame(rows, cols):
    return pd.DataFrame(rows, columns=cols).astype(str)


NODE_COLS = ["node_id", "kind", "name", "compartment", "compartments", "modifier_compartments",
             "loc_space", "loc_kind", "loc_outer", "sbo", "sbo_name", "sbo_path", "sbo_path_names",
             "ec", "go", "xrefs", "notes", "reversible", "n_reactants", "n_products", "n_modifiers",
             "degree", "boundary"]
COMP_COLS = ["compartment_id", "name", "canonical", "go_id", "kind", "parent_id", "encloses_id",
             "virtual", "depth", "source", "n_species"]
META_COLS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]


@pytest.fixture(scope="module")
def toy():
    """18 nodes in one compartment: six reactions r1..r6, twelve species s01..s12.

    Three blocks of six nodes each: {r1,r2,s01..s04}, {r3,r4,s05..s08}, {r5,r6,s09..s12}. Inside a block
    with reactions a, b and species p, q, u, v the edges are a-p, a-q, b-q, b-u, b-v, so the block is one
    connected piece (q bridges the two reactions) and no edge leaves it.
    """
    nodes, edges = [], []
    for i in range(1, 7):
        nodes.append([f"r{i}", "reaction", f"r{i}", "cyt", "cyt", "", "cyt", "inside", "", "SBO:0000176",
                      "biochemical reaction", "SBO:0000176", "biochemical reaction", "", "", "", "",
                      "false", "1", "1", "0", "2", "false"])
    for i in range(1, 13):
        nodes.append([f"s{i:02d}", "species", f"s{i:02d}", "cyt", "cyt", "", "cyt", "inside", "",
                      "SBO:0000247", "simple chemical", "SBO:0000247", "simple chemical", "", "", "", "",
                      "", "0", "0", "0", "1", "false"])
    pairs = []
    for blk in range(3):                       # reactions 2k+1, 2k+2 with species 4k+1 .. 4k+4
        a, b, p = 2 * blk + 1, 2 * blk + 2, 4 * blk + 1
        pairs += [(a, p), (a, p + 1), (b, p + 1), (b, p + 2), (b, p + 3)]
    for k, (r, s) in enumerate(pairs):
        edges.append([f"r{r}", f"s{s:02d}", "reactant" if k % 2 == 0 else "product", "1",
                      "SBO:0000010", "reactant"])
    comps = [["cyt", "cytosol", "cytosol", "GO:0005829", "space", "", "", "false", "0", "test", "12"]]
    return (_frame(nodes, NODE_COLS),
            _frame(edges, ["reaction_id", "species_id", "role", "stoichiometry", "sbo", "sbo_name"]),
            _frame(comps, COMP_COLS))


def _blocks_grouping():
    """The three blocks of six, one flat level, each box named after its own first member."""
    assign = {"g1": ["r1", "r2", "s01", "s02", "s03", "s04"],
              "g2": ["r3", "r4", "s05", "s06", "s07", "s08"],
              "g3": ["r5", "r6", "s09", "s10", "s11", "s12"]}
    rows = [[n, g, "", "", "R1", ""] for g, members in assign.items() for n in members]
    groups = _frame(sorted(rows), ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([[g, "1", "", f"box {g}", "topology", "member", members[0], ""]
                   for g, members in assign.items()], META_COLS)
    return groups, meta


def test_toy_blocks_values_known_by_construction(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    v = dict(zip(res["table"].rule, res["table"].value))
    assert v["M09"] == pytest.approx(0.0)              # every edge inside its box
    assert v["M04"] == pytest.approx(1.0 / 3.0)        # three equal boxes of six of eighteen
    assert v["M02"] == pytest.approx(0.0)              # no singleton
    assert v["M01"] == pytest.approx(0.0)              # no box over 50 nodes
    assert v["M07"] == pytest.approx(1.0)              # equal sizes: effective boxes = boxes
    assert v["M10"] == pytest.approx(1.0)              # every box one connected piece
    assert v["M12"] == pytest.approx(1.0)              # every box a weak community
    assert v["M20"] == pytest.approx(1.0)              # every box named
    assert res["verdicts"]["M04"] == "pass"            # 1/3 is inside the pass zone
    assert res["verdicts"]["M09"] == "pass"


def test_toy_not_applicable_metrics_are_nan_not_zero(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    v = dict(zip(res["table"].rule, res["table"].value))
    assert isnan(v["M03"]) and res["verdicts"]["M03"] == "na"   # 18 nodes: N <= 20 reads as one picture
    for rid in ("M05", "M06"):                     # flat grouping: no parent has a child
        assert isnan(v[rid]) and res["verdicts"][rid] == "na"
    for rid in ("M15", "M16", "M17", "M18", "M19"):   # one compartment
        assert isnan(v[rid]) and res["verdicts"][rid] == "na"
    assert isnan(v["M14"]) and res["verdicts"]["M14"] == "na"    # no clone_in used
    assert isnan(v["M25"]) and res["verdicts"]["M25"] == "na"    # no reference given
    assert isnan(v["M27"]) and res["verdicts"]["M27"] == "na"    # no rerun, no verified flag
    assert isnan(res["blocks"]["C"])                             # whole block not applicable


def test_toy_block_c_is_dropped_from_the_score(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    rest = {b: s for b, s in res["blocks"].items() if not isnan(s)}
    den = sum(SG.BLOCK_WEIGHTS[b] for b in rest)
    mean = sum(SG.BLOCK_WEIGHTS[b] * s for b, s in rest.items()) / den
    assert res["score"] == pytest.approx(mean * res["gate_factor"])


def test_toy_species_apart_from_reactions_cuts_every_edge(toy):
    nodes, edges, comps = toy
    ids = list(nodes["node_id"])
    l1 = ["RX" if i.startswith("r") else "SP" for i in ids]
    groups = _frame(sorted([[i, g, "", "", "R1", ""] for i, g in zip(ids, l1)]),
                    ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([["RX", "1", "", "reactions", "reaction_type", "none", "", ""],
                   ["SP", "1", "", "species", "reaction_type", "none", "", ""]], META_COLS)
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    v = dict(zip(res["table"].rule, res["table"].value))
    assert v["M09"] == pytest.approx(1.0)              # no edge has both ends in one box
    assert v["M04"] == pytest.approx(12.0 / 18.0)      # the species box holds two thirds of the map
    assert res["verdicts"]["M04"] == "fail"
    assert res["score"] <= SG.GATE_CAP + 1e-12         # gate cap applies


def test_toy_one_group_fails_the_graph_rules_rather_than_skipping_them(toy):
    nodes, edges, comps = toy
    ids = list(nodes["node_id"])
    groups = _frame(sorted([[i, "all", "", "", "R1", ""] for i in ids]),
                    ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([["all", "1", "", "everything", "mixed", "none", "", ""]], META_COLS)
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    for rid in SG.FAIL_IF_ONE_GROUP:
        assert res["verdicts"][rid] == "fail"
        assert res["utilities"][rid] == pytest.approx(0.0)
    assert res["verdicts"]["M04"] == "fail"
    assert res["score"] <= SG.GATE_CAP + 1e-12


def test_toy_singletons_fail_the_singleton_gate(toy):
    nodes, edges, comps = toy
    ids = list(nodes["node_id"])
    groups = _frame(sorted([[i, "g_" + i, "", "", "R1", ""] for i in ids]),
                    ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([["g_" + i, "1", "", i, "member", "member", i, ""] for i in sorted(ids)], META_COLS)
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    v = dict(zip(res["table"].rule, res["table"].value))
    assert v["M02"] == pytest.approx(1.0)              # every group holds one node
    assert res["verdicts"]["M02"] == "fail"
    assert res["score"] <= SG.GATE_CAP + 1e-12


def test_fan_metric_counts_child_boxes(toy):
    nodes, edges, comps = toy
    ids = list(nodes["node_id"])
    # one top group with three children: all parents are within the limit of twelve
    l2 = [f"k{(k % 3) + 1}" for k in range(len(ids))]
    groups = _frame(sorted([[i, "top", g, "", "R1", ""] for i, g in zip(ids, l2)]),
                    ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([["top", "1", "", "top", "mixed", "none", "", ""]]
                  + [[f"k{j}", "2", "top", f"k{j}", "topology", "none", "", ""] for j in (1, 2, 3)],
                  META_COLS)
    assert SG.fan_max12_share(nodes, edges, comps, groups, meta) == pytest.approx(1.0)
    # one top group with eighteen children: the only parent is over the limit
    l2b = ["k_" + i for i in ids]
    groups_b = _frame(sorted([[i, "top", g, "", "R1", ""] for i, g in zip(ids, l2b)]),
                      ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta_b = _frame([["top", "1", "", "top", "mixed", "none", "", ""]]
                    + [["k_" + i, "2", "top", i, "member", "member", i, ""] for i in sorted(ids)],
                    META_COLS)
    assert SG.fan_max12_share(nodes, edges, comps, groups_b, meta_b) == pytest.approx(0.0)


def test_single_child_rule_catches_a_dummy_level(toy):
    nodes, edges, comps = toy
    ids = list(nodes["node_id"])
    groups = _frame(sorted([[i, "P", "Q", "", "R1", ""] for i in ids]),
                    ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    meta = _frame([["P", "1", "", "P", "mixed", "none", "", ""],
                   ["Q", "2", "P", "Q", "mixed", "none", "", ""]], META_COLS)
    res = SG.score_from_tables(nodes, edges, comps, groups, meta)
    assert dict(zip(res["table"].rule, res["table"].value))["M06"] == pytest.approx(1.0)
    assert res["verdicts"]["M06"] == "fail"


def test_determinism_gate_zeroes_the_score(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    good = SG.score_from_tables(nodes, edges, comps, groups, meta, determinism=1.0)
    bad = SG.score_from_tables(nodes, edges, comps, groups, meta, determinism=0.0)
    assert good["verdicts"]["M27"] == "pass" and good["score"] > 0
    assert bad["verdicts"]["M27"] == "fail" and bad["score"] == 0.0


def test_scoring_is_deterministic(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    a = SG.score_from_tables(nodes, edges, comps, groups, meta)
    b = SG.score_from_tables(nodes, edges, comps, groups, meta)
    assert a["score"] == b["score"]
    pd.testing.assert_frame_equal(a["table"], b["table"])


def test_reference_rule_only_changes_block_e(toy):
    nodes, edges, comps = toy
    groups, meta = _blocks_grouping()
    ref = pd.DataFrame({"node_id": ["s01", "s02", "s03", "s05", "s06"],
                        "label": ["p1", "p1", "p1", "p2", "p2"],
                        "source": ["test"] * 5}).astype(str)
    with_ref = SG.score_from_tables(nodes, edges, comps, groups, meta, reference=ref)
    without = SG.score_from_tables(nodes, edges, comps, groups, meta)
    assert with_ref["verdicts"]["M25"] != "na"
    assert with_ref["score_core"] == pytest.approx(without["score"])
    for b in "ABCD":
        assert (isnan(with_ref["blocks"][b]) and isnan(without["blocks"][b])) or \
               with_ref["blocks"][b] == pytest.approx(without["blocks"][b])


# ============================================================ part 3: degenerate groupings on the panel
def _read(model, name):
    return SG.read_tsv(os.path.join(PANEL, model, name))


def _fair_meta(groups):
    """groups_meta for a synthetic grouping, built the way a careful strategy would declare it:
    a group whose id is `c_<space>` is a localisation group of that compartment, every other group is a
    topology group named after its first member (both are name sources the contract allows). Without
    this, `MG._meta` leaves every group unnamed and with basis `topology`, which would charge the
    compartment baselines for defects of the test fixture rather than of the grouping."""
    first_member, level_of, parent_of = {}, {}, {}
    for row in groups.to_dict("records"):
        up = ""
        for lv, col in enumerate(("L1", "L2", "L3"), start=1):
            g = str(row.get(col, "") or "")
            if not g:
                break
            level_of[g] = lv
            parent_of[g] = up
            first_member.setdefault(g, str(row["node_id"]))
            up = g
    rows = []
    for g in sorted(level_of):
        if g.startswith("c_") and level_of[g] == 1:
            cid = g[2:]
            rows.append([g, "1", parent_of[g], cid, "localisation", "compartment", cid, cid])
        else:
            rows.append([g, str(level_of[g]), parent_of[g], first_member[g], "topology", "member",
                         first_member[g], ""])
    return _frame(rows, META_COLS)


@pytest.mark.skipif(not os.path.isdir(PANEL), reason="panel/features not found")
@pytest.mark.parametrize("model", ["e_coli_core", "skm_ja"])
def test_no_degenerate_grouping_ranks_first(model):
    nodes, edges = _read(model, "nodes.tsv"), _read(model, "edges.tsv")
    comps = _read(model, "compartments.tsv")
    syn = MG.synthetic(nodes, edges)
    cases = {}
    for name in ("one_group", "singletons", "compartment", "random_sized",
                 "optimiser_nested", "compartment_optimiser"):
        g = syn[name]
        meta = _fair_meta(g)
        cases[name] = SG.score_from_tables(nodes, edges, comps, g, meta)["score"]
    order = sorted(cases, key=lambda k: -cases[k])
    assert order[0] not in ("one_group", "singletons", "compartment", "random_sized"), cases
    for degenerate in ("one_group", "singletons", "random_sized"):
        assert cases[degenerate] < cases[order[0]], cases


@pytest.mark.skipif(not os.path.isdir(PANEL), reason="panel/features not found")
@pytest.mark.parametrize("model", ["e_coli_core", "skm_ja"])
def test_degenerate_groupings_are_capped_by_the_gates(model):
    nodes, edges = _read(model, "nodes.tsv"), _read(model, "edges.tsv")
    comps = _read(model, "compartments.tsv")
    syn = MG.synthetic(nodes, edges)
    for name in ("one_group", "singletons"):
        res = SG.score_from_tables(nodes, edges, comps, syn[name], _fair_meta(syn[name]))
        assert res["score"] <= SG.GATE_CAP + 1e-12, (name, res["score"])


@pytest.mark.skipif(not os.path.isdir(PANEL), reason="panel/features not found")
def test_one_group_per_compartment_is_not_degenerate_but_does_not_win():
    """One box per compartment is what every strategy gets for free: it must score above the degenerate
    groupings and below a grouping that also splits by modules."""
    model = "e_coli_core"
    nodes, edges = _read(model, "nodes.tsv"), _read(model, "edges.tsv")
    comps = _read(model, "compartments.tsv")
    syn = MG.synthetic(nodes, edges)
    comp = SG.score_from_tables(nodes, edges, comps, syn["compartment"],
                                _fair_meta(syn["compartment"]))["score"]
    one = SG.score_from_tables(nodes, edges, comps, syn["one_group"], _fair_meta(syn["one_group"]))["score"]
    both = SG.score_from_tables(nodes, edges, comps, syn["compartment_optimiser"],
                                _fair_meta(syn["compartment_optimiser"]))["score"]
    assert one < comp < both, (one, comp, both)


@pytest.mark.skipif(not os.path.isdir(PANEL), reason="panel/features not found")
def test_panel_rows_have_no_zero_standing_in_for_na():
    """On a real model every rule is either applicable with a value or NaN; a not-applicable rule never
    reports 0."""
    model = "skm_ja"
    nodes, edges = _read(model, "nodes.tsv"), _read(model, "edges.tsv")
    comps = _read(model, "compartments.tsv")
    g = MG.synthetic(nodes, edges)["optimiser"]
    res = SG.score_from_tables(nodes, edges, comps, g, _fair_meta(g))
    tbl = res["table"]
    for _, row in tbl.iterrows():
        if res["verdicts"][row.rule] == "na":
            assert isnan(row.value) and isnan(row.utility)


@pytest.mark.skipif(not os.path.isdir(os.path.join(HERE, "groupings")), reason="strategies not found")
def test_cli_single_grouping_runs_and_writes_json(tmp_path):
    out = tmp_path / "one.json"
    rc = SG.main(["--features", os.path.join(PANEL, "skm_ja"),
                  "--grouping", os.path.join(HERE, "groupings", "hybrid", "skm_ja"),
                  "--reference", os.path.join(_ROOT, "heldout", "skm_ja.reference.tsv"),
                  "--determinism", "1", "--json", str(out)])
    assert rc == 0 and out.is_file()
    import json
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert len(payload["rules"]) == 27
    assert 0.0 <= payload["score"] <= 1.0
    assert payload["verdicts"]["M27"] == "pass"
