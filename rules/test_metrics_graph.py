r"""Tests for metrics_graph.py: small synthetic cases with values known by construction, and the
degenerate groupings on the panel models.

Run (from the folder that holds metrics_graph.py):
    python -m pytest test_metrics_graph.py -q
    # the panel tests look for panel/features next to, above, or in ./bundle of this folder; or say where:
    $env:GROUPING_PANEL = "C:\path\to\bundle\panel\features"; python -m pytest test_metrics_graph.py -q     # PowerShell
    GROUPING_PANEL=/path/to/bundle/panel/features python -m pytest test_metrics_graph.py -q                # bash
Panel tests are skipped (not failed) when the panel is not found.
"""
import math
import os
import pathlib

import pandas as pd
import pytest

import metrics_graph as mg

HERE = pathlib.Path(__file__).resolve().parent
PANEL_MODELS = ["brenda_merged", "e_coli_core", "R-HSA-77289", "skm_pss", "skm_ja", "skm_sa"]


# ------------------------------------------------------------------------------------------ builders
def model(reactions: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """reactions: {reaction_id: [species ids]} -> nodes, edges, comps (all strings, like the contract)."""
    species = sorted({s for part in reactions.values() for s in part})
    nodes = pd.DataFrame({"node_id": species + sorted(reactions),
                          "kind": ["species"] * len(species) + ["reaction"] * len(reactions)})
    nodes["loc_space"] = "c"
    nodes["sbo"] = ""
    edges = pd.DataFrame([(r, s, "reactant", "1.0", "", "") for r, part in sorted(reactions.items()) for s in part],
                         columns=["reaction_id", "species_id", "role", "stoichiometry", "sbo", "sbo_name"])
    comps = pd.DataFrame({"compartment_id": ["c"], "parent_id": [""]})
    return nodes.astype(str), edges.astype(str), comps.astype(str)


def grouping(nodes, l1: dict, l2: dict | None = None, l3: dict | None = None, clone: dict | None = None):
    ids = list(nodes["node_id"])
    g = pd.DataFrame({"node_id": ids, "L1": [l1[i] for i in ids],
                      "L2": [(l2 or {}).get(i, "") for i in ids], "L3": [(l3 or {}).get(i, "") for i in ids],
                      "rule": "R1", "clone_in": [(clone or {}).get(i, "") for i in ids]}).astype(str)
    return g, mg._meta(g)


def run(reactions, l1, l2=None, l3=None, clone=None, reference=None):
    nodes, edges, comps = model(reactions)
    g, meta = grouping(nodes, l1, l2, l3, clone)
    return mg.metrics(nodes, edges, comps, g, meta, reference)


def by_group(**groups) -> dict:
    return {n: g for g, members in groups.items() for n in members.split()}


TWO_PATHS = {"r1": ["a", "b"], "r2": ["c", "d"]}
CHAIN = {"r1": ["a", "b"], "r2": ["b", "c"]}


def star(n=20, extra=True):
    """Hub h in n reactions, each with two private species; optionally one reaction without the hub."""
    rx = {f"r{i:02d}": ["h", f"p{i:02d}", f"q{i:02d}"] for i in range(1, n + 1)}
    if extra:
        rx["rzz"] = ["x", "y"]
    return rx


def star_groups(n=20, extra=True):
    l1 = {}
    for i in range(1, n + 1):
        for v in (f"r{i:02d}", f"p{i:02d}", f"q{i:02d}"):
            l1[v] = f"g{i:02d}"
    l1["h"] = "g01"
    if extra:
        l1.update({"rzz": "gzz", "x": "gzz", "y": "gzz"})
    return l1


# --------------------------------------------------------------------------- values by construction
def test_two_separate_paths_perfect_grouping():
    m = run(TWO_PATHS, by_group(A="r1 a b", B="r2 c d"))
    assert m["n_nodes"] == 6 and m["n_edges"] == 4 and m["depth"] == 1
    assert m["cut_w_leaf"] == 0 and m["cut_raw_leaf"] == 0
    assert m["Q_w_L1"] == pytest.approx(0.5) and m["Q_raw_L1"] == pytest.approx(0.5)
    assert m["Q_ref"] == pytest.approx(0.5) and m["Q_eff"] == pytest.approx(1.0)
    assert m["connected_node_share_leaf"] == 1 and m["disconnected_group_share_leaf"] == 0
    assert m["weak_node_share_leaf"] == 1 and m["conductance_mean_leaf"] == 0
    assert m["largest_share_leaf"] == pytest.approx(0.5) and m["singleton_group_share_leaf"] == 0
    assert m["n_eff_groups_leaf"] == pytest.approx(2.0)
    # every edge sits in a box of 3 of 6 nodes; random filling: 4 of 8 pairs inside (0.5), 4 outside (1)
    assert m["edge_span_w"] == pytest.approx(0.5)
    assert m["edge_span_null"] == pytest.approx(0.75)
    assert m["edge_span_ratio"] == pytest.approx(2 / 3)
    assert math.isnan(m["home_plurality_share_leaf"])          # no species with two reactions
    assert math.isnan(m["Q_w_L2"]) and math.isnan(m["cut_w_L3"])  # unused levels
    assert math.isnan(m["clone_misuse_share"]) and m["n_clones"] == 0
    assert math.isnan(m["ref_bcubed_f1_best"])


def test_chain_one_crossing_edge_weights_and_modularity():
    # a - r1 - b - r2 - c ; b has two reactions, so its edges weigh 1/2; total weight W = 3 species
    m = run(CHAIN, by_group(A="a r1 b", B="r2 c"))
    assert m["cut_w_leaf"] == pytest.approx(0.5 / 3)
    assert m["cut_raw_leaf"] == pytest.approx(1 / 4)
    assert m["Q_w_L1"] == pytest.approx((1.5 / 3 - 1.5 * 2 / 9) + (1 / 3 - 1.5 * 1 / 9))     # = 1/3
    assert m["Q_raw_L1"] == pytest.approx((2 / 4 - 2 * 3 / 16) + (1 / 4 - 2 * 1 / 16))       # = 1/4
    assert m["conductance_mean_leaf"] == pytest.approx(0.2)     # both groups: cut 0.5 / min(vol, 2W - vol) 2.5
    assert m["weak_node_share_leaf"] == 1
    assert m["home_plurality_share_leaf"] == 1                  # b: one reaction at home, one away (tie)
    assert m["edge_span_w"] == pytest.approx((1 * 3 + 0.5 * 3 + 0.5 * 5 + 1 * 2) / (5 * 3))  # = 0.6
    assert m["edge_span_null"] == pytest.approx((2 * 3 / 5 + 1 * 2 / 5 + 3 * 1) / 6)
    assert m["hub_degree_cutoff"] >= 2 and m["n_hubs"] == 0


def test_species_in_the_wrong_box_fails_home_plurality():
    # b takes part in r1, r2, r3; it sits with r1 only, r2 and r3 share another box
    m = run({"r1": ["a", "b"], "r2": ["b", "c"], "r3": ["b", "d"]}, by_group(A="a r1 b", B="r2 c r3 d"))
    assert m["home_plurality_share_leaf"] == 0
    assert m["cut_w_leaf"] == pytest.approx((2 / 3) / 4)        # two of b's three edges cross, W = 4


def test_one_group_and_singletons_on_a_small_model():
    nodes, _, _ = model(CHAIN)
    one = run(CHAIN, {n: "all" for n in nodes["node_id"]})
    assert one["n_groups_leaf"] == 1 and one["largest_share_leaf"] == 1
    for key in ("cut_w_leaf", "Q_w_L1", "Q_w_best", "Q_eff", "weak_node_share_leaf", "conductance_mean_leaf"):
        assert math.isnan(one[key]), key
    assert one["edge_span_w"] == pytest.approx(1.0) and one["edge_span_ratio"] == pytest.approx(1.0)
    assert one["connected_node_share_leaf"] == 1
    single = run(CHAIN, {n: f"s_{n}" for n in nodes["node_id"]})
    assert single["cut_w_leaf"] == 1 and single["cut_raw_leaf"] == 1
    assert single["singleton_group_share_leaf"] == 1 and single["singleton_node_share_leaf"] == 1
    assert single["Q_w_L1"] == pytest.approx(0.0) and single["Q_eff"] == pytest.approx(0.0)
    assert single["weak_node_share_leaf"] == 0 and single["conductance_mean_leaf"] == pytest.approx(1.0)
    assert single["edge_span_w"] == pytest.approx(1.0) and single["edge_span_ratio"] == pytest.approx(1.0)
    assert math.isnan(single["connected_node_share_leaf"])      # no group with two nodes
    assert single["readable_node_share_leaf"] == 0
    v1, v2 = mg.judge(one), mg.judge(single)
    assert v1["G1"] == v1["G2"] == v1["G5"] == v1["G6"] == v1["G7"] == v1["G9"] == "fail"
    assert v2["G1"] == v2["G3"] == v2["G5"] == v2["G6"] == v2["G7"] == "fail" and v2["G4"] == "na"
    assert mg.score(one) <= mg.GATE_CAP and mg.score(single) <= mg.GATE_CAP


def test_disconnected_group_is_counted():
    rx = {"r1": ["a", "b"], "r2": ["c", "d"], "r3": ["e", "f"]}
    m = run(rx, by_group(A="r1 a b r2 c d", B="r3 e f"))
    assert m["connected_node_share_leaf"] == pytest.approx((3 + 3) / 9)
    assert m["disconnected_group_share_leaf"] == pytest.approx(1 / 2)
    assert m["cut_w_leaf"] == 0                                 # a bag of two pathways still cuts no edge


def test_groups_of_one_kind_are_not_communities():
    m = run(TWO_PATHS, by_group(S="a b c d", R="r1 r2"))
    assert m["cut_w_leaf"] == 1 and m["weak_node_share_leaf"] == 0
    assert m["connected_node_share_leaf"] == pytest.approx(2 / 6)
    assert m["Q_w_L1"] == pytest.approx(0.0)                   # Barber: no reaction-species pair inside a group


def test_hub_cutoff_is_data_derived():
    # 42 species with one reaction and one with 20: lambda of the non-hubs is 1, alpha = 0.05 / 43
    # P(Poisson(1) > 4) = 0.0037 > alpha >= P(Poisson(1) > 5) = 0.0006  ->  cutoff 5
    m = run(star(), star_groups())
    assert m["hub_degree_cutoff"] == 5 and m["n_hubs"] == 1
    assert m["hub_edge_share_raw"] == pytest.approx(20 / 62)
    assert mg.hub_cutoff([1] * 42 + [20]) == 5
    assert mg.hub_cutoff([3, 3, 3, 3]) >= 3                     # equal degrees: nobody is a hub


def test_hub_dominates_raw_cut_but_not_weighted_cut():
    m = run(star(), star_groups())
    assert m["cut_raw_leaf"] == pytest.approx(19 / 62)          # 19 of the hub's 20 edges cross
    assert m["cut_w_leaf"] == pytest.approx((19 / 20) / 43)     # the hub has one vote among 43 species
    assert m["connected_node_share_leaf"] == 1 and m["connected_nohub_node_share_leaf"] == 1


def test_clones_hub_copies_count_other_copies_do_not():
    l1 = star_groups()
    hub_clones = "|".join(f"g{i:02d}" for i in range(2, 21))
    m = run(star(), l1, clone={"h": hub_clones})
    assert m["n_clones"] == 19 and m["cut_w_leaf"] == 0 and m["cut_raw_leaf"] == 0
    assert m["cut_w_home_leaf"] == pytest.approx((19 / 20) / 43)
    assert m["clone_absorbed_w"] == pytest.approx((19 / 20) / 43)
    assert m["clone_idle_share"] == 0 and m["clone_nonhub_share"] == 0 and m["clone_misuse_share"] == 0
    assert mg.judge(m)["G11"] == "pass"
    # the modularity ignores clones altogether
    assert m["Q_w_L1"] == pytest.approx(run(star(), l1)["Q_w_L1"])
    # an idle copy of the hub (gzz has no reaction of h) and a copy of a non-hub species
    m2 = run(star(), l1, clone={"h": "gzz", "p01": "g02"})
    assert m2["n_clones"] == 2
    assert m2["clone_idle_share"] == 1.0                        # p01 has no reaction in g02 either
    assert m2["clone_nonhub_share"] == pytest.approx(0.5) and m2["clone_misuse_share"] == 1.0
    assert m2["cut_w_leaf"] == pytest.approx(m2["cut_w_home_leaf"])
    assert mg.judge(m2)["G11"] == "fail"


def test_non_hub_clone_never_lowers_the_cut():
    # b is shared by r1 and r2 and is no hub: a copy next to r2 is drawn but the edge still counts
    plain = run(CHAIN, by_group(A="a r1 b", B="r2 c"))
    cloned = run(CHAIN, by_group(A="a r1 b", B="r2 c"), clone={"b": "B"})
    assert cloned["n_clones"] == 1 and cloned["clone_nonhub_share"] == 1
    assert cloned["cut_w_leaf"] == pytest.approx(plain["cut_w_leaf"])
    assert cloned["edge_span_w"] == pytest.approx(plain["edge_span_w"])


NESTED = {"r1": ["a1", "b1", "x"], "r2": ["a2", "b2", "x", "y"], "r3": ["a3", "b3", "y"], "r4": ["a4", "b4"]}
NESTED_L1 = by_group(A="r1 a1 b1 x r2 a2 b2 y", B="r3 a3 b3 r4 a4 b4")
NESTED_L2 = by_group(P1="r1 a1 b1 x", P2="r2 a2 b2 y", P3="r3 a3 b3", P4="r4 a4 b4")


def test_two_levels_cut_grows_with_depth_and_span_integrates_levels():
    m = run(NESTED, NESTED_L1, NESTED_L2)
    assert m["depth"] == 2 and m["n_groups_L1"] == 2 and m["n_groups_L2"] == 4 and m["n_groups_leaf"] == 4
    assert m["cut_w_L1"] == pytest.approx(0.5 / 10)             # (r3, y)
    assert m["cut_w_L2"] == pytest.approx(1.0 / 10)             # (r3, y) and (r2, x)
    assert m["cut_w_leaf"] == m["cut_w_L2"] and m["cut_w_L1"] <= m["cut_w_L2"]
    assert m["cut_w_mean"] == pytest.approx(0.075)
    assert math.isnan(m["cut_w_L3"])
    assert m["Q_w_best"] == pytest.approx(max(m["Q_w_L1"], m["Q_w_L2"]))
    # boxes: P1 4, P2 4, P3 3, P4 3, A 8, whole model 14
    assert m["edge_span_w"] == pytest.approx((2.5 * 4 + 2.5 * 4 + 0.5 * 8 + 2 * 3 + 0.5 * 14 + 2 * 3) / (14 * 10))
    assert m["max_siblings"] == 2 and m["max_children"] == 2
    assert m["single_child_share"] == 0 and m["residual_node_share"] == 0
    flat = run(NESTED, NESTED_L2)
    assert flat["cut_w_leaf"] == pytest.approx(m["cut_w_leaf"])
    assert m["edge_span_w"] < flat["edge_span_w"]               # the nesting keeps (r2, x) inside box A


def test_node_placed_directly_in_a_parent_forms_a_residual_group():
    l2 = dict(NESTED_L2)
    l2["x"] = ""
    m = run(NESTED, NESTED_L1, l2)
    assert m["n_groups_L2"] == 5 and m["residual_node_share"] == pytest.approx(1 / 14)
    assert m["cut_w_L2"] == pytest.approx(1.5 / 10)             # both edges of x now leave a sub-box
    assert m["cut_w_L1"] == pytest.approx(0.5 / 10)
    assert m["max_children"] == 3                               # P1, P2 and the loose node


def test_redundant_nesting_is_flagged():
    nodes, _, _ = model(NESTED)
    copy = {n: "sub_" + g for n, g in NESTED_L2.items()}
    m = run(NESTED, NESTED_L2, copy)                            # every L1 group has exactly one child
    assert m["single_child_share"] == 1.0 and mg.judge(m)["G10"] == "fail"
    dummy = run(NESTED, {n: "root" for n in nodes["node_id"]}, NESTED_L2)
    assert dummy["single_child_share"] == pytest.approx(0.5)    # the root has one child: the dummy
    assert math.isnan(dummy["cut_w_L1"]) and dummy["cut_w_leaf"] == pytest.approx(0.1)
    v = mg.judge(dummy)
    assert v["G9"] == "fail" and v["G10"] == "fail"
    assert mg.score(dummy) < mg.score(run(NESTED, NESTED_L2))


def test_three_levels():
    l3 = {"r1": "P1a", "a1": "P1a", "b1": "P1b", "x": "P1b"}
    m = run(NESTED, NESTED_L1, NESTED_L2, l3)
    assert m["depth"] == 3 and m["n_groups_L3"] == 5
    assert m["cut_w_L1"] <= m["cut_w_L2"] <= m["cut_w_L3"]
    assert m["cut_w_leaf"] == m["cut_w_L3"]
    assert m["cut_w_L3"] == pytest.approx((1.0 + 1.0 + 0.5) / 10)   # + (r1, b1) and (r1, x)


def test_reference_bcubed_and_reference_cut():
    ref = pd.DataFrame({"node_id": ["a", "b", "c", "d"], "label": ["X", "X", "Y", "Y"]})
    good = run(TWO_PATHS, by_group(A="r1 a b", B="r2 c d"), reference=ref)
    assert good["ref_n_labelled"] == 4 and good["ref_coverage"] == pytest.approx(4 / 6)
    assert good["ref_bcubed_precision_leaf"] == 1 and good["ref_bcubed_recall_leaf"] == 1
    assert good["ref_bcubed_f1_best"] == 1
    assert good["ref_cut_w"] == 0                               # reactions inherit X and Y: no edge crosses
    nodes, _, _ = model(TWO_PATHS)
    one = run(TWO_PATHS, {n: "all" for n in nodes["node_id"]}, reference=ref)
    assert one["ref_bcubed_precision_leaf"] == pytest.approx(0.5) and one["ref_bcubed_recall_leaf"] == 1
    assert one["ref_bcubed_f1_leaf"] == pytest.approx(2 / 3)
    single = run(TWO_PATHS, {n: f"s_{n}" for n in nodes["node_id"]}, reference=ref)
    assert single["ref_bcubed_precision_leaf"] == 1 and single["ref_bcubed_recall_leaf"] == pytest.approx(0.5)
    # a node with two labels agrees with both of its neighbours
    multi = pd.DataFrame({"node_id": ["a", "b", "b", "c"], "label": ["X", "X", "Y", "Y"]})
    mm = run(CHAIN, by_group(A="a r1 b", B="r2 c"), reference=multi)
    assert mm["ref_bcubed_precision_leaf"] == 1
    assert mm["ref_bcubed_recall_leaf"] == pytest.approx((1 + 2 / 3 + 1 / 2) / 3)
    # a reference that splits a reaction from one of its species has a cut of its own
    bad = pd.DataFrame({"node_id": ["a", "b", "c"], "label": ["X", "Y", "Y"]})
    mb = run(CHAIN, by_group(A="a r1 b", B="r2 c"), reference=bad)
    assert mb["ref_cut_w"] == pytest.approx(0.5 / 3)            # r1 takes X (tie, sorted first): (r1, b) crosses


def test_judge_thresholds_gates_and_not_applicable():
    base = {"n_groups_leaf": 10.0, "n_groups_L1": 5.0}
    m = dict(base, cut_w_leaf=0.25, largest_share_leaf=0.2, singleton_group_share_leaf=0.0,
             connected_node_share_leaf=1.0, Q_eff=0.9, weak_node_share_leaf=0.9, edge_span_ratio=0.2,
             readable_node_share_leaf=0.9, max_siblings=9.0, single_child_share=float("nan"),
             clone_misuse_share=float("nan"), ref_bcubed_f1_best=float("nan"))
    v = mg.judge(m)
    assert all(v[r] == "pass" for r in ("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"))
    assert v["G10"] == v["G11"] == v["G12"] == "na"
    assert mg.score(m) == 1.0
    warn = dict(m, cut_w_leaf=0.4)
    assert mg.judge(warn)["G1"] == "warn" and mg.score(warn) == pytest.approx(1 - 0.5 * 3 / 16)
    assert mg.judge(dict(m, cut_w_leaf=0.51))["G1"] == "fail"
    gated = dict(m, singleton_group_share_leaf=0.48)            # the current pipeline's SKM figure
    assert mg.judge(gated)["G3"] == "fail" and mg.score(gated) == mg.GATE_CAP
    assert mg.judge(dict(m, largest_share_leaf=1.0))["G2"] == "fail"
    assert mg.judge(dict(m, max_siblings=10.0))["G9"] == "warn"
    assert mg.judge(dict(m, max_siblings=51.0))["G9"] == "fail"


def test_input_checks_and_determinism():
    nodes, edges, comps = model(CHAIN)
    g, meta = grouping(nodes, by_group(A="a r1 b", B="r2 c"))
    with pytest.raises(ValueError):
        mg.metrics(nodes, edges, comps, g.iloc[:-1], meta)
    with pytest.raises(ValueError):
        mg.metrics(nodes, edges, comps, pd.concat([g, g.iloc[:1]]), meta)
    a, b = mg.metrics(nodes, edges, comps, g, meta), mg.metrics(nodes, edges, comps, g.iloc[::-1], meta)
    assert a.keys() == b.keys()
    for k in a:
        assert (math.isnan(a[k]) and math.isnan(b[k])) or a[k] == b[k], k
        assert isinstance(a[k], float), k
    # a duplicated participant row (two roles of one species in one reaction) is one edge
    dup = pd.concat([edges, edges.iloc[:1].assign(role="modifier")])
    assert mg.metrics(nodes, dup, comps, g, meta)["n_edges"] == a["n_edges"]


# ------------------------------------------------------------------------------------- panel models
def _panel():
    cands = [os.environ.get("GROUPING_PANEL", "")] + [str(p / q) for p in (HERE, HERE.parent, HERE.parent.parent)
                                                      for q in ("panel/features", "bundle/panel/features")]
    for c in cands:
        if c and (pathlib.Path(c) / "skm_ja" / "nodes.tsv").exists():
            return pathlib.Path(c)
    return None


PANEL = _panel()
_cache = {}


def panel_scores(name):
    if name not in _cache:
        f = PANEL / name
        nodes, edges, comps = (mg._read(f / t) for t in ("nodes.tsv", "edges.tsv", "compartments.tsv"))
        out = {}
        for label, grp in mg.synthetic(nodes, edges).items():
            m = mg.metrics(nodes, edges, comps, grp, mg._meta(grp))
            out[label] = (m, mg.judge(m), mg.score(m))
        _cache[name] = out
    return _cache[name]


needs_panel = pytest.mark.skipif(PANEL is None, reason="panel/features not found; set GROUPING_PANEL")


@needs_panel
@pytest.mark.parametrize("name", PANEL_MODELS)
def test_panel_degenerate_groupings_never_win(name):
    s = panel_scores(name)
    best = max(v[2] for v in s.values())
    for bad in ("one_group", "singletons", "compartment", "random_sized"):
        assert s[bad][2] < best, bad
    for bad in ("one_group", "singletons", "random_sized", "random_clone_all", "species_apart_from_reactions",
                "per_reaction"):
        assert s[bad][2] <= mg.GATE_CAP + 1e-9, bad
        for good in ("optimiser", "optimiser_nested", "compartment_optimiser"):
            assert s[bad][2] < s[good][2], (bad, good)
    assert s["compartment"][2] < min(s["optimiser"][2], s["optimiser_nested"][2])


@needs_panel
@pytest.mark.parametrize("name", PANEL_MODELS)
def test_panel_each_degenerate_fails_the_rule_meant_for_it(name):
    s = panel_scores(name)
    one, single, rand = s["one_group"], s["singletons"], s["random_sized"]
    assert one[0]["n_groups_leaf"] == 1 and math.isnan(one[0]["cut_w_leaf"])
    assert [one[1][r] for r in ("G1", "G2", "G5", "G6", "G7", "G9")] == ["fail"] * 6
    assert single[0]["cut_w_leaf"] == 1.0 and single[0]["singleton_group_share_leaf"] == 1.0
    assert [single[1][r] for r in ("G1", "G3", "G5", "G6", "G7", "G8", "G9")] == ["fail"] * 7
    assert one[0]["edge_span_ratio"] == pytest.approx(1.0) and single[0]["edge_span_ratio"] == pytest.approx(1.0)
    assert abs(rand[0]["Q_eff"]) < 0.1 and rand[0]["edge_span_ratio"] > 0.9 and rand[0]["cut_w_leaf"] > 0.8
    assert [rand[1][r] for r in ("G1", "G4", "G5", "G6", "G7")] == ["fail"] * 5
    apart = s["species_apart_from_reactions"]
    assert apart[0]["cut_w_leaf"] == 1.0 and apart[1]["G4"] == "fail"


@needs_panel
@pytest.mark.parametrize("name", PANEL_MODELS)
def test_panel_gaming_is_caught(name):
    s = panel_scores(name)
    # cloning every species into every group it touches: the cut barely moves, the clone rule fails
    gamed, rand = s["random_clone_all"], s["random_sized"]
    assert gamed[1]["G11"] == "fail" and gamed[1]["G1"] == "fail"
    assert gamed[0]["cut_w_home_leaf"] == pytest.approx(rand[0]["cut_w_leaf"])
    assert gamed[2] <= rand[2] + 1e-9
    # a dummy top level over a good partition scores below the same partition without it
    assert s["dummy_root_optimiser"][2] < s["optimiser"][2]
    assert s["dummy_root_optimiser"][1]["G10"] == "fail"
    # hub copies help the optimiser's partition and are not misuse
    hubs = s["optimiser_clone_hubs"]
    assert hubs[0]["cut_w_leaf"] <= s["optimiser"][0]["cut_w_leaf"] + 1e-12
    assert hubs[0]["clone_misuse_share"] == 0 and hubs[1]["G11"] == "pass"
    assert hubs[0]["Q_w_best"] == pytest.approx(s["optimiser"][0]["Q_w_best"])


@needs_panel
@pytest.mark.parametrize("name", PANEL_MODELS)
def test_panel_level_behaviour(name):
    s = panel_scores(name)
    for label in ("optimiser_nested", "compartment_optimiser"):
        m = s[label][0]
        if m["depth"] == 2 and not math.isnan(m["cut_w_L1"]):
            assert m["cut_w_L1"] <= m["cut_w_L2"] + 1e-12
            assert m["largest_share_L1"] >= m["largest_share_L2"] - 1e-12
    nested, flat = s["optimiser_nested"][0], s["optimiser"][0]
    assert nested["edge_span_w"] <= flat["edge_span_w"] + 1e-12     # same leaves, nested boxes: never worse
    assert nested["max_siblings"] <= flat["max_siblings"]
    assert 0 < flat["Q_ref"] <= 1 and flat["Q_eff"] == pytest.approx(1.0)


@needs_panel
def test_panel_hubs_are_the_expected_kind_of_species():
    f = PANEL / "e_coli_core"
    nodes, edges = mg._read(f / "nodes.tsv"), mg._read(f / "edges.tsv")
    g = mg._Graph(nodes, edges)
    hubs = set(nodes["node_id"][g.is_hub])
    assert {"M_h_c", "M_h2o_c", "M_h_e", "M_atp_c"} == hubs       # protons, water, ATP: found from degrees alone
    assert g.hub_cut == 12
