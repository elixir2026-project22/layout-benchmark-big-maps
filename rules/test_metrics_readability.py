"""Tests of metrics_readability.py.

Run:
    python -m pytest test_metrics_readability.py -q
    python -m pytest test_metrics_readability.py -q -k panel          (only the panel-model tests)

Part 1: small synthetic cases whose correct values are known by construction.
Part 2: degenerate and baseline groupings (everything in one group, every node its own group, one group
per compartment, a random grouping with realistic sizes, a Louvain grouping) on all six panel models.
The panel features are looked up in $PANEL_FEATURES, then in ./panel/features, ./bundle/panel/features
and the same two one and two folders up; part 2 is skipped when none exists.
"""
import math
import os
import pathlib

import numpy as np
import pandas as pd
import pytest

import metrics_readability as mr

GC = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
MC = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]
NC = ["node_id", "kind", "name", "compartment", "loc_space", "loc_kind", "sbo_path", "ec", "go", "xrefs"]


def df(rows, cols):
    return pd.DataFrame([[str(x) for x in r] for r in rows], columns=cols)


def nodes_df(rows):
    """rows: (node_id, kind, loc_space, loc_kind[, ec[, sbo_path[, go[, xrefs]]]])"""
    out = []
    for r in rows:
        r = list(r) + [""] * (8 - len(r))
        out.append([r[0], r[1], r[0], r[2], r[2], r[3], r[5], r[4], r[6], r[7]])
    return df(out, NC)


def meta_row(gid, level=1, parent="", label="", basis="mixed", source=None, xref="", comp=""):
    return [gid, level, parent, label, basis, source or ("none" if not label else "notes"), xref, comp]


COMPS = df([["e", "extracellular", ""], ["c", "cytosol", "e"]], ["compartment_id", "name", "parent_id"])


# ------------------------------------------------------------------------------------------ toy model
# A_e --T1--> A_c --R1--> B_c --R2--> C_c --R3--> D_c ; X_c modifies R1 and R2. 10 nodes, 10 edges.
def toy():
    nodes = nodes_df([("A_e", "species", "e", "inside"), ("A_c", "species", "c", "inside"),
                      ("B_c", "species", "c", "inside"), ("C_c", "species", "c", "inside"),
                      ("D_c", "species", "c", "inside"), ("X_c", "species", "c", "inside"),
                      ("T1", "reaction", "c", "boundary"), ("R1", "reaction", "c", "inside"),
                      ("R2", "reaction", "c", "inside"), ("R3", "reaction", "c", "inside")])
    edges = df([("T1", "A_e", "reactant"), ("T1", "A_c", "product"),
                ("R1", "A_c", "reactant"), ("R1", "B_c", "product"), ("R1", "X_c", "modifier"),
                ("R2", "B_c", "reactant"), ("R2", "C_c", "product"), ("R2", "X_c", "modifier"),
                ("R3", "C_c", "reactant"), ("R3", "D_c", "product")], ["reaction_id", "species_id", "role"])
    return nodes, edges


def toy_groups(clone=None):
    """P = {A_e, T1, A_c, R1, B_c}, Q = {R2, C_c, R3, D_c, X_c}."""
    clone = clone or {}
    part = {"A_e": "P", "T1": "P", "A_c": "P", "R1": "P", "B_c": "P",
            "R2": "Q", "C_c": "Q", "R3": "Q", "D_c": "Q", "X_c": "Q"}
    g = df([[n, p, "", "", "S", clone.get(n, "")] for n, p in sorted(part.items())], GC)
    m = df([meta_row("P"), meta_row("Q")], MC)
    return g, m


def test_toy_edges_balance_and_compartment_pieces():
    nodes, edges = toy()
    g, m = toy_groups()
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["n_nodes"] == 10 and r["n_groups"] == 2 and r["placed_share"] == 1.0
    assert r["top_boxes"] == 2
    assert math.isnan(r["top_groups"])                       # 10 nodes: the rule does not apply
    # R1-X_c and R2-B_c leave their L1 group: 2 of 10 edges
    assert r["long_edge_share"] == pytest.approx(0.2)
    assert r["long_edge_share_raw"] == pytest.approx(0.2)
    assert r["cut_edge_share"] == pytest.approx(0.2)
    assert r["boundaries_crossed_mean"] == pytest.approx(0.4)   # each long edge crosses two L1 borders
    # 4 edges inside each group, degree sum 10 each: Q = 2 * (4/10 - (10/20)^2) = 0.3
    assert r["modularity_L1"] == pytest.approx(0.3)
    assert r["size_balance"] == pytest.approx(1.0)            # 5 + 5
    assert r["largest_box_share"] == pytest.approx(0.5)
    assert r["singleton_share"] == 0.0
    assert r["oversize_node_share"] == 0.0
    # P mixes e and c, Q is pure c: e is drawn once (in P), c twice (in P and on the map) -> 1.5
    assert r["comp_redraw"] == pytest.approx(1.5) and r["comp_redraw_max"] == 2
    assert r["mixed_box_node_share"] == pytest.approx(0.5)
    assert r["named_share"] == 0.0 and r["distinct_share"] == 0.0 and r["backed_share"] == 0.0
    assert math.isnan(r["fan_ok_share"])                      # flat grouping
    assert r["depth_mean"] == 1.0 and r["depth_excess"] == 0.0
    # clones needed to make every edge local: (X_c in P), (B_c in Q) -> 2 per 6 species
    assert r["clone_demand"] == pytest.approx(2 / 6)
    assert r["clone_load"] == 0.0 and r["cloned_species_share"] == 0.0
    assert math.isnan(r["clone_useless_share"])


def test_toy_clones():
    nodes, edges = toy()
    g, m = toy_groups(clone={"X_c": "P"})                     # a useful clone: R1 is in P
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["long_edge_share"] == pytest.approx(0.1)         # R1-X_c is now local
    assert r["long_edge_share_raw"] == pytest.approx(0.2)     # the grouping itself did not change
    assert r["clone_load"] == pytest.approx(1 / 6)
    assert r["cloned_species_share"] == pytest.approx(1 / 6)
    assert r["clone_useless_share"] == 0.0
    g, m = toy_groups(clone={"X_c": "P", "D_c": "P"})         # D_c has no reaction in P
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["clone_useless_share"] == pytest.approx(0.5)
    assert r["cloned_species_share"] == pytest.approx(2 / 6)
    assert r["long_edge_share"] == pytest.approx(0.1)


def test_toy_degenerate_groupings():
    nodes, edges = toy()
    one = (df([[n, "ALL", "", "", "S", ""] for n in nodes.node_id], GC), df([meta_row("ALL")], MC))
    r = mr.metrics(nodes, edges, COMPS, *one)
    assert r["long_edge_share"] == 0.0 and r["clone_demand"] == 0.0
    assert r["comp_redraw"] == 1.0                            # each compartment appears in one place
    assert math.isnan(r["size_balance"]) and r["largest_box_share"] == 1.0
    single = (df([[n, "g" + n, "", "", "S", ""] for n in nodes.node_id], GC),
              df([meta_row("g" + n) for n in nodes.node_id], MC))
    r = mr.metrics(nodes, edges, COMPS, *single)
    assert r["singleton_share"] == 1.0 and r["singleton_node_share"] == 1.0
    assert r["long_edge_share"] == 1.0 and r["size_balance"] == pytest.approx(1.0)
    assert r["comp_redraw"] == 1.0                            # every box is pure: drawn inside one compartment
    assert r["top_boxes"] == 10


def test_localisation_levels_boundary_merge_and_singleton_exemption():
    nodes, edges = toy()
    part = {"A_e": "Le", "T1": "Lcb", "A_c": "Lc", "R1": "Lc", "B_c": "Lc",
            "R2": "Lc", "C_c": "Lc", "R3": "Lc", "D_c": "Lc", "X_c": "Lc"}
    g = df([[n, p, "", "", "S", ""] for n, p in sorted(part.items())], GC)
    m = df([meta_row("Le", label="extracellular", basis="localisation", source="compartment", xref="e", comp="e"),
            meta_row("Lc", label="cytosol", basis="localisation", source="compartment", xref="c", comp="c"),
            meta_row("Lcb", label="cytosol", basis="localisation", source="compartment", xref="c", comp="c@boundary")],
           MC)
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["n_groups"] == 3
    assert r["top_boxes"] == 2                                # c and c@boundary are one drawn compartment
    assert r["top_all_localisation"] == 1.0
    assert r["singleton_share"] == 0.0                        # Le and Lcb: the only box of their place
    assert r["singleton_node_share"] == pytest.approx(0.2)    # reported all the same
    assert r["named_share"] == 1.0 and r["backed_share"] == 1.0
    assert r["sibling_dup_share"] == 0.0                      # "cytosol" twice is the same compartment
    assert r["distinct_share"] == 1.0
    assert r["comp_redraw"] == 1.0
    assert r["depth_mean"] == 0.0                             # localisation levels are not counted
    assert r["long_edge_share"] == pytest.approx(0.2)         # T1-A_e and T1-A_c cross compartment borders


def test_localisation_below_function_is_drawn_once_per_parent():
    nodes, edges = toy()
    part = {"A_e": ("P", "Pe"), "T1": ("P", "Pc"), "A_c": ("P", "Pc"), "R1": ("P", "Pc"), "B_c": ("P", "Pc"),
            "R2": ("Q", "Qc"), "C_c": ("Q", "Qc"), "R3": ("Q", "Qc"), "D_c": ("Q", "Qc"), "X_c": ("Q", "Qc")}
    g = df([[n, a, b, "", "S", ""] for n, (a, b) in sorted(part.items())], GC)
    loc = lambda gid, par, c: meta_row(gid, 2, par, label=c, basis="localisation", source="compartment", xref=c, comp=c)
    m = df([meta_row("P", label="uptake"), meta_row("Q", label="conversion"),
            loc("Pe", "P", "e"), loc("Pc", "P", "c"), loc("Qc", "Q", "c")], MC)
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["loc_box_repeat"] == pytest.approx(1.5)          # e under P; c under P and under Q
    assert r["comp_redraw"] == pytest.approx(1.5)             # same answer from the node placement alone
    assert r["depth_mean"] == 1.0 and r["depth_levels_mean"] == 2.0 and r["depth_max"] == 2
    assert r["fan_ok_share"] == 1.0                           # Q's only child is a compartment: allowed
    assert r["only_child_share"] == pytest.approx(0.5)
    assert r["global_dup_share"] == pytest.approx(2 / 5)      # "c" under P and under Q: not siblings
    assert r["sibling_dup_share"] == 0.0
    assert r["singleton_share"] == pytest.approx(0.0)         # Pe is the only box of e
    # each c box holds about half of the c nodes of the MAP but all c nodes of its PARENT
    assert r["name_low_recall_share"] == 0.0


def test_unplaced_and_unknown_nodes_are_reported_not_fatal():
    nodes, edges = toy()
    g, m = toy_groups()
    g = pd.concat([g[g.node_id != "D_c"], df([["ghost", "P", "", "", "S", ""]], GC)], ignore_index=True)
    r = mr.metrics(nodes, edges, COMPS, g, m)
    assert r["n_nodes"] == 9 and r["placed_share"] == pytest.approx(0.9)
    empty = mr.metrics(nodes, edges, COMPS, g.iloc[0:0], m.iloc[0:0])
    assert empty["n_nodes"] == 0 and all(math.isnan(empty[spec[0]]) for spec in mr.RULES.values())


# ------------------------------------------------------------------------------------------ names
def named_model():
    """8 reactions: four EC 1.1.1.x, two EC 2.7.1.x, two without EC; 24 nodes in all."""
    rows = [(f"R{i}", "reaction", "c", "inside", ec, "SBO:0000000>SBO:0000176") for i, ec in
            enumerate(["1.1.1.1", "1.1.1.2", "1.1.1.3", "1.1.1.4", "2.7.1.1", "2.7.1.2", "", ""])]
    rows += [(f"S{i}", "species", "c", "inside") for i in range(16)]
    nodes = nodes_df(rows)
    edges = df([(f"R{i}", f"S{2 * i}", "reactant") for i in range(8)]
               + [(f"R{i}", f"S{2 * i + 1}", "product") for i in range(8)], ["reaction_id", "species_id", "role"])
    return nodes, edges


def named_groups(meta_rows, assign):
    part = {}
    for i in range(8):
        for n in (f"R{i}", f"S{2 * i}", f"S{2 * i + 1}"):
            part[n] = assign[i]
    return df([[n, p, "", "", "S", ""] for n, p in sorted(part.items())], GC), df(meta_rows, MC)


ONE_COMP = df([["c", "cytosol", ""]], ["compartment_id", "name", "parent_id"])


def test_names_duplicates_recall_and_backing():
    nodes, edges = named_model()
    # U holds 1 of the 4 EC-1 reactions, V the other 3; both are called by an EC-1 class.
    g, m = named_groups([meta_row("U", label="EC 1", basis="reaction_type", source="ec", xref="ec-code:1"),
                         meta_row("V", label="EC 1", basis="reaction_type", source="ec", xref="ec-code:1.1"),
                         meta_row("W", label="EC 2.7", basis="reaction_type", source="ec", xref="ec-code:2.7"),
                         meta_row("Z")],
                        ["U", "V", "V", "V", "W", "W", "Z", "Z"])
    r = mr.metrics(nodes, edges, ONE_COMP, g, m)
    assert r["top_groups"] == 4                               # 24 nodes: the rule applies
    assert r["named_share"] == pytest.approx(3 / 4)
    assert r["named_node_share"] == pytest.approx(18 / 24)
    assert r["sibling_dup_share"] == pytest.approx(2 / 3)     # "EC 1" twice
    assert r["name_low_recall_share"] == pytest.approx(1 / 3)  # U holds 1/4 of the EC-1 reactions
    assert r["name_recall_mean"] == pytest.approx((0.25 + 0.75 + 1.0) / 3)
    assert r["distinct_share"] == pytest.approx(1 / 3)        # only W is unique and informative
    assert r["backed_share"] == 1.0                           # each class is carried by a member
    assert math.isnan(r["comp_redraw"])                       # one compartment: rule not applicable
    # same partition, unique names, but one name is backed by nothing in the group
    g, m = named_groups([meta_row("U", label="alcohol dehydrogenase", source="member", xref="R0"),
                         meta_row("V", label="oxidoreductases", basis="reaction_type", source="ec", xref="ec-code:1.1"),
                         meta_row("W", label="glycolysis", basis="function", source="go", xref="GO:0006096"),
                         meta_row("Z", label="process", basis="reaction_type", source="sbo", xref="SBO:0000176")],
                        ["U", "V", "V", "V", "W", "W", "Z", "Z"])
    r = mr.metrics(nodes, edges, ONE_COMP, g, m)
    assert r["named_share"] == 1.0 and r["sibling_dup_share"] == 0.0
    # backed: V (EC carried) and Z (SBO carried). U is a member name, W's GO id is on no node.
    assert r["backed_share"] == pytest.approx(2 / 4)
    assert r["member_named_share"] == pytest.approx(1 / 4)
    # informative: U (member, not checked), V (3/4). W recall 0, Z holds 2 of 8 SBO:0000176 reactions.
    assert r["name_low_recall_share"] == pytest.approx(2 / 4)
    assert r["distinct_share"] == pytest.approx(2 / 4)


def test_name_equal_to_parent_and_fan_out():
    nodes, edges = named_model()
    part, meta = {}, [meta_row("A", label="metabolism"), meta_row("B", label="other")]
    for i in range(8):
        l1 = "A" if i < 6 else "B"
        l2 = f"A{i}" if i < 6 else "B0"
        for n in (f"R{i}", f"S{2 * i}", f"S{2 * i + 1}"):
            part[n] = (l1, l2)
    meta += [meta_row(f"A{i}", 2, "A", label=("metabolism" if i == 0 else f"step {i}")) for i in range(6)]
    meta += [meta_row("B0", 2, "B", label="rest")]
    g = df([[n, a, b, "", "S", ""] for n, (a, b) in sorted(part.items())], GC)
    r = mr.metrics(nodes, edges, ONE_COMP, g, df(meta, MC))
    assert r["n_groups"] == 9 and r["n_groups_L2"] == 7
    assert r["same_as_ancestor_share"] == pytest.approx(1 / 9)   # A0 repeats its parent's name
    assert r["distinct_share"] == pytest.approx(8 / 9)
    assert r["fan_ok_share"] == pytest.approx(0.5)            # A has 6 children (ok), B has 1 (not ok)
    assert r["only_child_share"] == pytest.approx(0.5)
    assert r["fan_max"] == 6
    assert r["depth_mean"] == 2.0 and r["depth_needed"] == 1 and r["depth_excess"] == 1.0
    assert r["oversize_node_share"] == 0.0 and r["box_size_max"] == 6


def test_fan_out_upper_limit_and_box_size():
    n = 13 * 4 + 60                                            # 13 small L2 boxes and one box of 60
    rows = [(f"N{i:03d}", "species", "c", "inside") for i in range(n)]
    nodes = nodes_df(rows)
    edges = df([], ["reaction_id", "species_id", "role"])
    g, meta = [], [meta_row("A"), meta_row("BIG")]
    for i in range(52):
        g.append([f"N{i:03d}", "A", f"A{i // 4:02d}", "", "S", ""])
    for i in range(52, n):
        g.append([f"N{i:03d}", "BIG", "", "", "S", ""])
    meta += [meta_row(f"A{j:02d}", 2, "A") for j in range(13)]
    r = mr.metrics(nodes, edges, ONE_COMP, df(g, GC), df(meta, MC))
    assert r["fan_ok_share"] == 0.0 and r["fan_max"] == 13     # 13 children > 12
    assert r["oversize_node_share"] == pytest.approx(60 / 112)
    assert r["box_over100_node_share"] == 0.0
    assert r["box_size_max"] == 60 and r["box_size_median"] == 4
    assert math.isnan(r["long_edge_share"])                    # no edges
    # effective number of boxes: 13 boxes of 4 and one of 60
    p = [4 / 112] * 13 + [60 / 112]
    assert r["size_balance"] == pytest.approx(math.exp(-sum(x * math.log(x) for x in p)) / 14)


def span_model(box):
    """Two boxes a, b of `box` reactions + `box` species each. Every reaction has one species in its own
    box; reactions 0..5 of box a also use species 0..5 of box b. 2*box + 6 edges."""
    rows, e = [], []
    for side in "ab":
        for i in range(box):
            rows += [(f"R{side}{i:02d}", "reaction", "c", "inside"), (f"S{side}{i:02d}", "species", "c", "inside")]
            e.append((f"R{side}{i:02d}", f"S{side}{i:02d}", "reactant"))
    e += [(f"Ra{i:02d}", f"Sb{i:02d}", "product") for i in range(6)]
    return nodes_df(rows), df(e, ["reaction_id", "species_id", "role"])


def test_long_edge_needs_a_common_box_of_readable_span():
    for box, expect in ((30, 6 / 66), (25, 0.0)):             # parent of 120 nodes vs parent of 100 nodes
        nodes, edges = span_model(box)
        g = df([[n, "A", "A" + n[1], "", "S", ""] for n in sorted(nodes.node_id)], GC)
        m = df([meta_row("A"), meta_row("Aa", 2, "A"), meta_row("Ab", 2, "A")], MC)
        r = mr.metrics(nodes, edges, ONE_COMP, g, m)
        assert r["long_edge_share"] == pytest.approx(expect)  # a/b edges: common box A only
        assert r["l1_cross_share"] == 0.0
        assert r["cut_edge_share"] == pytest.approx(6 / (2 * box + 6))
        assert r["clone_demand"] == pytest.approx((6 if expect else 0) / (2 * box))
        flat = df([[n, "A" + n[1], "", "", "S", ""] for n in sorted(nodes.node_id)], GC)
        r = mr.metrics(nodes, edges, ONE_COMP, flat, df([meta_row("Aa"), meta_row("Ab")], MC))
        assert r["long_edge_share"] == pytest.approx(6 / (2 * box + 6))   # no common box at all
        one = df([[n, "A", "", "", "S", ""] for n in sorted(nodes.node_id)], GC)
        r = mr.metrics(nodes, edges, ONE_COMP, one, df([meta_row("A")], MC))
        assert r["long_edge_share"] == (1.0 if 4 * box > mr.LONG_SPAN else 0.0)


def test_depth_needed_and_scores():
    assert [mr.depth_needed(x) for x in (10, 175, 176, 1225, 1226)] == [1, 1, 2, 2, 3]
    s = mr._score_one
    assert s(0.05, "low", 0.10, 0.30, 0.70) == 1.0
    assert s(0.30, "low", 0.10, 0.30, 0.70) == pytest.approx(0.5)
    assert s(0.50, "low", 0.10, 0.30, 0.70) == pytest.approx(0.25)
    assert s(0.90, "low", 0.10, 0.30, 0.70) == 0.0
    assert s(0.75, "high", 0.90, 0.60, 0.20) == pytest.approx(0.75)
    rng = ((3, 12), (2, 20), (1, 40))
    assert [s(k, "range", *rng) for k in (1, 2, 3, 12, 20, 40)] == pytest.approx([0, 0.5, 1, 1, 0.5, 0])
    assert math.isnan(s(float("nan"), "low", 0, 1, 2))
    v = mr.verdicts({"top_groups": 1.0, "oversize_node_share": 0.2, "singleton_share": 0.0})
    assert v["RD01_TOP"] == "fail" and v["RD02_BOXSIZE"] == "warn" and v["RD03_SINGLETON"] == "pass"
    assert v["RD04_NAMED"] == "na"
    # a top level made only of compartments has no lower limit on their number
    assert mr.verdicts({"top_groups": 1.0, "top_all_localisation": 1.0})["RD01_TOP"] == "pass"
    assert mr.verdicts({"top_groups": 30.0, "top_all_localisation": 1.0})["RD01_TOP"] == "fail"


def test_reference_is_optional_and_measures_scatter():
    nodes, edges = toy()
    g, m = toy_groups()
    assert math.isnan(mr.metrics(nodes, edges, COMPS, g, m)["ref_scatter"])
    ref = df([("A_e", "up"), ("A_c", "up"), ("B_c", "up"),                 # all in P
              ("C_c", "down"), ("D_c", "down"), ("X_c", "down"),           # all in Q
              ("A_c", "split"), ("B_c", "split"), ("C_c", "split"), ("D_c", "split")], ["node_id", "label"])
    r = mr.metrics(nodes, edges, COMPS, g, m, reference=ref)
    assert r["ref_labelled_share"] == pytest.approx(0.6)
    # up -> 1 group, down -> 1 group, split -> 2 groups evenly; weights 3, 3, 4
    assert r["ref_scatter"] == pytest.approx((3 * 1 + 3 * 1 + 4 * 2) / 10)
    # P holds up x3 + split x2, Q holds down x3 + split x2: the same effective number in both
    eff = math.exp(-(0.6 * math.log(0.6) + 0.4 * math.log(0.4)))
    assert r["ref_mix"] == pytest.approx(eff)


# ------------------------------------------------------------------------------------------ panel
MODELS = ["brenda_merged", "e_coli_core", "R-HSA-77289", "skm_pss", "skm_ja", "skm_sa"]


def _panel():
    here = pathlib.Path(__file__).resolve().parent
    cands = [os.environ.get("PANEL_FEATURES", "")]
    for base in (pathlib.Path.cwd(), here, here.parent, here.parent.parent):
        cands += [base / "panel" / "features", base / "bundle" / "panel" / "features"]
    for c in cands:
        if c and (pathlib.Path(c) / MODELS[1] / "nodes.tsv").exists():
            return pathlib.Path(c)
    return None


PANEL = _panel()
needs_panel = pytest.mark.skipif(PANEL is None, reason="panel/features not found (set PANEL_FEATURES)")


def _frames(rows, meta):
    return (pd.DataFrame(rows, columns=GC).astype(str).sort_values("node_id").reset_index(drop=True),
            pd.DataFrame(meta, columns=MC).astype(str).drop_duplicates("group_id").sort_values("group_id")
            .reset_index(drop=True))


def all_in_one(nodes):
    return _frames([[n, "ALL", "", "", "S", ""] for n in nodes.node_id], [meta_row("ALL")])


def singletons(nodes):
    return _frames([[n, "g_" + n, "", "", "S", ""] for n in nodes.node_id],
                   [meta_row("g_" + n) for n in nodes.node_id])


def per_compartment(nodes, comps):
    name = dict(zip(comps.compartment_id, comps.name))
    return _frames([[n, "loc_" + c, "", "", "S", ""] for n, c in zip(nodes.node_id, nodes.loc_space)],
                   [meta_row("loc_" + c, label=name.get(c, c), basis="localisation", source="compartment",
                             xref=c, comp=c) for c in sorted(set(nodes.loc_space))])


def random_realistic(nodes, seed=0, mean_size=25):
    """About N/25 groups with log-normal sizes (a few large, many small), nodes assigned at random."""
    rng = np.random.default_rng(seed)
    k = max(2, round(len(nodes) / mean_size))
    w = rng.lognormal(0.0, 1.0, k)
    lab = rng.choice(k, size=len(nodes), p=w / w.sum())
    return _frames([[n, f"r{x:03d}", "", "", "S", ""] for n, x in zip(nodes.node_id, lab)],
                   [meta_row(f"r{x:03d}", basis="topology") for x in sorted(set(lab))])


def louvain_named(nodes, edges, seed=1, split=50):
    """A plausible grouping: Louvain communities of the species-reaction graph, communities above 50
    nodes split once more, every group named after its highest-degree member."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(sorted(nodes.node_id))
    G.add_edges_from(sorted(zip(edges.reaction_id, edges.species_id)))
    deg, nm = dict(G.degree()), dict(zip(nodes.node_id, nodes.name))

    def comm(sub):
        cs = nx.community.louvain_communities(G.subgraph(sub), seed=seed)
        return [sorted(c) for c in sorted(cs, key=lambda c: sorted(c)[0])]

    rows, meta = {}, []

    def name(gid, level, parent, members, used):
        cand = sorted(members, key=lambda v: (-deg[v], v))
        hub = next((v for v in cand if nm.get(v, v) not in used), cand[0])
        used.add(nm.get(hub, hub))
        meta.append(meta_row(gid, level, parent, label=nm.get(hub, hub) or hub, basis="topology",
                             source="member", xref=hub))

    used = set()
    for i, c in enumerate(comm(sorted(G))):
        gid = f"t{i:02d}"
        name(gid, 1, "", c, used)
        if len(c) > split:
            used2 = set()
            for j, cc in enumerate(comm(c)):
                name(f"{gid}_{j:02d}", 2, gid, cc, used2)
                for v in cc:
                    rows[v] = [v, gid, f"{gid}_{j:02d}", "", "S", ""]
        else:
            for v in c:
                rows[v] = [v, gid, "", "", "S", ""]
    return _frames([rows[v] for v in sorted(rows)], meta)


@pytest.fixture(scope="module", params=MODELS)
def panel_model(request):
    d = PANEL / request.param
    rd = lambda p: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    nodes, edges, comps = rd(d / "nodes.tsv"), rd(d / "edges.tsv"), rd(d / "compartments.tsv")
    G = {"one": all_in_one(nodes), "single": singletons(nodes), "comp": per_compartment(nodes, comps),
         "random": random_realistic(nodes), "louvain": louvain_named(nodes, edges)}
    M = {k: mr.metrics(nodes, edges, comps, g, m) for k, (g, m) in G.items()}
    return request.param, nodes, edges, comps, M


@needs_panel
def test_panel_all_in_one(panel_model):
    _name, nodes, edges, _comps, M = panel_model
    m, v = M["one"], mr.verdicts(M["one"])
    assert m["n_nodes"] == len(nodes) and m["n_groups"] == 1 and m["top_groups"] == 1
    assert m["l1_cross_share"] == 0.0 and m["cut_edge_share"] == 0.0     # no edge leaves the box ...
    assert m["long_edge_share"] == 1.0          # ... yet every edge is long: the only common box has > 100 nodes
    assert m["oversize_node_share"] == 1.0 and m["largest_box_share"] == 1.0
    assert m["singleton_share"] == 0.0 and m["named_share"] == 0.0 and m["cloned_species_share"] == 0.0
    for rid in ("RD03_SINGLETON", "RD06_CLONE", "RD12_DEPTH"):           # rewarded here ...
        assert v[rid] == "pass"
    assert v["RD10_COMPDRAW"] in ("pass", "na")
    for rid in ("RD01_TOP", "RD02_BOXSIZE", "RD04_NAMED", "RD05_LONGEDGE"):   # ... and stopped here
        assert v[rid] == "fail"


@needs_panel
def test_panel_singletons(panel_model):
    _name, nodes, edges, _comps, M = panel_model
    m, v = M["single"], mr.verdicts(M["single"])
    assert m["n_groups"] == len(nodes) and m["top_groups"] == len(nodes)
    assert m["singleton_share"] == 1.0 and m["long_edge_share"] == 1.0
    assert m["oversize_node_share"] == 0.0 and m["size_balance"] == pytest.approx(1.0)
    assert v["RD02_BOXSIZE"] == "pass" and v["RD11_BALANCE"] == "pass"   # rewarded here ...
    assert v["RD01_TOP"] == "fail" and v["RD03_SINGLETON"] == "fail"     # ... and stopped here
    assert v["RD05_LONGEDGE"] == "fail"
    species = int((nodes.kind == "species").sum())
    pairs = len(set(zip(edges.reaction_id, edges.species_id)))
    assert m["clone_demand"] == pytest.approx(pairs / species)           # one clone per distinct edge


@needs_panel
def test_panel_per_compartment(panel_model):
    _name, nodes, _edges, comps, M = panel_model
    m = M["comp"]
    n_spaces = nodes.loc_space.nunique()
    assert m["n_groups"] == n_spaces and m["top_boxes"] == n_spaces
    assert m["named_share"] == 1.0 and m["backed_share"] == 1.0 and m["distinct_share"] == 1.0
    assert m["depth_mean"] == 0.0 and m["top_all_localisation"] == 1.0
    if n_spaces > 1:
        assert m["comp_redraw"] == 1.0 and m["mixed_box_node_share"] == 0.0
    else:
        assert math.isnan(m["comp_redraw"])
    biggest = nodes.loc_space.value_counts().iloc[0]
    assert m["box_size_max"] == biggest
    assert m["oversize_node_share"] == pytest.approx(
        nodes.loc_space.value_counts().loc[lambda s: s > mr.BOX_MAX].sum() / len(nodes))


@needs_panel
def test_panel_random_is_poor(panel_model):
    _name, _nodes, _edges, _comps, M = panel_model
    m = M["random"]
    assert m["long_edge_share"] > 0.6 and m["modularity_L1"] < 0.1
    assert m["named_share"] == 0.0
    assert mr.verdicts(m)["RD05_LONGEDGE"] in ("warn", "fail")


@needs_panel
def test_panel_no_degenerate_grouping_ranks_first(panel_model):
    name, _nodes, _edges, _comps, M = panel_model
    score = {k: mr.readability_score(m) for k, m in M.items()}
    assert all(0.0 <= s <= 1.0 for s in score.values())
    best = max(score, key=score.get)
    assert best not in ("one", "single"), (name, score)
    assert score["one"] < 0.3 and score["single"] < 0.3, (name, score)
    for k in ("comp", "random", "louvain"):
        if k == "random" and name == "skm_pss":
            continue                        # 73 random top-level groups are as unreadable as the degenerate cases
        assert score[k] > max(score["one"], score["single"]), (name, k, score)
    assert score["louvain"] > score["random"], (name, score)


@needs_panel
def test_panel_cloning_everything_does_not_pay(panel_model):
    """Cloning every species into the home box of every reaction it is long-connected to removes the
    long edges a clone can remove; RD06 takes the gain back."""
    name, nodes, edges, comps, M = panel_model
    g, meta = louvain_named(nodes, edges)
    home = {n: (b or a) for n, a, b in zip(g.node_id, g.L1, g.L2)}
    path = {n: {a, b} - {""} for n, a, b in zip(g.node_id, g.L1, g.L2)}
    want = {}
    for r, s in zip(edges.reaction_id, edges.species_id):
        if home[r] not in path[s]:
            want.setdefault(s, set()).add(home[r])
    g2 = g.copy()
    g2["clone_in"] = [("|".join(sorted(want[n])) if n in want else "") for n in g2.node_id]
    m1, m2 = M["louvain"], mr.metrics(nodes, edges, comps, g2, meta)
    assert m2["clone_useless_share"] == 0.0
    assert m2["long_edge_share_raw"] == pytest.approx(m1["long_edge_share"])   # the grouping is the same
    assert m2["long_edge_share"] <= 0.02 < m1["long_edge_share"]
    assert m2["cut_edge_share"] == 0.0
    assert m2["clone_load"] >= m1["clone_demand"] > 0
    s1, s2 = mr.rule_scores(m1), mr.rule_scores(m2)
    gain = mr.RULES["RD05_LONGEDGE"][5] * (s2["RD05_LONGEDGE"] - s1["RD05_LONGEDGE"])
    loss = mr.RULES["RD06_CLONE"][5] * (s1["RD06_CLONE"] - s2["RD06_CLONE"])
    assert gain <= loss + 1e-9, (name, gain, loss)
    assert mr.readability_score(m2) <= mr.readability_score(m1) + 1e-9
