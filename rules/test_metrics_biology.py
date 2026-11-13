"""Run:
    pytest -q test_metrics_biology.py
    GROUPING_PANEL=path/to/panel/features pytest -q test_metrics_biology.py
    PowerShell:  $env:GROUPING_PANEL = "path\\to\\panel\\features"; pytest -q test_metrics_biology.py

Tests of metrics_biology.py.

Part 1 uses small hand-made models whose correct metric values are known by construction.
Part 2 builds degenerate and control groupings on all six panel models (one group, every node its own group,
every node its own compartment box, a random grouping with realistic sizes, one group per compartment, and
three non-trivial controls) and checks that no degenerate grouping passes the biology rules or comes out on
top. The panel is looked for in $GROUPING_PANEL, then in panel/features and bundle/panel/features next to this
file, in the working directory and up to three parents; part 2 is skipped with a message if it is not found.
networkx is needed for part 2 (Louvain communities with a fixed seed, for the controls).
"""
import math
import os
import pathlib
import random
from collections import Counter, defaultdict

import pandas as pd
import pytest

import metrics_biology as mb

NODE_COLS = ["node_id", "kind", "name", "compartment", "compartments", "modifier_compartments", "loc_space", "loc_kind",
             "loc_outer", "sbo", "sbo_name", "sbo_path", "sbo_path_names", "ec", "go", "xrefs", "notes", "reversible",
             "n_reactants", "n_products", "n_modifiers", "degree", "boundary"]
COMP_COLS = ["compartment_id", "name", "canonical", "go_id", "kind", "parent_id", "encloses_id", "virtual", "depth",
             "source", "n_species"]
GCOLS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
MCOLS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]


# ------------------------------------------------------------------------------------- table builders
def nodes_df(rows):
    return pd.DataFrame([dict(dict.fromkeys(NODE_COLS, ""), **r) for r in rows], columns=NODE_COLS).astype(str)


def edges_df(pairs):
    return pd.DataFrame([dict(reaction_id=r, species_id=s, role="reactant", stoichiometry="1", sbo="", sbo_name="")
                         for r, s in pairs],
                        columns=["reaction_id", "species_id", "role", "stoichiometry", "sbo", "sbo_name"]).astype(str)


def comps_df(rows):
    return pd.DataFrame([dict(dict.fromkeys(COMP_COLS, ""), kind="space", **r) for r in rows], columns=COMP_COLS).astype(str)


def frames(rows, meta):
    """rows: node id -> (L1, L2, L3, rule[, clone_in]); meta: group id -> dict of meta columns (level 1 by default)."""
    g = pd.DataFrame([dict(node_id=n, L1=v[0], L2=v[1], L3=v[2], rule=v[3], clone_in=v[4] if len(v) > 4 else "")
                      for n, v in sorted(rows.items())], columns=GCOLS).astype(str)
    used = {x for v in rows.values() for x in v[:3] if x}
    m = pd.DataFrame([dict(dict(group_id=k, level="1", parent="", label="", basis="mixed", name_source="none",
                                name_xref="", compartment_id=""), **v) for k, v in sorted(meta.items()) if k in used],
                     columns=MCOLS).astype(str)
    return g, m


def flat(assign, meta=None):
    """assign: node id -> group id, one level."""
    meta = dict(meta or {})
    for g in set(assign.values()):
        meta.setdefault(g, {})
    return frames({n: (g, "", "", "R1") for n, g in assign.items()}, meta)


ONE_COMP = comps_df([dict(compartment_id="c", name="cell", source="single_compartment")])


def rx(i, sbo_path="", ec="", go="", space="c", kind="inside", outer="", comps=None, xrefs="", notes=""):
    return dict(node_id=i, kind="reaction", compartment=space, compartments=comps or space, loc_space=space,
                loc_kind=kind, loc_outer=outer, sbo=sbo_path.split(">")[-1], sbo_path=sbo_path, ec=ec, go=go,
                xrefs=xrefs, notes=notes, degree="1")


def sp(i, space="c", kind="inside", outer="", comp=None, xrefs="", notes=""):
    return dict(node_id=i, kind="species", compartment=comp or space, compartments=comp or space, loc_space=space,
                loc_kind=kind, loc_outer=outer, xrefs=xrefs, notes=notes, degree="1")


def run(nodes, groups, meta, edges=None, comps=None, **kw):
    return mb.metrics(nodes, edges if edges is not None else edges_df([]), comps if comps is not None else ONE_COMP,
                      groups, meta, **kw)


# =========================================================================== part 1: values known by construction
# ---- BIO-2 homogeneity
FOUR = nodes_df([rx("r1", "X>A"), rx("r2", "X>A"), rx("r3", "X>B"), rx("r4", "X>B")])


def test_homogeneity_known_value():
    # 4 reactions, terms A, A, B, B. IC(A) = IC(B) = ln(4/2)/ln(4) = 0.5, IC(X) = 0.
    # same-type groups: S_in = 0.5; chance: S_all = (0.5 + 0 + 0)/3 = 1/6; hom = (1/2 - 1/6)/(1 - 1/6) = 0.4
    m = run(FOUR, *flat(dict(r1="g1", r2="g1", r3="g2", r4="g2")))
    assert m["hom"] == pytest.approx(0.4)
    assert m["hom_sbo"] == pytest.approx(0.4)
    assert math.isnan(m["hom_go"]) and math.isnan(m["hom_ec"])


def test_homogeneity_mixed_groups_are_below_chance():
    # every group pairs an A with a B: S_in = 0 -> (0 - 1/6)/(5/6) = -0.2
    m = run(FOUR, *flat(dict(r1="g1", r3="g1", r2="g2", r4="g2")))
    assert m["hom"] == pytest.approx(-0.2)


def test_homogeneity_degenerate_groupings_are_exactly_chance():
    assert run(FOUR, *flat(dict.fromkeys(["r1", "r2", "r3", "r4"], "g")))["hom"] == pytest.approx(0.0)
    assert run(FOUR, *flat({n: "g" + n for n in ["r1", "r2", "r3", "r4"]}))["hom"] == pytest.approx(0.0)


def test_homogeneity_ec_hierarchy_and_unannotated_nodes():
    # 4 reactions with EC + 4 species without annotation: eligible nodes N = 8.
    # shared class ec:1 (2 of 8 nodes): IC = ln(8/2)/ln(8) = 2/3. same-class groups: S_in = 2/3,
    # S_all = (2/3)/3 = 2/9 -> hom = (2/3 - 2/9)/(1 - 2/9) = 4/7. Species do not enter (no term).
    nodes = nodes_df([rx("r1", ec="1.1.1.1"), rx("r2", ec="1.1.1.2"), rx("r3", ec="2.7.1.1"), rx("r4", ec="2.7.1.2"),
                      sp("s1"), sp("s2"), sp("s3"), sp("s4")])
    g = flat(dict(r1="g1", r2="g1", s1="g1", s2="g1", r3="g2", r4="g2", s3="g2", s4="g2"))
    m = run(nodes, *g)
    assert m["hom"] == pytest.approx(4 / 7)
    assert m["hom_ec"] == pytest.approx(4 / 7)
    assert m["annotated_share"] == pytest.approx(0.5)


def test_homogeneity_uninformative_term_gives_nan_or_zero():
    # every reaction carries the same term: IC = 0 everywhere, nothing to share -> 0 (not 1)
    nodes = nodes_df([rx("r1", "X>A"), rx("r2", "X>A"), rx("r3", "X>A")])
    assert run(nodes, *flat(dict(r1="g1", r2="g1", r3="g2")))["hom"] == pytest.approx(0.0)


def test_homogeneity_uses_leaf_level():
    g, m = frames(dict(r1=("top", "a", "", "R1"), r2=("top", "a", "", "R1"), r3=("top", "b", "", "R1"),
                       r4=("top", "b", "", "R1")),
                  dict(top={}, a=dict(level="2", parent="top"), b=dict(level="2", parent="top")))
    out = run(FOUR, g, m)
    assert out["hom"] == pytest.approx(0.4) and out["hom_L1"] == pytest.approx(0.0)


# ---- BIO-1 agreement with a reference
ABCD = nodes_df([sp("a"), sp("b"), sp("c"), sp("d")])


def ref(pairs):
    return pd.DataFrame(pairs, columns=["node_id", "label"]).astype(str)


def test_reference_perfect_and_baselines():
    r = ref([("a", "X"), ("b", "X"), ("c", "Y"), ("d", "Y")])
    m = run(ABCD, *flat(dict(a="g1", b="g1", c="g2", d="g2")), reference=r)
    assert m["ref_f1"] == pytest.approx(1.0)
    # one group: P = 1/2, R = 1 -> 2/3; singletons: P = 1, R = 1/2 -> 2/3; one compartment = one group
    assert m["ref_f1_one_group"] == pytest.approx(2 / 3)
    assert m["ref_f1_singletons"] == pytest.approx(2 / 3)
    assert m["ref_baseline_f1"] == pytest.approx(2 / 3)
    assert m["ref_skill"] == pytest.approx(1.0)
    assert m["ref_labelled_share"] == pytest.approx(1.0)


def test_reference_wrong_grouping_has_negative_skill():
    r = ref([("a", "X"), ("b", "X"), ("c", "Y"), ("d", "Y")])
    m = run(ABCD, *flat(dict(a="g1", c="g1", b="g2", d="g2")), reference=r)
    assert m["ref_precision"] == pytest.approx(0.5) and m["ref_recall"] == pytest.approx(0.5)
    assert m["ref_skill"] == pytest.approx((0.5 - 2 / 3) / (1 - 2 / 3))


def test_reference_degenerate_groupings_have_no_skill():
    r = ref([("a", "X"), ("b", "X"), ("c", "Y"), ("d", "Y")])
    assert run(ABCD, *flat(dict.fromkeys("abcd", "g")), reference=r)["ref_skill"] == pytest.approx(0.0)
    assert run(ABCD, *flat({n: "g" + n for n in "abcd"}), reference=r)["ref_skill"] == pytest.approx(0.0)


def test_reference_multilabel_needs_a_clone_for_full_recall():
    # b belongs to X and Y. Home with a, a copy with c: every shared label is matched by a shared group.
    nodes = nodes_df([sp("a"), sp("b"), sp("c"), sp("u")])
    r = ref([("a", "X"), ("b", "X"), ("b", "Y"), ("c", "Y")])            # u is unlabelled
    g, m = frames(dict(a=("g1", "", "", "R1"), b=("g1", "", "", "R1", "g2"), c=("g2", "", "", "R1"),
                       u=("g2", "", "", "R1")), dict(g1={}, g2={}))
    out = run(nodes, g, m, reference=r)
    assert out["ref_precision"] == pytest.approx(1.0) and out["ref_recall"] == pytest.approx(1.0)
    assert out["ref_labelled_share"] == pytest.approx(0.75)
    # without the copy: recall a = 1, b = (1 + 1/2 + 0)/3 = 1/2, c = (1 + 0)/2 = 1/2 -> 2/3; precision 1; F = 0.8
    g2, m2 = flat(dict(a="g1", b="g1", c="g2", u="g2"))
    out = run(nodes, g2, m2, reference=r)
    assert out["ref_precision"] == pytest.approx(1.0)
    assert out["ref_recall"] == pytest.approx(2 / 3)
    assert out["ref_f1"] == pytest.approx(0.8)


def test_reference_takes_the_best_level():
    r = ref([("a", "X"), ("b", "X"), ("c", "Y"), ("d", "Y")])
    g, m = frames(dict(a=("top", "g1", "", "R1"), b=("top", "g1", "", "R1"), c=("top", "g2", "", "R1"),
                       d=("top", "g2", "", "R1")),
                  dict(top={}, g1=dict(level="2", parent="top"), g2=dict(level="2", parent="top")))
    out = run(ABCD, g, m, reference=r)
    assert out["ref_f1_L1"] == pytest.approx(2 / 3) and out["ref_f1_L2"] == pytest.approx(1.0)
    assert out["ref_best_level"] == 2.0 and out["ref_skill"] == pytest.approx(1.0)


def test_reference_absent_gives_nan():
    out = run(ABCD, *flat(dict.fromkeys("abcd", "g")))
    assert math.isnan(out["ref_skill"]) and math.isnan(out["ref_f1"])


# ---- LOC localisation: e (extracellular) around c (cytosol) around n (nucleus)
CELL = comps_df([dict(compartment_id="e", name="extracellular", go_id="GO:0005576"),
                 dict(compartment_id="c", name="cytosol", go_id="GO:0005829", parent_id="e"),
                 dict(compartment_id="n", name="nucleus", go_id="GO:0005634", parent_id="c")])
CELL_NODES = nodes_df([sp("se", "e"), sp("sc1", "c"), sp("sc2", "c"), sp("sn", "n"),
                       rx("rc", "X>A", space="c"), rx("rn", "X>A", space="n"),
                       rx("t_ec", "X>T", space="c", kind="boundary", outer="e", comps="c|e"),
                       rx("t_cn", "X>T", space="n", kind="boundary", outer="c", comps="c|n")])
CELL_EDGES = edges_df([("rc", "sc1"), ("rc", "sc2"), ("rn", "sn"), ("t_ec", "se"), ("t_ec", "sc1"),
                       ("t_cn", "sc2"), ("t_cn", "sn")])


def loc_meta(space, boundary=False, **kw):
    return dict(basis="localisation", label=space + (" boundary" if boundary else ""), name_source="compartment",
                name_xref=space, compartment_id=space + ("@boundary" if boundary else ""), **kw)


def per_compartment_small(boundary=True):
    assign, meta = {}, {}
    for r in CELL_NODES.itertuples(index=False):
        b = boundary and r.loc_kind == "boundary"
        g = f"loc_{r.loc_space}" + ("_b" if b else "")
        assign[r.node_id] = g
        meta[g] = loc_meta(r.loc_space, b)
    return flat(assign, meta)


def test_localisation_one_box_per_compartment_is_perfect():
    m = run(CELL_NODES, *per_compartment_small(), edges=CELL_EDGES, comps=CELL)
    for key in ("loc_expressed", "loc_boundary", "loc_nesting", "loc_repeat", "loc_span_ok", "loc_pure_leaf"):
        assert m[key] == pytest.approx(1.0), key
    assert all(v == "pass" for k, v in mb.judge(m).items() if k.startswith("LOC"))


def test_localisation_transport_mixed_into_the_interior():
    m = run(CELL_NODES, *per_compartment_small(boundary=False), edges=CELL_EDGES, comps=CELL)
    assert m["loc_boundary"] == pytest.approx(0.0) and m["loc_boundary_in_interior"] == pytest.approx(1.0)
    assert m["loc_expressed"] == pytest.approx(1.0)
    assert mb.judge(m)["LOC-3"] == "fail"


def test_localisation_not_expressed_in_one_group():
    m = run(CELL_NODES, *flat(dict.fromkeys(CELL_NODES.node_id, "g")), edges=CELL_EDGES, comps=CELL)
    assert m["loc_expressed"] == 0.0 and m["loc_nesting"] == 0.0 and m["loc_boundary"] == 0.0
    assert m["loc_span_ok"] == 0.0            # the group spans three compartments and is not split by place
    assert m["loc_span_connected"] == 1.0     # ... although its own transport reactions link them
    assert math.isnan(m["loc_repeat"])
    assert m["loc_pure_leaf"] == 0.0


def test_localisation_does_not_apply_to_one_compartment():
    m = run(FOUR, *flat(dict.fromkeys(["r1", "r2", "r3", "r4"], "g")))
    for key in ("loc_expressed", "loc_boundary", "loc_nesting", "loc_repeat", "loc_span_ok"):
        assert math.isnan(m[key]), key
    assert all(v == "na" for k, v in mb.judge(m).items() if k.startswith("LOC"))


def function_first():
    """F1 = {sc1, sc2, rc, t_cn, sn, rn} with boxes c, n and n@boundary; F2 = {se, t_ec} with boxes e and c@boundary."""
    rows = dict(sc1=("F1", "F1_c", "", "R1"), sc2=("F1", "F1_c", "", "R1"), rc=("F1", "F1_c", "", "R1"),
                sn=("F1", "F1_n", "", "R1"), rn=("F1", "F1_n", "", "R1"), t_cn=("F1", "F1_nb", "", "R1"),
                se=("F2", "F2_e", "", "R1"), t_ec=("F2", "F2_cb", "", "R1"))
    meta = dict(F1=dict(basis="function"), F2=dict(basis="function"),
                F1_c=loc_meta("c", level="2", parent="F1"), F1_n=loc_meta("n", level="2", parent="F1"),
                F1_nb=loc_meta("n", True, level="2", parent="F1"), F2_e=loc_meta("e", level="2", parent="F2"),
                F2_cb=loc_meta("c", True, level="2", parent="F2"))
    return frames(rows, meta)


def test_localisation_function_first_nesting_and_span():
    m = run(CELL_NODES, *function_first(), edges=CELL_EDGES, comps=CELL)
    assert m["loc_expressed"] == 1.0 and m["loc_boundary"] == 1.0
    # enclosing spaces to draw: c-nodes (sc1, sc2, rc, t_ec) need e: 4; n-nodes (sn, rn, t_cn) need c and e: 6.
    # drawn in the node's own scope: F1 has a c box but no e box -> sn, rn, t_cn get c (3), nobody in F1 gets e;
    # F2 has an e box -> t_ec gets e (1). 4 of 10.
    assert m["loc_nesting"] == pytest.approx(0.4)
    assert m["loc_repeat"] == pytest.approx(1.0)       # no compartment_id is used by two boxes
    assert m["loc_span_share"] == 1.0 and m["loc_span_connected"] == 1.0 and m["loc_span_ok"] == 1.0


def test_localisation_span_without_transport_is_not_accepted():
    # F = {sc1, sn} split into boxes c and n, but the transport t_cn that links c and n is in another group
    rows = {n: ("rest", "", "", "R1") for n in CELL_NODES.node_id}
    rows.update(sc1=("F", "F_c", "", "R1"), sn=("F", "F_n", "", "R1"))
    meta = dict(rest={}, F=dict(basis="function"), F_c=loc_meta("c", level="2", parent="F"),
                F_n=loc_meta("n", level="2", parent="F"))
    m = run(CELL_NODES, *frames(rows, meta), edges=CELL_EDGES, comps=CELL)
    # the 6 nodes of "rest" span e, c, n without boxes (bad); sc1 and sn sit in F, which lacks its transport (bad)
    assert m["loc_span_ok"] == pytest.approx(0.0)
    assert m["loc_span_connected"] == pytest.approx(6 / 8)


def test_localisation_repeated_compartment():
    # cytosol interior (sc1, sc2, rc) drawn as two boxes of 2 and 1 nodes; the other four boxes are whole.
    g, meta = per_compartment_small()
    g = g.copy()
    g.loc[g.node_id == "rc", "L1"] = "loc_c_2"
    meta = pd.concat([meta, pd.DataFrame([dict(group_id="loc_c_2", level="1", parent="", **loc_meta("c"))])])[MCOLS]
    meta = meta.sort_values("group_id").astype(str)
    m = run(CELL_NODES, g, meta, edges=CELL_EDGES, comps=CELL)
    # share of its compartment that a node's box holds: sc1, sc2 -> 2/3, rc -> 1/3, the other five -> 1
    together = (2 / 3 + 2 / 3 + 1 / 3 + 5) / 8
    assert m["loc_repeat"] == pytest.approx(1 / together)
    assert m["loc_boxes_max"] == 2.0
    # the nucleus nodes can no longer be nested in "the" cytosol box: 3 of their 6 enclosing spaces are ambiguous
    assert m["loc_nesting"] == pytest.approx(7 / 10)


def test_localisation_singletons_declared_as_compartments_fail_the_counter_rule():
    assign = {r.node_id: "g_" + r.node_id for r in CELL_NODES.itertuples(index=False)}
    meta = {"g_" + r.node_id: loc_meta(r.loc_space) for r in CELL_NODES.itertuples(index=False)}
    m = run(CELL_NODES, *flat(assign, meta), edges=CELL_EDGES, comps=CELL)
    assert m["loc_expressed"] == 1.0                       # LOC-1 alone is satisfied ...
    # ... LOC-2 is not: e holds 1 node (share 1), c holds 4 (share 1/4 each), n holds 3 (share 1/3 each)
    assert m["loc_repeat"] == pytest.approx(8 / (1 + 4 * 0.25 + 3 / 3))
    assert mb.judge(m)["LOC-2"] == "fail"
    assert math.isnan(m["name_lift"]) and m["name_grounded"] == 0.0   # names on single nodes do not count


# ---- FAIR names
def test_name_lift_true_informative_name():
    g, m = flat(dict(r1="g1", r2="g1", r3="g2", r4="g2"),
                dict(g1=dict(basis="reaction_type", label="A", name_source="sbo", name_xref="A"),
                     g2=dict(basis="reaction_type", label="B", name_source="sbo", name_xref="B")))
    out = run(FOUR, g, m)
    assert out["name_support"] == 1.0 and out["name_lift"] == 1.0 and out["name_traceable"] == 1.0


def test_name_lift_half_true_name():
    # "A" on a group holding one A and one B: support 1/2, prevalence of A = 1/2 -> lift 0
    g, m = flat(dict(r1="g1", r3="g1", r2="g2", r4="g2"),
                dict(g1=dict(label="A", name_source="sbo", name_xref="A")))
    out = run(FOUR, g, m)
    assert out["name_support"] == pytest.approx(0.5) and out["name_lift"] == pytest.approx(0.0)


def test_name_that_holds_for_every_node_is_vacuous():
    g, m = flat(dict.fromkeys(["r1", "r2", "r3", "r4"], "g"), dict(g=dict(label="X", name_source="sbo", name_xref="X")))
    out = run(FOUR, g, m)
    assert math.isnan(out["name_lift"]) and out["name_vacuous_share"] == 1.0
    assert out["name_coverage"] == 1.0 and out["name_grounded"] == 0.0
    assert math.isnan(out["stranded_poor"])                    # no node without annotation: FAIR-4 does not apply
    assert mb.score(out)["block_fair"] == 0.0                  # so the vacuous name earns nothing


def test_name_not_carried_by_any_member_is_not_traceable():
    g, m = flat(dict.fromkeys(["r1", "r2", "r3", "r4"], "g"),
                dict(g=dict(label="kinase", name_source="go", name_xref="GO:0016301")))
    out = run(FOUR, g, m)
    assert out["name_traceable"] == 0.0 and out["name_grounded"] == 0.0
    assert out["name_support"] == 0.0 and out["name_lift"] == 0.0     # carried by nobody: support 0, prevalence 0


def test_name_carried_only_outside_the_group_has_negative_lift():
    # "B" on the group of the two A reactions: support 0, prevalence 1/2 -> lift (0 - 1/2)/(1 - 1/2) = -1
    g, m = flat(dict(r1="g1", r2="g1", r3="g2", r4="g2"), dict(g1=dict(label="B", name_source="sbo", name_xref="B")))
    out = run(FOUR, g, m)
    assert out["name_lift"] == pytest.approx(-1.0) and out["name_traceable"] == 0.0
    assert mb.judge(out)["FAIR-2"] == "fail"


def test_member_name_counts_half_and_resolvable_name_counts_full():
    g, m = flat(dict(r1="g1", r2="g1", r3="g2", r4="g2"),
                dict(g1=dict(label="around r1", name_source="member", name_xref="r1"),
                     g2=dict(label="B", name_source="sbo", name_xref="SBO:0000002")))
    nodes = nodes_df([rx("r1", "SBO:0000000>SBO:0000001"), rx("r2", "SBO:0000000>SBO:0000001"),
                      rx("r3", "SBO:0000000>SBO:0000002"), rx("r4", "SBO:0000000>SBO:0000002")])
    out = run(nodes, g, m)
    assert out["name_resolvable"] == pytest.approx(0.5)
    assert out["name_grounded"] == pytest.approx(0.5 + 0.5 * 0.5)
    assert out["name_coverage"] == 1.0


def test_name_from_notes_code_prefix():
    nodes = nodes_df([sp("s1", notes='{"mapman": "10.5.1_ethylene.biosynthesis"}'),
                      sp("s2", notes='{"mapman": "10.5.2_ethylene.perception;50.2_EC"}'),
                      sp("s3", notes='{"mapman": "14.5_RNA"}'), sp("s4", notes='{"mapman": "14.5_RNA"}')])
    g, m = flat(dict(s1="g1", s2="g1", s3="g2", s4="g2"),
                dict(g1=dict(label="ethylene", name_source="notes", name_xref="mapman:10.5"),
                     g2=dict(label="RNA", name_source="notes", name_xref="mapman:14.5")))
    out = run(nodes, g, m)
    assert out["name_support"] == 1.0 and out["name_lift"] == 1.0
    # and the same codes drive the "class" homogeneity channel: IC(10.5) = IC(14.5) = ln(4/2)/ln(4) = 0.5
    assert out["hom_class"] == pytest.approx(0.4)


# ---- FAIR-4 nodes without annotation, FAIR-1 determinism
def test_stranded_nodes_without_annotation():
    nodes = nodes_df([rx("r1", "X>A"), rx("r2", "X>B"), sp("s1"), sp("s2"), sp("s3")])
    edges = edges_df([("r1", "s1"), ("r1", "s2"), ("r2", "s3")])
    g = flat(dict(r1="g1", s1="g1", s2="g2", r2="g2", s3="g3"))      # s1 with its reaction; s2, s3 apart
    out = run(nodes, *g, edges=edges)
    assert out["poor_share"] == pytest.approx(3 / 5)
    assert out["stranded_poor"] == pytest.approx(2 / 3)
    assert out["stranded_annotated"] == pytest.approx(0.5)           # r1 has s1; r2 has nobody
    # a declared copy of s3 next to its reaction r2 counts
    g2, m2 = frames(dict(r1=("g1", "", "", "R1"), s1=("g1", "", "", "R1"), s2=("g2", "", "", "R1"),
                         r2=("g2", "", "", "R1"), s3=("g3", "", "", "R1", "g2")), dict(g1={}, g2={}, g3={}))
    assert run(nodes, g2, m2, edges=edges)["stranded_poor"] == pytest.approx(1 / 3)


def test_determinism():
    g, m = flat(dict(r1="g1", r2="g1", r3="g2", r4="g2"))
    assert math.isnan(run(FOUR, g, m)["fair_deterministic"])
    assert run(FOUR, g, m, rerun=(g.copy(), m.copy()))["fair_deterministic"] == 1.0
    g2, m2 = flat(dict(r1="g1", r2="g2", r3="g2", r4="g2"))
    out = run(FOUR, g, m, rerun=(g2, m2))
    assert out["fair_deterministic"] == 0.0 and out["fair_rerun_same_rows"] == pytest.approx(0.75)
    assert mb.judge(out)["FAIR-1"] == "fail" and mb.score(out)["score"] == 0.0


def test_metrics_are_reproducible_and_row_order_free():
    g, m = function_first()
    a = run(CELL_NODES, g, m, edges=CELL_EDGES, comps=CELL)
    b = run(CELL_NODES.iloc[::-1].reset_index(drop=True), g.iloc[::-1].reset_index(drop=True),
            m.iloc[::-1].reset_index(drop=True), edges=CELL_EDGES.iloc[::-1].reset_index(drop=True), comps=CELL)
    for k in a:
        if k == "fair_sorted":
            continue
        assert (math.isnan(a[k]) and math.isnan(b[k])) or a[k] == pytest.approx(b[k]), k


def test_utility_goes_through_the_thresholds():
    rule = {r["id"]: r for r in mb.RULES}
    assert mb.unit(rule["BIO-2"], 0.0) == 0.0 and mb.unit(rule["BIO-2"], -0.3) == 0.0
    assert mb.unit(rule["BIO-2"], 0.05) == pytest.approx(0.25) and mb.unit(rule["BIO-2"], 0.21) == pytest.approx(0.75)
    assert mb.unit(rule["BIO-2"], 1.0) == 1.0
    assert mb.unit(rule["LOC-2"], 1.0) == 1.0 and mb.unit(rule["LOC-2"], 2.0) == pytest.approx(0.25)
    assert mb.unit(rule["LOC-2"], 1e9) == pytest.approx(0.0, abs=1e-6)
    assert mb.unit(rule["FAIR-4"], 0.0) == 1.0 and mb.unit(rule["FAIR-4"], 0.10) == pytest.approx(0.75)
    assert mb.unit(rule["FAIR-4"], 0.30) == pytest.approx(0.25) and mb.unit(rule["FAIR-4"], 1.0) == 0.0
    assert mb.unit(rule["LOC-1"], 0.9) == pytest.approx(0.25) and mb.unit(rule["LOC-1"], 1.0) == 1.0


# ============================================================ part 2: degenerate groupings on the panel models
def _find_panel():
    cands = [os.environ.get("GROUPING_PANEL", "")]
    for base in (pathlib.Path(__file__).resolve().parent, pathlib.Path.cwd()):
        for up in (base, *list(base.parents)[:3]):
            cands += [up / "panel" / "features", up / "bundle" / "panel" / "features"]
    for c in cands:
        if c and (pathlib.Path(c) / "e_coli_core" / "nodes.tsv").exists():
            return pathlib.Path(c)
    return None


PANEL = _find_panel()
MODELS = ["brenda_merged", "e_coli_core", "R-HSA-77289", "skm_pss", "skm_ja", "skm_sa"]
DEGENERATE = ["one_group", "singletons", "singletons_as_localisation", "random"]
needs_panel = pytest.mark.skipif(PANEL is None, reason="panel/features not found: set GROUPING_PANEL")


def _loc_key(r, boundary=True):
    b = boundary and r.loc_kind == "boundary"
    return f"loc_{r.loc_space}" + ("_boundary" if b else ""), loc_meta(r.loc_space, b)


def _neighbours(edges):
    nb = defaultdict(set)
    for rid, sid in zip(edges.reaction_id, edges.species_id):
        nb[rid].add(sid)
        nb[sid].add(rid)
    return nb


def _communities(node_ids, nb, seed=0):
    nx = pytest.importorskip("networkx")
    g = nx.Graph()
    g.add_nodes_from(sorted(node_ids))
    g.add_edges_from((a, b) for a in sorted(node_ids) for b in sorted(nb.get(a, ())) if b in node_ids and a < b)
    return [sorted(c) for c in sorted(nx.community.louvain_communities(g, seed=seed), key=lambda c: sorted(c)[0])]


def build(nodes, edges):
    """name -> (groups, meta): the degenerate groupings, the compartment baseline and three controls."""
    nb = _neighbours(edges)
    node = {r.node_id: r for r in nodes.itertuples(index=False)}
    deg = {n: len(nb.get(n, ())) for n in node}
    coms = _communities(set(node), nb)
    out = {}
    out["one_group"] = flat(dict.fromkeys(node, "G"))
    out["singletons"] = flat({n: "g_" + n for n in node}, {"g_" + n: dict(basis="topology") for n in node})
    out["singletons_as_localisation"] = flat({n: "g_" + n for n in node},
                                             {"g_" + n: loc_meta(r.loc_space) for n, r in node.items()})
    ids = sorted(node)
    random.Random(0).shuffle(ids)
    assign, i = {}, 0
    for k, com in enumerate(sorted(coms, key=len, reverse=True)):      # realistic sizes: those of the communities
        for n in ids[i:i + len(com)]:
            assign[n] = f"r{k:04d}"
        i += len(com)
    out["random"] = flat(assign, {g: dict(basis="topology") for g in set(assign.values())})
    assign, meta = {}, {}
    for n, r in node.items():
        assign[n], md = _loc_key(r)
        meta[assign[n]] = md
    out["compartments"] = flat(assign, meta)
    # control 1: compartment, then communities inside it, named after their hub
    rows, meta2, by = {}, dict(meta), defaultdict(list)
    for n, g in assign.items():
        by[g].append(n)
    for g, mem in sorted(by.items()):
        for k, com in enumerate(_communities(set(mem), nb)):
            hub = sorted(com, key=lambda n: (-deg[n], n))[0]
            meta2[f"{g}_c{k:03d}"] = dict(level="2", parent=g, basis="topology", label="around " + hub,
                                          name_source="member", name_xref=hub)
            rows.update({n: (g, f"{g}_c{k:03d}", "", "R2") for n in com})
    out["compartments>communities"] = frames(rows, meta2)
    # control 2: communities of the whole graph, then compartment boxes inside each
    rows, meta3 = {}, {}
    for k, com in enumerate(coms):
        meta3[f"c{k:03d}"] = dict(basis="topology")
        for n in com:
            g, md = _loc_key(node[n])
            meta3[f"c{k:03d}_{g}"] = dict(md, level="2", parent=f"c{k:03d}")
            rows[n] = (f"c{k:03d}", f"c{k:03d}_{g}", "", "R2")
    out["communities>compartments"] = frames(rows, meta3)
    # control 3: compartment, then the most specific annotation a reaction has (EC sub-subclass, else SBO term);
    # a species follows the commonest class of its reactions
    cls = {}
    for n, r in node.items():
        if r.kind == "reaction":
            ec = sorted(filter(None, r.ec.split("|")))
            cls[n] = ("ec", ".".join(ec[0].split(".")[:3])) if ec else ("sbo", r.sbo) if r.sbo else None
    for n, r in node.items():
        if r.kind == "species":
            c = Counter(cls[x] for x in nb.get(n, ()) if cls.get(x))
            cls[n] = sorted(c, key=lambda t: (-c[t], t))[0] if c else None
    rows, meta4 = {}, dict(meta)
    for n, g in assign.items():
        if cls.get(n):
            src, term = cls[n]
            sub = f"{g}_{src}_{term}"
            meta4[sub] = dict(level="2", parent=g, basis="reaction_type", label=term, name_source=src,
                              name_xref=("ec-code:" + term) if src == "ec" else term)
            rows[n] = (g, sub, "", "R2")
        else:
            rows[n] = (g, "", "", "R3")
    out["compartments>type"] = frames(rows, meta4)
    return out, coms


_CACHE = {}


def panel(model):
    if model not in _CACHE:
        f = PANEL / model
        nodes, edges, comps = (mb.read_tsv(f / x) for x in ("nodes.tsv", "edges.tsv", "compartments.tsv"))
        groupings, coms = build(nodes, edges)
        # a planted reference: the communities are the "curated pathways"; a species that takes part in reactions
        # of other communities carries their labels too (multi-label); a fixed 30 % of nodes is unlabelled; in the
        # skm models only species are labelled.
        home = {n: f"P{k:03d}" for k, com in enumerate(coms) for n in com}
        nb = _neighbours(edges)
        kind = dict(zip(nodes.node_id, nodes.kind))
        labelled = [n for i, n in enumerate(sorted(home)) if i % 10 >= 3
                    and (kind[n] == "species" or not model.startswith("skm"))]
        labs = {n: {home[n]} | ({home[x] for x in nb.get(n, ())} if kind[n] == "species" else set()) for n in labelled}
        reference = pd.DataFrame([(n, x) for n in sorted(labs) for x in sorted(labs[n])], columns=["node_id", "label"])
        # the grouping that equals the planted reference: communities, with a copy of each multi-label species
        rows = {n: (home[n].replace("P", "c"), "", "", "R1",
                    "|".join(sorted(x.replace("P", "c") for x in labs.get(n, set()) - {home[n]})))
                for n in sorted(home)}
        groupings["planted"] = frames(rows, {g: dict(basis="topology") for g in {v[0] for v in rows.values()}})
        res = {}
        for name, (g, m) in groupings.items():
            out = mb.metrics(nodes, edges, comps, g, m, reference=reference, rerun=(g.copy(), m.copy()))
            out.update(mb.score(out))
            res[name] = out
        _CACHE[model] = (nodes, edges, comps, groupings, res)
    return _CACHE[model]


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_synthetic_groupings_keep_the_contract(model):
    try:
        import sys
        sys.path.insert(0, str(PANEL.parent.parent / "common"))
        import contract
    except ImportError:
        pytest.skip("common/contract.py not found next to the panel")
    import tempfile
    groupings = panel(model)[3]
    for name, (g, m) in groupings.items():
        d = pathlib.Path(tempfile.mkdtemp())
        g.to_csv(d / "groups.tsv", sep="\t", index=False, lineterminator="\n")
        m.to_csv(d / "groups_meta.tsv", sep="\t", index=False, lineterminator="\n")
        assert contract.validate(PANEL / model, d) == [], name


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_degenerate_groupings_score_at_chance_on_biology(model):
    res = panel(model)[4]
    for name in ("one_group", "singletons", "singletons_as_localisation"):
        assert res[name]["hom"] == pytest.approx(0.0, abs=1e-9), name
        assert res[name]["ref_skill"] <= 1e-9, name
    assert abs(res["random"]["hom"]) < 0.05
    assert res["random"]["ref_skill"] < 0.05
    assert res["compartments"]["ref_skill"] <= 1e-9          # a baseline by construction
    for name in DEGENERATE:
        v = mb.judge(res[name])
        assert v["BIO-1"] == "fail" and v["BIO-2"] == "fail", name


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_planted_reference_is_recovered(model):
    res = panel(model)[4]
    assert res["planted"]["ref_f1"] == pytest.approx(1.0)
    assert res["planted"]["ref_skill"] == pytest.approx(1.0)
    assert mb.judge(res["planted"])["BIO-1"] == "pass"
    # a hierarchy that holds the reference partition at its first level is found there
    assert res["communities>compartments"]["ref_best_level"] == 1.0
    assert res["communities>compartments"]["ref_skill"] > 0.25


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_degenerate_groupings_fail_localisation_and_names(model):
    nodes, _e, _c, _g, res = panel(model)
    multi = nodes.loc_space.nunique() > 1
    one, single, sloc, comp = (res[k] for k in ("one_group", "singletons", "singletons_as_localisation", "compartments"))
    if not multi:
        assert all(math.isnan(one[k]) for k in ("loc_expressed", "loc_repeat", "loc_nesting", "loc_span_ok"))
    else:
        assert one["loc_expressed"] == 0.0 and single["loc_expressed"] == 0.0
        assert one["loc_span_ok"] == 0.0
        assert mb.judge(one)["LOC-1"] == "fail" and mb.judge(single)["LOC-1"] == "fail"
        assert sloc["loc_expressed"] == 1.0 and mb.judge(sloc)["LOC-2"] == "fail"     # the counter-rule
        assert sloc["loc_nesting"] < 1.0
        assert {v for k, v in mb.judge(comp).items() if k.startswith("LOC")} == {"pass"}
    for r in (one, single, sloc):
        assert r["name_grounded"] == 0.0
    assert single["stranded_poor"] == 1.0 and one["stranded_poor"] == 0.0


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_no_degenerate_grouping_comes_out_best(model):
    res = panel(model)[4]
    controls = ["compartments>communities", "communities>compartments", "compartments>type", "planted"]
    worst_control = min(res[c]["score"] for c in controls)
    for name in DEGENERATE:
        assert res[name]["score"] < 0.25, (name, res[name]["score"])
        assert res[name]["score"] < worst_control, (name, res[name]["score"], worst_control)
    # the compartment baseline never beats every control
    assert res["compartments"]["score"] < max(res[c]["score"] for c in controls)


@needs_panel
@pytest.mark.parametrize("model", MODELS)
def test_panel_no_degenerate_grouping_comes_out_best_without_reference(model):
    nodes, edges, comps, groupings, _res = panel(model)
    sc = {}
    for name, (g, m) in groupings.items():
        if name == "planted":
            continue
        sc[name] = mb.score(mb.metrics(nodes, edges, comps, g, m))["score"]
    controls = ["compartments>communities", "communities>compartments", "compartments>type"]
    for name in DEGENERATE:
        assert sc[name] < 0.25 and sc[name] < min(sc[c] for c in controls), (name, sc)
    assert sc["compartments"] <= max(sc[c] for c in controls) + 1e-9
