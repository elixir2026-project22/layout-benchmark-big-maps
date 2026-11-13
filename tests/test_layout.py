"""The layout stage: the options, the contract, the recursive placement, determinism.

These tests build their own four-species model, so they run in a clean clone,
where ``panel/models/`` and ``panel/features/`` are not present.
"""
from __future__ import annotations

import csv

import pytest

EXPECTED = [
    ("graphviz-dot", "layered"),
    ("igraph-sugiyama", "layered (stands in for tree)"),
    ("graphviz-neato", "stress"),
    ("graphviz-sfdp", "force-directed"),
    ("igraph-drl", "for large graphs"),
    ("graphviz-twopi", "radial"),
    ("graphviz-circo", "circular"),
    ("igraph-mds", "baseline"),
]


def test_one_option_per_family():
    import options as OPT

    assert [(n, f) for n, (_fn, f, _link) in OPT.OPTIONS.items()] == EXPECTED
    assert len({f for _n, f in EXPECTED}) == len(EXPECTED), "one option per family"


def test_every_option_is_in_the_committed_gallery(root):
    import options as OPT

    with open(root / "layout" / "examples" / "gallery.tsv", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 288, "8 options x 9 groupings x 4 models"
    assert {r["option"] for r in rows} == set(OPT.OPTIONS)
    assert not [r for r in rows if r["error"]], "the committed gallery has no failures"


def test_the_grouping_is_contract_valid(tiny_features, tiny_grouping):
    import contract

    assert contract.validate(str(tiny_features), str(tiny_grouping)) == []


@pytest.mark.parametrize("option", [n for n, _f in EXPECTED])
def test_placement_places_every_node(option, tiny_features, tiny_grouping):
    import options as OPT
    import placement as PL

    if not OPT.available(option):
        pytest.skip(f"{option}: its tool is not on this machine")
    place = PL.Placement(str(tiny_features), str(tiny_grouping), option, seed=0)
    pos = place.positions_rows()
    with open(tiny_features / "nodes.tsv", encoding="utf-8", newline="") as fh:
        n_nodes = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
    assert len({r["node_id"] for r in pos}) == n_nodes
    assert all(float(r["width"]) > 0 and float(r["height"]) > 0 for r in pos)


def test_sibling_boxes_do_not_overlap(tiny_features, tiny_grouping):
    import placement as PL

    place = PL.Placement(str(tiny_features), str(tiny_grouping), "igraph-mds", seed=0)
    boxes = place.boxes_rows()
    by_parent = {}
    for b in boxes:
        by_parent.setdefault(b["parent"], []).append(b)
    for siblings in by_parent.values():
        for i, a in enumerate(siblings):
            for b in siblings[i + 1:]:
                ax, ay, aw, ah = (float(a[k]) for k in ("x", "y", "width", "height"))
                bx, by, bw, bh = (float(b[k]) for k in ("x", "y", "width", "height"))
                overlap = max(0.0, min(ax + aw, bx + bw) - max(ax, bx)) * \
                          max(0.0, min(ay + ah, by + bh) - max(ay, by))
                assert overlap == 0.0, f"{a['box_id']} overlaps {b['box_id']}"


def test_one_seed_gives_one_placement(tiny_features, tiny_grouping):
    """The determinism contract: the same seed twice gives the same tables."""
    import placement as PL

    first = PL.Placement(str(tiny_features), str(tiny_grouping), "igraph-mds", seed=0)
    second = PL.Placement(str(tiny_features), str(tiny_grouping), "igraph-mds", seed=0)
    assert first.positions_rows() == second.positions_rows()
    assert first.boxes_rows() == second.boxes_rows()


def test_a_group_of_one_entity_is_not_a_group(tiny_features, tmp_path):
    """Rule 2 of PROTOCOL.md: the box is dropped, its node goes to the parent."""
    import placement as PL

    with open(tiny_features / "nodes.tsv", encoding="utf-8", newline="") as fh:
        nodes = sorted(r["node_id"] for r in csv.DictReader(fh, delimiter="\t"))
    out = tmp_path / "singletons"
    out.mkdir()
    with open(out / "groups.tsv", "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["node_id", "L1", "L2", "L3", "rule", "clone_in"],
                            delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for i, nid in enumerate(nodes):
            wr.writerow({"node_id": nid, "L1": f"G{i}", "L2": "", "L3": "", "rule": "R1", "clone_in": ""})
    with open(out / "groups_meta.tsv", "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["group_id", "level", "parent", "label", "basis",
                                            "name_source", "name_xref", "compartment_id"],
                            delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for i in range(len(nodes)):
            wr.writerow({"group_id": f"G{i}", "level": "1", "parent": "", "label": "", "basis": "mixed",
                         "name_source": "none", "name_xref": "", "compartment_id": ""})
    place = PL.Placement(str(tiny_features), str(out), "igraph-mds", seed=0)
    assert len(place.dropped_singletons) == len(nodes)
    assert [b for b in place.boxes_rows() if b["kind"] == "group"] == []
