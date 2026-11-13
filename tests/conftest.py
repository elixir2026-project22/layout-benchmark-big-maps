"""Test setup: the repository root, and the repository's own scripts on sys.path.

The working code is run as scripts (``python layout/placement.py ...``), so the
tests import it the way the scripts import each other.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
for rel in ("src", "common", "examples", "layout", "rules", "scoring"):
    p = str(ROOT / rel)
    if p not in sys.path:
        sys.path.insert(0, p)


@pytest.fixture(scope="session")
def root() -> pathlib.Path:
    """The repository checkout under test."""
    return ROOT


@pytest.fixture(scope="session")
def tiny_model(tmp_path_factory) -> pathlib.Path:
    """A four-species, two-reaction SBML model written from scratch.

    Nothing in ``panel/`` is needed for it, so the layout stage can be tested in
    a clean clone, where the models themselves are not present.
    """
    libsbml = pytest.importorskip("libsbml")
    doc = libsbml.SBMLDocument(3, 2)
    m = doc.createModel()
    m.setId("tiny")
    for cid, name in (("c", "cytosol"), ("e", "extracellular")):
        c = m.createCompartment()
        c.setId(cid)
        c.setName(name)
        c.setConstant(True)
        c.setSize(1.0)
    for sid, comp in (("a", "c"), ("b", "c"), ("d", "c"), ("x", "e")):
        s = m.createSpecies()
        s.setId(sid)
        s.setName(sid.upper())
        s.setCompartment(comp)
        for setter in (s.setConstant, s.setBoundaryCondition, s.setHasOnlySubstanceUnits):
            setter(False)
    for rid, reactants, products in (("r1", ["a"], ["b"]), ("r2", ["b"], ["d"]), ("r3", ["x"], ["a"])):
        r = m.createReaction()
        r.setId(rid)
        r.setReversible(False)
        r.setFast(False)
        for sid in reactants:
            sr = r.createReactant()
            sr.setSpecies(sid)
            sr.setConstant(True)
            sr.setStoichiometry(1.0)
        for sid in products:
            sr = r.createProduct()
            sr.setSpecies(sid)
            sr.setConstant(True)
            sr.setStoichiometry(1.0)
    path = tmp_path_factory.mktemp("tiny") / "tiny.xml"
    libsbml.writeSBMLToFile(doc, str(path))
    return path


@pytest.fixture(scope="session")
def tiny_features(tiny_model, tmp_path_factory) -> pathlib.Path:
    """``panel/features/<model>``-style tables for :func:`tiny_model`."""
    import sbml_features as F

    out = tmp_path_factory.mktemp("features")
    nodes, edges, comps = F.extract(str(tiny_model))
    F.write_tsv(out / "nodes.tsv", nodes, F.NODE_COLUMNS)
    F.write_tsv(out / "edges.tsv", edges, F.EDGE_COLUMNS)
    F.write_tsv(out / "compartments.tsv", comps, F.COMPARTMENT_COLUMNS)
    return out


@pytest.fixture(scope="session")
def tiny_grouping(tiny_features, tmp_path_factory) -> pathlib.Path:
    """A two-group grouping of :func:`tiny_features`, in the contract's format."""
    import csv

    with open(tiny_features / "nodes.tsv", encoding="utf-8", newline="") as fh:
        nodes = list(csv.DictReader(fh, delimiter="\t"))
    out = tmp_path_factory.mktemp("grouping")
    with open(out / "groups.tsv", "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["node_id", "L1", "L2", "L3", "rule", "clone_in"],
                            delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for n in sorted(nodes, key=lambda n: n["node_id"]):
            wr.writerow({"node_id": n["node_id"], "L1": "G1" if n["kind"] == "species" else "G2",
                         "L2": "", "L3": "", "rule": "R1", "clone_in": ""})
    with open(out / "groups_meta.tsv", "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=["group_id", "level", "parent", "label", "basis",
                                            "name_source", "name_xref", "compartment_id"],
                            delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for gid in ("G1", "G2"):
            wr.writerow({"group_id": gid, "level": "1", "parent": "", "label": "", "basis": "mixed",
                         "name_source": "none", "name_xref": "", "compartment_id": ""})
    return out
