"""SBML -> node, edge and compartment tables: the shared input of every grouping strategy.

Every strategy starts from these tables, so strategies differ in grouping logic only, never in how
they read the model. Nothing here groups anything.

Run:
    python common/sbml_features.py MODEL.xml -o panel/features/MODEL

Writes into the output folder:
    nodes.tsv         one row per species and per reaction
    edges.tsv         one row per reaction participant (reactant, product, modifier)
    compartments.tsv  the compartment containment tree (see compartments.py)

nodes.tsv columns
    node_id, kind (species|reaction), name
    compartment       species: its compartment. reaction: the SBML `compartment` attribute, else the one
                      compartment of its reactants and products, else empty
    compartments      "|"-joined compartments of the reactants and products (species: its own)
    modifier_compartments   reaction only: compartments of its modifiers
    loc_space, loc_kind, loc_outer   where the node is drawn in the cell (compartments.Tree.locate):
                      inside a space, on the boundary of a space (transport across it, or a membrane
                      participant; loc_outer is the space on the other side), or between sibling spaces
                      (loc_space is then their common parent). Computed here once, identically for all.
    sbo, sbo_name, sbo_path, sbo_path_names   SBO term, and its ancestors root > ... > term (first parent)
    ec, go            "|"-joined EC numbers and GO ids from the annotation
    xrefs             "|"-joined "qualifier=db:id" for every annotation resource
    notes             JSON object of the "key: value" lines in the SBML notes ("text" holds free text)
    reversible, n_reactants, n_products, n_modifiers   reaction only
    degree            species: number of reactions it takes part in. reaction: number of participants
    boundary          species only: boundaryCondition

edges.tsv columns
    reaction_id, species_id, role (reactant|product|modifier), stoichiometry, sbo, sbo_name
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import re
import sys

import libsbml

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compartments import COLUMNS as COMPARTMENT_COLUMNS  # noqa: E402
from compartments import Tree, build_tree  # noqa: E402

NODE_COLUMNS = ["node_id", "kind", "name", "compartment", "compartments", "modifier_compartments", "loc_space",
                "loc_kind", "loc_outer", "sbo", "sbo_name", "sbo_path", "sbo_path_names", "ec", "go", "xrefs", "notes",
                "reversible", "n_reactants", "n_products", "n_modifiers", "degree", "boundary"]
EDGE_COLUMNS = ["reaction_id", "species_id", "role", "stoichiometry", "sbo", "sbo_name"]
_SBO = None


def sbo_terms() -> dict[str, dict]:
    """SBO from the OBO file next to this script: id -> {"name", "is_a": [parent ids]}."""
    global _SBO
    if _SBO is None:
        terms, cur = {}, None
        for line in (HERE / "SBO_OBO.obo").read_text(encoding="utf-8", errors="replace").splitlines() + [""]:
            line = line.strip()
            if line == "[Term]":
                cur = {"is_a": []}
            elif not line:
                if cur and "id" in cur:
                    terms[cur["id"]] = cur
                cur = None
            elif cur is not None:
                if line.startswith("id: "):
                    cur["id"] = line[4:]
                elif line.startswith("name: "):
                    cur["name"] = line[6:]
                elif line.startswith("is_a: "):
                    cur["is_a"].append(line[6:].split("!")[0].strip())
        _SBO = terms
    return _SBO


def sbo_path(term: str) -> list[str]:
    """Root > ... > term, following the first parent (SBO terms used for reactions have one)."""
    terms = sbo_terms()
    if term not in terms:
        return [term] if term else []
    out = [term]
    while terms.get(out[-1], {}).get("is_a") and terms[out[-1]]["is_a"][0] not in out:
        out.append(terms[out[-1]]["is_a"][0])
    return out[::-1]


def sbo_name(term: str) -> str:
    return sbo_terms().get(term, {}).get("name", "")


def _resource(uri: str) -> tuple[str, str]:
    """An annotation resource as (database, identifier)."""
    u = uri.strip()
    m = re.search(r"chebiId=(CHEBI:\d+)", u)
    if m:
        return "chebi", m.group(1)
    m = re.search(r"rhea-db\.org/rhea/(\d+)", u)
    if m:
        return "rhea", m.group(1)
    m = re.match(r"urn:miriam:([^:]+):(.+)$", u)
    if m:
        return m.group(1).lower(), m.group(2)
    m = re.search(r"identifiers\.org/(.+)$", u)
    if m:
        rest = m.group(1)
        if "/" in rest and ":" not in rest.split("/")[0]:
            db, ident = rest.split("/", 1)
            return db.lower(), ident
        if ":" in rest:
            db, ident = rest.split(":", 1)
            return db.lower(), (rest if db.upper() in ("GO", "CHEBI") else ident)
    m = re.match(r"([A-Za-z][\w.]*):(.+)$", u)
    if m and not u.lower().startswith("http"):
        return m.group(1).lower(), m.group(2)
    return "uri", u


def _xrefs(el) -> list[tuple[str, str, str]]:
    out = []
    for i in range(el.getNumCVTerms()):
        cv = el.getCVTerm(i)
        if cv.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER:
            qual = "bqbiol:" + libsbml.BiolQualifierType_toString(cv.getBiologicalQualifierType())
        else:
            qual = "bqmodel:" + libsbml.ModelQualifierType_toString(cv.getModelQualifierType())
        for j in range(cv.getNumResources()):
            db, ident = _resource(cv.getResourceURI(j))
            out.append((qual, db, ident))
    return out


def _notes(el) -> dict[str, str]:
    raw = el.getNotesString() or ""
    if not raw:
        return {}
    out, free = {}, []
    for chunk in re.split(r"</p>|<br\s*/?>|\n", raw):
        text = re.sub(r"<[^>]+>", " ", chunk)
        text = re.sub(r"\s+", " ", text).strip()
        if not text:
            continue
        m = re.match(r"([A-Za-z][A-Za-z0-9_ ]{1,40}):\s*(.*)$", text)
        if m:
            out[m.group(1).strip()] = m.group(2).strip()[:1000]
        else:
            free.append(text)
    if free:
        out["text"] = " ".join(free)[:1000]
    return out


def _join(values) -> str:
    return "|".join(sorted({v for v in values if v}))


def extract(path) -> tuple[list[dict], list[dict], list[dict]]:
    """(node rows, edge rows, compartment rows) of one SBML file."""
    doc = libsbml.readSBMLFromFile(str(path))
    model = doc.getModel()
    if model is None:
        raise RuntimeError(f"{path}: not readable as SBML")

    species = [model.getSpecies(i) for i in range(model.getNumSpecies())]
    reactions = [model.getReaction(i) for i in range(model.getNumReactions())]
    sp_comp = {s.getId(): s.getCompartment() for s in species}

    edges, degree = [], collections.Counter()
    r_comps, r_mod_comps = {}, {}
    for r in reactions:
        parts = ([("reactant", r.getReactant(i)) for i in range(r.getNumReactants())]
                 + [("product", r.getProduct(i)) for i in range(r.getNumProducts())]
                 + [("modifier", r.getModifier(i)) for i in range(r.getNumModifiers())])
        seen = set()
        for role, ref in parts:
            sid = ref.getSpecies()
            stoich = ref.getStoichiometry() if role != "modifier" and ref.isSetStoichiometry() else ""
            edges.append(dict(reaction_id=r.getId(), species_id=sid, role=role, stoichiometry=stoich,
                              sbo=ref.getSBOTermID(), sbo_name=sbo_name(ref.getSBOTermID())))
            seen.add(sid)
        for sid in seen:
            degree[sid] += 1
        degree[r.getId()] = len(seen)
        r_comps[r.getId()] = {sp_comp.get(ref.getSpecies(), "") for role, ref in parts if role != "modifier"} - {""}
        r_mod_comps[r.getId()] = {sp_comp.get(ref.getSpecies(), "") for role, ref in parts if role == "modifier"} - {""}

    # compartment containment tree
    pairs = collections.Counter()
    for rid, comps in r_comps.items():
        allc = sorted(comps | r_mod_comps[rid])
        for i, a in enumerate(allc):
            for b in allc[i + 1:]:
                pairs[frozenset((a, b))] += 1
    n_in = collections.Counter(sp_comp.values())
    comp_in = []
    for i in range(model.getNumCompartments()):
        c = model.getCompartment(i)
        gos = [ident for _q, db, ident in _xrefs(c) if db == "go"]
        outside = c.getOutside() if hasattr(c, "isSetOutside") and c.isSetOutside() else ""
        comp_in.append(dict(id=c.getId(), name=c.getName() or c.getId(), go=gos, outside=outside,
                            n_species=n_in.get(c.getId(), 0)))
    comp_rows = build_tree(comp_in, pairs)
    tree = Tree(comp_rows)

    def common(el, kind):
        x = _xrefs(el)
        term = el.getSBOTermID()
        path_ = sbo_path(term)
        return dict(node_id=el.getId(), kind=kind, name=el.getName() or "", sbo=term, sbo_name=sbo_name(term),
                    sbo_path=">".join(path_), sbo_path_names=" > ".join(sbo_name(t) or t for t in path_),
                    ec=_join(ident for _q, db, ident in x if db in ("ec-code", "ec")),
                    go=_join(ident for _q, db, ident in x if db == "go"),
                    xrefs="|".join(sorted({f"{q}={db}:{ident}" for q, db, ident in x})),
                    notes=json.dumps(_notes(el), ensure_ascii=False, sort_keys=True) if _notes(el) else "",
                    degree=degree.get(el.getId(), 0))

    nodes = []
    for s in species:
        row = common(s, "species")
        space, kind, outer = tree.locate([s.getCompartment()])
        row.update(compartment=s.getCompartment(), compartments=s.getCompartment(), modifier_compartments="",
                   loc_space=space, loc_kind=kind, loc_outer=outer, reversible="", n_reactants="", n_products="",
                   n_modifiers="", boundary=int(s.getBoundaryCondition()))
        nodes.append(row)
    for r in reactions:
        row = common(r, "reaction")
        comps = r_comps[r.getId()]
        own = r.getCompartment() if r.isSetCompartment() else (next(iter(comps)) if len(comps) == 1 else "")
        where = comps or ({own} if own else set()) or r_mod_comps[r.getId()]
        space, kind, outer = tree.locate(where)
        row.update(compartment=own, compartments=_join(comps), modifier_compartments=_join(r_mod_comps[r.getId()]),
                   loc_space=space, loc_kind=kind, loc_outer=outer, reversible=int(r.getReversible()),
                   n_reactants=r.getNumReactants(), n_products=r.getNumProducts(), n_modifiers=r.getNumModifiers(),
                   boundary="")
        nodes.append(row)
    return nodes, edges, comp_rows


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("sbml")
    ap.add_argument("-o", "--out", required=True, help="output folder")
    a = ap.parse_args(argv)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    nodes, edges, comps = extract(a.sbml)
    write_tsv(out / "nodes.tsv", nodes, NODE_COLUMNS)
    write_tsv(out / "edges.tsv", edges, EDGE_COLUMNS)
    write_tsv(out / "compartments.tsv", comps, COMPARTMENT_COLUMNS)
    n_sp = sum(r["kind"] == "species" for r in nodes)
    print(f"{pathlib.Path(a.sbml).name}: {n_sp} species, {len(nodes) - n_sp} reactions, {len(edges)} edges, "
          f"{len(comps)} compartments -> {out}")


if __name__ == "__main__":
    main()
