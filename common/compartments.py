"""Compartment containment tree: which compartment lies inside which.

Localisation must match cell biology: the extracellular region surrounds the cytosol, organelles lie
inside the cytosol, lumens inside their organelle, and a membrane is the boundary of the space it
encloses. SBML Level 3 has no containment attribute, and GO's own part_of does not place the nucleus
inside the cytoplasm, so the tree comes from an explicit vocabulary of cell compartments (VOCAB
below; GO ids in go_cellular_component.json were retrieved from OLS, not typed by hand).

Each model compartment is matched to a vocabulary entry, in this fixed order (recorded in `source`):
  1. sbml_outside        SBML Level 2 `outside` attribute, when a model has it
  2. go_annotation       the compartment's GO cellular-component cross-reference
  3. name_vocabulary     the compartment name; the longest matching pattern wins
     name_vocabulary_duplicate  a second compartment matching an entry already taken by a larger
                         one is a sub-compartment of it (e.g. "Inner mitochondria" inside "Mitochondria")
  4. transport_adjacency no match: parent = the space it exchanges most reactions with
  5. single_compartment  a model with one unmatched compartment: it is the root

Vocabulary ancestors that the model lacks are added as *virtual* spaces only where needed to keep the
biology (e.g. a model with nucleus and chloroplast but no cytosol gets a virtual cytosol). A virtual
space with a single child is collapsed into that child.

Run (normally called through sbml_features.py):
    python common/compartments.py MODEL.xml
"""
from __future__ import annotations

import collections
import json
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent

# key: (GO labels, parent key, optional, name pattern)
# optional = never created as a virtual space (skipped in the ancestor chain when the model lacks it)
SPACES = {
    "extracellular": (["extracellular region", "extracellular space"], None, False,
                      r"extra[- ]?cellular|external|^medium$|environment|^outside$|exterior"),
    "apoplast": (["apoplast"], "extracellular", True, r"apoplast"),
    "cell_wall": (["cell wall"], "apoplast", True, r"cell wall"),
    "periplasm": (["periplasmic space"], "cell_wall", True, r"periplasm"),
    "cytosol": (["cytosol", "cytoplasm"], "periplasm", False, r"cytosol|cytoplasm|^intracellular$|^cell$"),
    "nucleus": (["nucleus"], "cytosol", False, r"nucleus|^nuclear$"),
    "nucleoplasm": (["nucleoplasm"], "nucleus", False, r"nucleoplasm"),
    "mitochondrion": (["mitochondrion"], "cytosol", False, r"mitochondri"),
    "mito_intermembrane": (["mitochondrial intermembrane space"], "mitochondrion", True, r"intermembrane"),
    "mito_matrix": (["mitochondrial matrix"], "mito_intermembrane", False, r"mitochondri\w* matrix"),
    "chloroplast": (["chloroplast", "plastid"], "cytosol", False, r"chloroplast|plastid"),
    "chloroplast_stroma": (["chloroplast stroma"], "chloroplast", True, r"(chloroplast|plastid)\w* stroma|^stroma$"),
    "thylakoid": (["thylakoid"], "chloroplast_stroma", False, r"thylakoid"),
    "thylakoid_lumen": (["thylakoid lumen"], "thylakoid", False, r"thylakoid lumen"),
    "er": (["endoplasmic reticulum"], "cytosol", False, r"endoplasmic reticulum|^er$"),
    "er_lumen": (["endoplasmic reticulum lumen"], "er", False, r"(endoplasmic reticulum|\ber) lumen"),
    "golgi": (["Golgi apparatus"], "cytosol", False, r"golgi"),
    "golgi_lumen": (["Golgi lumen"], "golgi", False, r"golgi\w* lumen"),
    "peroxisome": (["peroxisome"], "cytosol", False, r"peroxisom|glyoxysom"),
    "peroxisomal_matrix": (["peroxisomal matrix"], "peroxisome", False, r"peroxisom\w* matrix"),
    "vacuole": (["vacuole"], "cytosol", False, r"vacuol"),
    "vacuolar_lumen": (["vacuolar lumen"], "vacuole", False, r"vacuol\w* lumen"),
    "lysosome": (["lysosome"], "cytosol", False, r"lysosom"),
    "lysosomal_lumen": (["lysosomal lumen"], "lysosome", False, r"lysosom\w* lumen"),
    "endosome": (["endosome"], "cytosol", False, r"endosom"),
    "endosome_lumen": (["endosome lumen"], "endosome", False, r"endosom\w* lumen"),
}
# key: (GO labels, the space it encloses, name pattern)
MEMBRANES = {
    "plasma_membrane": (["plasma membrane"], "cytosol", r"plasma membrane|cell membrane|plasmalemma|cytoplasmic membrane"),
    "nuclear_envelope": (["nuclear envelope"], "nucleus", r"nuclear (envelope|membrane)"),
    "mito_outer_membrane": (["mitochondrial outer membrane"], "mitochondrion",
                            r"mitochondri\w* outer membrane|outer mitochondri\w* membrane|mitochondri\w* (membrane|envelope)"),
    "mito_inner_membrane": (["mitochondrial inner membrane"], "mito_matrix",
                            r"mitochondri\w* inner membrane|inner mitochondri\w* membrane"),
    "er_membrane": (["endoplasmic reticulum membrane"], "er", r"(endoplasmic reticulum|\ber) membrane"),
    "golgi_membrane": (["Golgi membrane"], "golgi", r"golgi\w* membrane"),
    "peroxisomal_membrane": (["peroxisomal membrane"], "peroxisome", r"peroxisom\w* membrane"),
    "vacuolar_membrane": (["vacuolar membrane"], "vacuole", r"vacuol\w* membrane|tonoplast"),
    "lysosomal_membrane": (["lysosomal membrane"], "lysosome", r"lysosom\w* membrane"),
    "chloroplast_envelope": (["chloroplast envelope"], "chloroplast", r"(chloroplast|plastid)\w* (envelope|membrane)"),
    "thylakoid_membrane": (["thylakoid membrane"], "thylakoid", r"thylakoid membrane"),
    "endosome_membrane": (["endosome membrane"], "endosome", r"endosom\w* membrane"),
}
COLUMNS = ["compartment_id", "name", "canonical", "go_id", "kind", "parent_id", "encloses_id", "virtual",
           "depth", "source", "n_species"]


def _go_table() -> dict[str, str]:
    """GO label (lower case) -> GO id, from the OLS lookup stored next to this file."""
    path = HERE / "go_cellular_component.json"
    if not path.exists():
        return {}
    return {k.lower(): v for k, v in json.loads(path.read_text(encoding="utf-8")).items()}


def _vocab_by_go() -> dict[str, str]:
    go = _go_table()
    out = {}
    for table in (SPACES, MEMBRANES):
        for key, spec in table.items():
            for label in spec[0]:
                if label.lower() in go:
                    out[go[label.lower()]] = key
    return out


def _key_go_id(key: str) -> str:
    go = _go_table()
    labels = (SPACES.get(key) or MEMBRANES.get(key))[0]
    return next((go[label.lower()] for label in labels if label.lower() in go), "")


def match_name(name: str) -> str | None:
    """Vocabulary key for a compartment name: the longest matching pattern wins."""
    text = re.sub(r"[_\-]+", " ", (name or "").strip().lower())
    best, best_len = None, 0
    for table, idx in ((SPACES, 3), (MEMBRANES, 2)):
        for key, spec in table.items():
            m = re.search(spec[idx], text)
            if m and len(m.group(0)) > best_len:
                best, best_len = key, len(m.group(0))
    return best


def _space_chain(key: str) -> list[str]:
    out = []
    while key is not None:
        out.append(key)
        key = SPACES[key][1]
    return out  # key, parent, ..., root


def build_tree(compartments: list[dict], transport_pairs: collections.Counter) -> list[dict]:
    """Rows of compartments.tsv.

    compartments: [{"id", "name", "go": [GO ids], "outside": id or "", "n_species"}]
    transport_pairs: Counter of frozenset({compartment a, compartment b}) over reactions touching both.
    """
    by_go = _vocab_by_go()
    key_of, source = {}, {}
    taken, inside_first = {}, {}
    order = sorted(compartments, key=lambda c: (-c.get("n_species", 0), c["id"]))
    for c in order:
        key, src = None, None
        for g in c.get("go", []):
            if g in by_go:
                key, src = by_go[g], "go_annotation"
                break
        if key is None:
            key = match_name(c.get("name") or "")
            src = "name_vocabulary" if key else None
        if key is not None and key in taken:  # two compartments, one vocabulary entry: the smaller lies in the larger
            inside_first[c["id"]] = taken[key]
            key, src = None, None
        if key is not None:
            taken[key] = c["id"]
        key_of[c["id"]], source[c["id"]] = key, src

    # tree over vocabulary keys: present keys plus the ancestors needed
    present_spaces = {k for k in key_of.values() if k in SPACES}
    present_membranes = {k for k in key_of.values() if k in MEMBRANES}
    needed = set(present_spaces)
    anchored = set()  # virtual spaces that a present membrane encloses
    for mk in present_membranes:
        needed.add(MEMBRANES[mk][1])
        anchored.add(MEMBRANES[mk][1])
    nodes = set()
    for k in needed:
        for a in _space_chain(k):
            if a in present_spaces or a in needed or not SPACES[a][2]:
                nodes.add(a)

    def parent_key(k):
        for a in _space_chain(k)[1:]:
            if a in nodes:
                return a
        return None

    parent = {k: parent_key(k) for k in nodes}

    def children(k):
        return sorted(c for c, p in parent.items() if p == k)

    def is_virtual(k):
        return k not in present_spaces

    # drop virtual roots above the model, then collapse virtual spaces with one child
    changed = True
    while changed:
        changed = False
        for k in sorted(nodes):
            if not is_virtual(k) or k in anchored:
                continue
            kids = children(k)
            if parent[k] is None and len(kids) <= 1 or parent[k] is not None and len(kids) == 1 or not kids:
                for c in kids:
                    parent[c] = parent[k]
                nodes.discard(k)
                parent.pop(k)
                changed = True
                break
    roots = sorted(k for k in nodes if parent[k] is None)
    if len(roots) > 1:  # cannot happen with one vocabulary root, kept as a guard
        raise RuntimeError(f"more than one root: {roots}")

    node_id = {k: (taken[k] if k in taken else f"virtual_{k}") for k in nodes}
    for mk in present_membranes:
        node_id[mk] = taken[mk]
    rows = {}
    info = {c["id"]: c for c in compartments}
    for k in sorted(nodes):
        cid = node_id[k]
        rows[cid] = dict(compartment_id=cid, name=info[cid]["name"] if cid in info else SPACES[k][0][0], canonical=k,
                         go_id=_key_go_id(k), kind="space", parent_id=node_id[parent[k]] if parent[k] else "",
                         encloses_id="", virtual=int(is_virtual(k)),
                         source=source.get(cid) or "vocabulary_ancestor", n_species=info.get(cid, {}).get("n_species", 0))
    for mk in sorted(present_membranes):
        cid = node_id[mk]
        inner = MEMBRANES[mk][1]
        while inner not in nodes:  # the enclosed space was collapsed: the membrane bounds what replaced it
            kids = [c for c in SPACES if SPACES[c][1] == inner and (c in nodes or any(x in nodes for x in _descend(c)))]
            if not kids:
                break
            inner = kids[0] if kids[0] in nodes else next(x for x in _descend(kids[0]) if x in nodes)
        inner_id = node_id.get(inner, "")
        rows[cid] = dict(compartment_id=cid, name=info[cid]["name"], canonical=mk, go_id=_key_go_id(mk), kind="membrane",
                         parent_id=rows[inner_id]["parent_id"] if inner_id in rows else "", encloses_id=inner_id,
                         virtual=0, source=source[cid], n_species=info[cid].get("n_species", 0))

    # compartments the vocabulary does not know
    unmatched = [c["id"] for c in order if key_of[c["id"]] is None]
    spaces_now = [cid for cid, r in rows.items() if r["kind"] == "space"]
    root_id = next((cid for cid in spaces_now if rows[cid]["parent_id"] == ""), "")
    cyt_id = next((cid for cid in spaces_now if rows[cid]["canonical"] == "cytosol"), "")
    for cid in unmatched:
        c = info[cid]
        if c.get("outside") and c["outside"] in info:
            par, src = c["outside"], "sbml_outside"
        elif cid in inside_first:
            first = rows[inside_first[cid]]
            par = first["encloses_id"] if first["kind"] == "membrane" else first["compartment_id"]
            src = "name_vocabulary_duplicate"
        elif not rows:
            par, src = "", "single_compartment" if len(compartments) == 1 else "unmatched_root"
        else:
            cand = collections.Counter({o: n for pair, n in transport_pairs.items() if cid in pair
                                        for o in pair if o != cid and o in rows and rows[o]["kind"] == "space"})
            if cand:
                top = max(cand.values())
                par = sorted(o for o, n in cand.items() if n == top)[0]
                if cyt_id and cand.get(cyt_id) == top:
                    par = cyt_id
                src = "transport_adjacency"
            else:
                par, src = (cyt_id or root_id), "unmatched_default"
        rows[cid] = dict(compartment_id=cid, name=c["name"], canonical="", go_id="", kind="space", parent_id=par,
                         encloses_id="", virtual=0, source=src, n_species=c.get("n_species", 0))

    def depth(cid, seen=()):
        p = rows[cid]["parent_id"]
        return 0 if not p or p not in rows or cid in seen else 1 + depth(p, (*seen, cid))

    for cid in rows:
        rows[cid]["depth"] = depth(cid)
    return sorted(rows.values(), key=lambda r: (r["depth"], r["kind"], r["compartment_id"]))


def _descend(key):
    out = []
    for c in SPACES:
        if SPACES[c][1] == key:
            out.append(c)
            out += _descend(c)
    return out


class Tree:
    """Queries on compartments.tsv rows (list of dicts, or a path to the file)."""

    def __init__(self, rows):
        if isinstance(rows, (str, pathlib.Path)):
            import csv
            with open(rows, encoding="utf-8", newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        self.rows = {r["compartment_id"]: r for r in rows}

    def space_of(self, cid: str) -> tuple[str, bool]:
        """(space id, on its boundary?) - a membrane is the boundary of the space it encloses."""
        r = self.rows[cid]
        if r["kind"] == "membrane" and r["encloses_id"]:
            return r["encloses_id"], True
        return cid, False

    def ancestors(self, cid: str) -> list[str]:
        """cid, its parent, ..., the root."""
        out = [cid]
        while self.rows[out[-1]]["parent_id"] and self.rows[out[-1]]["parent_id"] not in out:
            out.append(self.rows[out[-1]]["parent_id"])
        return out

    def root(self) -> str:
        return next(c for c, r in self.rows.items() if r["kind"] == "space" and not r["parent_id"])

    def locate(self, compartment_ids) -> tuple[str, str, str]:
        """Where a node touching these compartments is drawn: (space id, kind, outer space id).

        inside:   every participant is in one space             -> (that space, "inside", "")
        boundary: the spaces lie on one root-to-leaf path, or a participant is in a membrane
                  -> (the innermost space, "boundary", the space on the other side)
        between:  the spaces are siblings or cousins            -> (their lowest common ancestor, "between", "")
        """
        cids = sorted({c for c in compartment_ids if c in self.rows})
        if not cids:
            return self.root(), "inside", ""
        pos = [self.space_of(c) for c in cids]
        spaces = sorted({s for s, _ in pos})
        on_membrane = any(b for _, b in pos)
        if len(spaces) == 1:
            s = spaces[0]
            return (s, "boundary", self.rows[s]["parent_id"]) if on_membrane else (s, "inside", "")
        deepest = max(spaces, key=lambda s: (len(self.ancestors(s)), s))
        chain = self.ancestors(deepest)
        if all(s in chain for s in spaces):
            outer = next(a for a in chain[1:] if a in spaces)
            return deepest, "boundary", outer
        common = [a for a in chain if all(a in self.ancestors(s) for s in spaces)]
        return (common[0] if common else self.root()), "between", ""


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    from sbml_features import extract  # noqa: E402

    _nodes, _edges, comps = extract(sys.argv[1])
    for r in comps:
        print("  " * int(r["depth"]) + f"{r['compartment_id']} [{r['name']}] {r['kind']}"
              + (f" encloses {r['encloses_id']}" if r["encloses_id"] else "") + (" (virtual)" if r["virtual"] else "")
              + f"  <- {r['source']}")
