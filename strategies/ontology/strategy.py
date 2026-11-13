"""Run (from the bundle root; one model, then check and export):
    python strategies/ontology/strategy.py --features panel/features/skm_ja --sbml panel/models/skm_ja.xml --out strategies/ontology/out/skm_ja
    python common/contract.py validate panel/features/skm_ja strategies/ontology/out/skm_ja
    python common/contract.py write-sbml panel/models/skm_ja.xml strategies/ontology/out/skm_ja --features panel/features/skm_ja --strategy ontology -o strategies/ontology/out/skm_ja/skm_ja_grouped.xml

All six models (PowerShell):
    foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") { python strategies/ontology/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/ontology/out/$m }

Ontology grouping strategy: boxes from the controlled identifiers an SBML model carries - named bins in
the notes (MapMan), the EC tree, pathway cross-references, GO terms, and SBO as the fallback.

Reads only nodes.tsv, edges.tsv and compartments.tsv of --features (the --sbml file is accepted, as the
contract requires, but everything needed is already in the feature tables). Writes groups.tsv and
groups_meta.tsv into --out. The rules R1..R15 named in the comments are described in RULES.md.
No random numbers are used and every choice is made over sorted keys, so two runs are byte-identical.
"""
from __future__ import annotations

import argparse
import collections
import csv
import html
import json
import pathlib
import re

# ---------------------------------------------------------------------------------------------------
# Parameters (RULES.md, parameter table). The same values for every model.
M_MIN = 3        # a box needs at least 3 reactions (function boxes) or 3 nodes (localisation boxes)
PAGE = 12        # at most 12 sibling boxes: the top level, and the children of any box
S_MAX = 0.5      # a level-1 box holding more than half of all nodes is a catch-all and is opened
CLONE_MIN_BOXES = 3   # a species is cloned when its reactions lie in at least 3 level-1 boxes

# Ontologies in the order a node tries them (R4). SBO is not in this list: it is the fallback (R9).
ORDER = ("bins", "ec", "pathway", "go")
HERITABLE = {"bins", "pathway", "go"}      # process terms pass through shared species (R3, R8); EC does not
UNASSIGNED = "unassigned"

# The seven top-level classes of the Enzyme Commission tree. The SBML carries EC numbers but no EC
# names; these seven class names are the only vocabulary embedded in this file. Deeper EC terms are
# labelled by their number alone.
EC_CLASS = {"1": "oxidoreductases", "2": "transferases", "3": "hydrolases", "4": "lyases",
            "5": "isomerases", "6": "ligases", "7": "translocases"}

# identifiers.org namespaces whose identifiers are pathways whatever the qualifier
PATHWAY_DBS = {"kegg.pathway", "wikipathways", "smpdb", "panther.pathway", "unipathway", "pw", "pathwaycommons"}
KEGG_MAP = re.compile(r"^(?:path:)?[a-z]{2,4}(\d{5})$")            # map00908, ko04075, ath00906 -> 00908 ...
CYC_PATHWAY = re.compile(r"(?:^|-)PWY\d*(?:-|$)")                   # PWY-4203, ETHYL-PWY, PWY0-1234
BIN = re.compile(r"^(\d+(?:\.\d+)*)_(.+)$")                         # "10.7.1_Phytohormone action.jasmonic acid.biosynthesis"
EC_NUM = re.compile(r"^\d+(?:\.\d+){0,3}$")
EC_IN_BIN = re.compile(r"^EC_(\d+(?:-\d+){0,3})\s*(.*)$")           # bin name "EC_2-7 transferase transferring ..."


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def natural(text):
    """Sort key: numbers inside an id compare as numbers (mapman:10.7 after mapman:9.2)."""
    return [(0, int(p)) if p.isdigit() else (1, p) for p in re.split(r"(\d+)", text) if p]


class Terms:
    """Every ontology term met in the model: kind, parent, name pieces, the identifier behind it."""

    def __init__(self):
        self.kind, self.parent, self.names, self.xref, self.source = {}, {}, {}, {}, {}
        self.has_child = set()

    def add(self, tid, kind, parent, names, xref, source):
        if tid not in self.kind:
            self.kind[tid], self.parent[tid], self.xref[tid], self.source[tid] = kind, parent, xref, source
            self.names[tid] = names
        elif names and not self.names[tid]:
            self.names[tid] = names
        if parent:
            self.has_child.add(parent)

    def base(self, gid):
        return gid[:-len("~other")] if gid.endswith("~other") else gid

    def kind_of(self, gid):
        return self.kind.get(self.base(gid), "")


# ---------------------------------------------------------------------------------------------------
# R1  term paths of a node, from its own annotation
def own_paths(node, terms: Terms) -> dict[str, list[tuple]]:
    """{ontology kind: [path, ...]}; a path is a tuple of term ids from the top of the tree to the leaf."""
    out = collections.defaultdict(set)
    notes = json.loads(node["notes"]) if node["notes"] else {}

    # bins: any notes key whose value is a ";"-list of "<dotted code>_<dotted names>" (MapMan style)
    for key in sorted(notes):
        parts = [html.unescape(p.strip()) for p in notes[key].split(";") if p.strip()]
        if not parts or not all(BIN.match(p) for p in parts):
            continue
        ns = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")
        for p in parts:
            code, name = BIN.match(p).groups()
            c, nm = code.split("."), name.split(".")
            path = []
            ec_from = next((i for i, seg in enumerate(nm) if EC_IN_BIN.match(seg)), None) if len(nm) == len(c) else None
            if ec_from is not None and all(EC_IN_BIN.match(seg) for seg in nm[ec_from:]):
                # the EC tree written as a bin ("Enzyme classification.EC_2 transferases.EC_2-7 ..."):
                # read as EC terms, named by the text the notes give after the EC number
                for seg in nm[ec_from:]:
                    num, text = EC_IN_BIN.match(seg).groups()
                    tid = "ec-code:" + num.replace("-", ".")
                    terms.add(tid, "ec", path[-1] if path else "", (text.strip(),) if text.strip() else (), tid, "ec")
                    path.append(tid)
                out["ec"].add(tuple(path))
                continue
            for i in range(1, len(c) + 1):
                tid = f"{ns}:{'.'.join(c[:i])}"
                names = tuple(nm[:i]) if len(nm) == len(c) else ((name,) if i == len(c) else ())
                terms.add(tid, "bins", path[-1] if path else "", names, tid, "notes")
                path.append(tid)
            out["bins"].add(tuple(path))

    # EC: annotation (ec column) and the notes key ec_numbers
    ecs = set(filter(None, node["ec"].split("|")))
    ecs |= {e.strip() for e in re.split(r"[;,]", notes.get("ec_numbers", "")) if e.strip()}
    for ec in sorted(ecs):
        digits = [d for d in ec.split(".") if d not in ("-", "")]
        code = ".".join(digits)
        if not EC_NUM.match(code):
            continue
        path = []
        for i in range(1, len(digits) + 1):
            tid = "ec-code:" + ".".join(digits[:i])
            terms.add(tid, "ec", path[-1] if path else "", (), tid, "ec")
            path.append(tid)
        out["ec"].add(tuple(path))

    # pathway cross-references (flat): recognised by the syntax of the database and identifier
    for x in filter(None, node["xrefs"].split("|")):
        qual, _, res = x.partition("=")
        db, _, ident = res.partition(":")
        tid = ""
        if db in PATHWAY_DBS:
            tid = f"{db}:{ident}"
        elif db == "kegg" and KEGG_MAP.match(ident):
            tid = "kegg.pathway:map" + KEGG_MAP.match(ident).group(1)
        elif db in ("biocyc", "metacyc") and CYC_PATHWAY.search(ident.split(":")[-1]):
            tid = "biocyc:" + ident.split(":")[-1]
        elif db == "uri" and re.search(r"[?&]object=([^&#]+)", ident):
            frame = re.search(r"[?&]object=([^&#]+)", ident).group(1)
            if CYC_PATHWAY.search(frame):
                tid = "biocyc:" + frame
        elif db == "reactome" and qual.endswith("isPartOf"):
            tid = f"reactome:{ident}"
        if tid:
            terms.add(tid, "pathway", "", (), tid, "xref")
            out["pathway"].add((tid,))

    # GO (flat: the bundle has no GO graph, so no depth can be chosen)
    for g in filter(None, node["go"].split("|")):
        terms.add(g, "go", "", (), g, "go")
        out["go"].add((g,))

    # SBO path of a reaction, without the ontology root (used by R9 only)
    if node["kind"] == "reaction" and node["sbo_path"]:
        ids = node["sbo_path"].split(">")
        nms = node["sbo_path_names"].split(" > ")
        path = []
        for i, t in enumerate(ids):
            if i == 0 and len(ids) > 1:
                continue
            terms.add(t, "sbo", path[-1] if path else "", (nms[i] if i < len(nms) else t,), t, "sbo")
            path.append(t)
        out["sbo"].add(tuple(path))
    return {k: sorted(v) for k, v in out.items()}


class Model:
    def __init__(self, features_dir):
        f = pathlib.Path(features_dir)
        self.nodes = {n["node_id"]: n for n in read_tsv(f / "nodes.tsv")}
        self.comps = {c["compartment_id"]: c for c in read_tsv(f / "compartments.tsv")}
        self.species = sorted(n for n, r in self.nodes.items() if r["kind"] == "species")
        self.reactions = sorted(n for n, r in self.nodes.items() if r["kind"] == "reaction")
        self.terms = Terms()
        self.paths = {n: own_paths(self.nodes[n], self.terms) for n in sorted(self.nodes)}
        self.role = {r: {"reactant": set(), "product": set(), "modifier": set()} for r in self.reactions}
        self.rx_of = {s: set() for s in self.species}
        for e in read_tsv(f / "edges.tsv"):
            if e["reaction_id"] in self.role and e["species_id"] in self.rx_of:
                self.role[e["reaction_id"]][e["role"]].add(e["species_id"])
                self.rx_of[e["species_id"]].add(e["reaction_id"])
        self.parts = {r: sorted(set().union(*self.role[r].values())) for r in self.reactions}

        # R1-R3: where a reaction's terms of each ontology come from
        self.src = {}          # reaction -> kind -> (rule id, [source node ids])
        for r in self.reactions:
            self.src[r] = {}
            mods = sorted(self.role[r]["modifier"])
            subs = sorted(self.role[r]["reactant"] | self.role[r]["product"])
            for kind in ORDER:
                if self.paths[r].get(kind):
                    self.src[r][kind] = ("R1", [r])
                elif any(self.paths[s].get(kind) for s in mods):
                    self.src[r][kind] = ("R2", [s for s in mods if self.paths[s].get(kind)])
                elif kind in HERITABLE and any(self.paths[s].get(kind) for s in subs):
                    self.src[r][kind] = ("R3", [s for s in subs if self.paths[s].get(kind)])

        # main metabolites named in a notes key of the reaction (R12): token -> the one species id ending with it
        self.main = {r: set() for r in self.reactions}
        for r in self.reactions:
            notes = json.loads(self.nodes[r]["notes"]) if self.nodes[r]["notes"] else {}
            for tok in (t.strip() for t in re.split(r"[;,]", notes.get("main_metabolites", "")) if t.strip()):
                hit = [s for s in self.parts[r] if s == tok or s.endswith("_" + tok)]
                if len(hit) == 1:
                    self.main[r].add(hit[0])

    def loc_key(self, n):
        row = self.nodes[n]
        return row["loc_space"] + ("@boundary" if row["loc_kind"] == "boundary" else "")


# ---------------------------------------------------------------------------------------------------
# R5/R7  which box a path leads to, given the terms that are opened and the boxes that are dropped
def box_of(path, start, opened, banned):
    """First term of path[start:] that is not opened and not dropped; else the "rest" box of an opened term."""
    chain = []
    for t in path[start:]:
        if t in opened:
            chain.append(t)
            continue
        if t not in banned:
            return t
        break
    for h in reversed(chain):
        if h + "~other" not in banned:
            return h + "~other"
    return None


def candidates(m: Model, node, kind, start_of, opened, banned):
    """Boxes the node's own paths of this ontology lead to. start_of(path) -> index, or None to skip the path."""
    out = set()
    for p in m.paths[node].get(kind, ()):
        i = start_of(p)
        if i is not None:
            b = box_of(p, i, opened, banned)
            if b:
                out.add(b)
    return out


def choose_boxes(m: Model, reactions, start_of, opened, banned):
    """R4: every reaction picks one box. Returns {reaction: (box, rule id)} for those that have a candidate."""
    cache = {}

    def cand(node, kind):
        if (node, kind) not in cache:
            cache[(node, kind)] = candidates(m, node, kind, start_of, opened, banned)
        return cache[(node, kind)]

    rule, top = {}, {}
    support = collections.Counter()
    for r in reactions:
        for kind in ORDER:                       # R4a: first ontology of the order that offers a box
            if kind not in m.src[r]:
                continue
            rid, sources = m.src[r][kind]
            votes = collections.Counter(b for s in sources for b in sorted(cand(s, kind)))
            if votes:                            # R4b: the box named by most of the reaction's sources
                best = max(votes.values())
                rule[r], top[r] = rid, sorted((b for b in votes if votes[b] == best), key=natural)
                support.update(votes.keys())
                break
    placed = {r: top[r][0] for r in top if len(top[r]) == 1}
    pending = [r for r in reactions if r in top and r not in placed]
    while pending:
        near = sorted({s for r in pending for s in m.parts[r]})
        homes = species_homes(m, placed, start_of, opened, banned, allowed=near)
        new = {}
        for r in pending:                        # R4c: the candidate in which most participants already live
            v = collections.Counter(homes[s][0] for s in m.parts[r] if s in homes and homes[s][0] in top[r])
            if v:
                best = max(v.values())
                win = [b for b in v if v[b] == best]
                if len(win) == 1:
                    new[r] = win[0]
        if not new:                              # R4d: nothing decides - the least ambiguous reactions take
            fewest = min(len(top[r]) for r in pending)          # the candidate that already holds more
            size = collections.Counter(placed.values())         # reactions, then R4c resumes
            for r in pending:
                if len(top[r]) == fewest:
                    new[r] = min(top[r], key=lambda b: (-size[b], -support[b], natural(b)))
        placed.update(new)
        pending = [r for r in pending if r not in placed]
    return {r: (placed[r], rule[r]) for r in reactions if r in placed}, support


def species_homes(m: Model, r_box, start_of, opened, banned, allowed=None):
    """R11/R12: {species: (box, rule)} from the boxes of its reactions (r_box may map to None = no box)."""
    exists = {b for b in r_box.values() if b}
    count = collections.Counter(b for b in r_box.values() if b)
    out = {}
    for s in m.species if allowed is None else allowed:
        rx = sorted(r for r in m.rx_of[s] if r in r_box)
        main = [r for r in rx if s in m.main[r]]
        weight = collections.defaultdict(float)
        for r in (main or rx):
            weight[r_box[r]] += 1.0 / max(1, len(m.parts[r]))
        own = set()
        for kind in ORDER:                       # R11: the species' own term, if that box exists
            own = candidates(m, s, kind, start_of, opened, banned) & exists
            if own:
                break
        if own:
            best = min(own, key=lambda b: (-round(weight.get(b, 0.0), 9), -count[b], natural(b)))
            out[s] = (best, "R11")
        elif weight:                             # R12: the box holding most of its reactions (weighted)
            best = min(weight, key=lambda b: (-round(weight[b], 9), b is None, natural(b or "")))
            out[s] = (best, "R12")
    return out


def drop_first(m: Model, small, n_rx, support):
    """Which too-small box is dropped next.

    A box can only gain reactions from boxes of the same or an earlier ontology and from terms below it,
    so boxes are judged in ontology order and deepest term first. Within that, a box that could not reach
    M_MIN even if every reaction naming it joined goes before one that could; then the smallest.
    """
    def key(b):
        t = m.terms
        kind = t.kind_of(b)
        rank = ORDER.index(kind) if kind in ORDER else len(ORDER)
        can_grow = max(support.get(b, 0), n_rx[b]) >= M_MIN
        return (rank, -len(_ancestors(t, t.base(b))), can_grow, n_rx[b], natural(b))
    return min(small, key=key)


def majority(m: Model, r, homes, ok):
    """Box most participants of r live in (only boxes passing ok); ties: more modifiers, then the id."""
    votes, mods = collections.Counter(), collections.Counter()
    for s in m.parts[r]:
        b = homes.get(s, (None,))[0]
        if b and ok(b):
            votes[b] += 1
            mods[b] += s in m.role[r]["modifier"]
    if not votes:
        return None
    return min(votes, key=lambda b: (-votes[b], -mods[b], natural(b)))


# ---------------------------------------------------------------------------------------------------
def place_level1(m: Model, opened, banned):
    top = lambda p: 0
    r_box, r_rule = {}, {}
    chosen, support = choose_boxes(m, m.reactions, top, opened, banned)
    for r, (b, rule) in chosen.items():
        r_box[r], r_rule[r] = b, rule

    def spread(ok, rule):
        while True:
            homes = species_homes(m, r_box, top, opened, banned)
            new = {}
            for r in m.reactions:
                if r not in r_box:
                    b = majority(m, r, homes, ok)
                    if b:
                        new[r] = b
            if not new:
                return
            for r, b in new.items():
                r_box[r], r_rule[r] = b, rule

    # R8: process boxes spread to reactions without a term, through shared species
    spread(lambda b: m.terms.kind_of(b) in HERITABLE, "R8")
    # R9: reaction type - the deepest SBO term of the reaction that is still a box
    for r in m.reactions:
        if r not in r_box:
            for p in m.paths[r].get("sbo", ()):
                for t in reversed(p):
                    if t not in banned:
                        r_box[r], r_rule[r] = t, "R9"
                        break
    # R10: last resort - the box of most participants, whatever its ontology; else "unassigned"
    spread(lambda b: True, "R10")
    for r in m.reactions:
        if r not in r_box:
            r_box[r], r_rule[r] = UNASSIGNED, "R10"
    homes = species_homes(m, r_box, top, opened, banned)
    s_box = {s: homes.get(s, (UNASSIGNED, "R10"))[0] for s in m.species}
    s_rule = {s: homes.get(s, (UNASSIGNED, "R10"))[1] for s in m.species}
    return r_box, r_rule, s_box, s_rule, support


def only_child(m: Model, box, members, start_of, opened, banned):
    """R5 compression: the single child term that every seeded member reaches through `box`, else None."""
    nxt = set()
    for r in members:
        for kind in ORDER:
            if kind not in m.src[r]:
                continue
            for s in m.src[r][kind][1]:
                for p in m.paths[s].get(kind, ()):
                    i = start_of(p)
                    if i is None or box_of(p, i, opened, banned) != box:
                        continue
                    j = p.index(box)
                    nxt.add(p[j + 1] if j + 1 < len(p) else None)
    return next(iter(nxt)) if len(nxt) == 1 and None not in nxt else None


def level1(m: Model):
    opened, banned = set(), set()
    top = lambda p: 0
    n_all = len(m.nodes)
    m.log = []                                   # (what happened, box, reactions in it) - for NOTES.md
    while True:
        r_box, r_rule, s_box, s_rule, support = place_level1(m, opened, banned)
        n_rx = collections.Counter(r_box.values())
        boxes = sorted(n_rx, key=natural)
        droppable = [b for b in boxes if b != UNASSIGNED]
        small = [b for b in droppable if n_rx[b] < M_MIN]                       # R5: fewer than M_MIN reactions
        if small:
            victim = drop_first(m, small, n_rx, support)
            m.log.append(("R5 too small", victim, n_rx[victim]))
            banned.add(victim)
            continue
        if len(boxes) > PAGE and droppable:                                     # R6: the page is full
            victim = min(droppable, key=lambda b: (n_rx[b], natural(b)))
            m.log.append(("R6 page full", victim, n_rx[victim]))
            banned.add(victim)
            continue
        tree = [b for b in boxes if m.terms.kind_of(b) in ("bins", "ec") and not b.endswith("~other")]
        chain = [b for b in tree if only_child(m, b, [r for r in m.reactions if r_box[r] == b and r_rule[r] in
                                                      ("R1", "R2", "R3")], top, opened, banned)]
        if chain:                                                               # R5: single-child chain
            m.log.append(("R5 single child", chain[0], n_rx[chain[0]]))
            opened.add(chain[0])
            continue
        n_nodes = collections.Counter(list(r_box.values()) + list(s_box.values()))
        big = [b for b in tree if n_nodes[b] > S_MAX * n_all and b in m.terms.has_child]
        if big:                                                                 # R7: catch-all
            victim = max(big, key=lambda b: (n_nodes[b], natural(b)))
            m.log.append(("R7 catch-all opened", victim, n_rx[victim]))
            opened.add(victim)
            continue
        return r_box, r_rule, s_box, s_rule, opened


def level2(m: Model, box, rx, sp, r_rule, opened1):
    """R13: inside one level-1 box of a tree ontology, the next level of that tree. {node: child box or None}."""
    if m.terms.kind_of(box) not in ("bins", "ec") or box.endswith("~other"):
        return {}
    seeds = [r for r in rx if r_rule[r] in ("R1", "R2", "R3")]
    opened, banned = set(), set()

    def below(p):
        return p.index(box) + 1 if box in p else None

    while True:
        chosen, support = choose_boxes(m, seeds, below, opened, banned)
        r_sub = {r: b for r, (b, _rule) in chosen.items()}
        for r in rx:
            r_sub.setdefault(r, None)
        homes = species_homes(m, r_sub, below, opened, banned, allowed=sp)
        for r in rx:                                         # reactions placed by R8/R10 follow their species
            if r not in seeds:
                r_sub[r] = majority(m, r, homes, lambda b: True)
        homes = species_homes(m, r_sub, below, opened, banned, allowed=sp)
        n_rx = collections.Counter(b for b in r_sub.values() if b)
        small = [b for b in n_rx if n_rx[b] < M_MIN]
        if small:
            banned.add(drop_first(m, small, n_rx, support))
            continue
        if len(n_rx) > PAGE:
            banned.add(min(n_rx, key=lambda b: (n_rx[b], natural(b))))
            continue
        chain = [b for b in sorted(n_rx, key=natural) if not b.endswith("~other")
                 and only_child(m, b, [r for r in seeds if r_sub[r] == b], below, opened, banned)]
        if chain:
            opened.add(chain[0])
            continue
        if len(n_rx) == 1 and max(n_rx.values()) <= S_MAX * len(rx):
            return {}                                        # a lone sub-box holding a minority divides nothing
        out = dict(r_sub)
        out.update({s: homes[s][0] if s in homes else None for s in sp})
        return out


# ---------------------------------------------------------------------------------------------------
def names_of(t: Terms, tid):
    """Name of an EC term: the text the SBML notes give for it, else the class name for the seven classes."""
    if t.names.get(tid):
        return t.names[tid][0]
    code = tid.split(":", 1)[1]
    return EC_CLASS.get(code, "") if "." not in code else ""


def label_of(m: Model, gid, parent_term=""):
    t = m.terms
    base = t.base(gid)
    kind = t.kind[base]
    if kind == "bins":
        names = t.names[base]
        skip = len(t.names.get(parent_term, ())) if parent_term and names[:len(t.names.get(parent_term, ()))] == \
            t.names.get(parent_term, ()) else 0
        text = ".".join(names[skip:]) if names else base
    elif kind == "ec":
        code = base.split(":", 1)[1]
        own = names_of(t, base)
        cls = EC_CLASS.get(code.split(".")[0], "")
        text = f"EC {code}" + (f" ({own})" if own else f" ({cls})" if cls and not parent_term else "")
    elif kind == "pathway":
        db, ident = base.split(":", 1)
        text = ("KEGG pathway " if db == "kegg.pathway" else db + " pathway ") + ident
    elif kind == "sbo":
        text = t.names[base][0]
    else:
        text = base
    return text


def run(features_dir, out_dir):
    m = Model(features_dir)
    r_box, r_rule, s_box, s_rule, opened1 = level1(m)
    l1 = dict(r_box)
    l1.update(s_box)
    rule = dict(r_rule)
    rule.update(s_rule)
    members = collections.defaultdict(list)
    for n in sorted(l1):
        members[l1[n]].append(n)

    l2 = {}
    for box in sorted(members, key=natural):
        rx = [n for n in members[box] if m.nodes[n]["kind"] == "reaction"]
        sp = [n for n in members[box] if m.nodes[n]["kind"] == "species"]
        l2.update(level2(m, box, rx, sp, r_rule, opened1))

    # function path of every node, then R14: localisation as the innermost level
    func = {n: (l1[n], l2[n]) if l2.get(n) else (l1[n],) for n in l1}
    has_kids = {f[0] for f in func.values() if len(f) == 2}
    by_func = collections.defaultdict(list)
    for n in sorted(func):
        by_func[func[n]].append(n)
    loc = {}
    for f in sorted(by_func, key=lambda x: [natural(p) for p in x]):
        if len(f) == 1 and f[0] in has_kids:
            continue                                         # nodes left directly in a box that has sub-boxes
        part = collections.defaultdict(list)
        for n in by_func[f]:
            part[m.loc_key(n)].append(n)
        if len(part) < 2:
            continue                                         # the box lies in one place: nothing to split
        for key in sorted(part):                             # a small boundary part joins its space
            space = key.split("@")[0]
            if key.endswith("@boundary") and len(part[key]) < M_MIN and space in part:
                part[space] += part.pop(key)
        keep = [key for key in sorted(part) if len(part[key]) >= M_MIN]
        if len(keep) < 2:
            continue                                         # one place only (plus strays): not divided
        for key in keep:
            for n in part[key]:
                loc[n] = key

    def gid_of(path):
        return "/".join(path)

    rows, meta = [], {}
    chain = {}
    for n in sorted(func):
        path = list(func[n]) + ([loc[n]] if n in loc else [])
        ids = [gid_of(path[:i + 1]) for i in range(len(path))]
        chain[n] = ids
        for i, gid in enumerate(ids):
            if gid in meta:
                continue
            is_loc = n in loc and i == len(path) - 1
            if is_loc:
                cid = loc[n].split("@")[0]
                name = m.comps[cid]["name"] if cid in m.comps else cid
                meta[gid] = dict(group_id=gid, level=i + 1, parent=ids[i - 1], basis="localisation",
                                 label=name + (" boundary" if loc[n].endswith("@boundary") else ""),
                                 name_source="compartment", name_xref=cid, compartment_id=loc[n])
            elif path[i] == UNASSIGNED:
                # R10 leftovers: named by the deepest SBO term all their reactions share, if they all have one
                sbo = [m.paths[r]["sbo"][0] for r in m.reactions if l1[r] == UNASSIGNED and m.paths[r].get("sbo")]
                common = []
                if sbo and len(sbo) == sum(1 for r in m.reactions if l1[r] == UNASSIGNED):
                    for level in zip(*sbo):
                        if len(set(level)) > 1:
                            break
                        common.append(level[0])
                if common:
                    meta[gid] = dict(group_id=gid, level=1, parent="", basis="reaction_type",
                                     label=m.terms.names[common[-1]][0], name_source="sbo", name_xref=common[-1],
                                     compartment_id="")
                else:
                    meta[gid] = dict(group_id=gid, level=1, parent="", basis="mixed", label="", name_source="none",
                                     name_xref="", compartment_id="")
            else:
                base = m.terms.base(path[i])
                kind = m.terms.kind[base]
                meta[gid] = dict(group_id=gid, level=i + 1, parent=ids[i - 1] if i else "",
                                 basis="reaction_type" if kind == "sbo" else "function",
                                 label=label_of(m, path[i], m.terms.base(path[0]) if i else ""),
                                 name_source=m.terms.source[base], name_xref=m.terms.xref[base], compartment_id="",
                                 _rest=path[i].endswith("~other"), _base=base, _parent_term=m.terms.parent[base])
    # a "rest" box says "(other)" only when a sibling box of the same opened term exists beside it
    for gid, row in meta.items():
        if row.get("_rest"):
            sib = [g for g, r in meta.items() if g != gid and r["parent"] == row["parent"] and r.get("_base")
                   and row["_base"] in _ancestors(m.terms, r["_base"])]
            if sib:
                row["label"] += " (other)"

    # R15: clones of species whose reactions lie in at least CLONE_MIN_BOXES level-1 boxes
    func_gid = {n: gid_of(func[n]) for n in func}
    loc_groups = {gid for gid, r in meta.items() if r["basis"] == "localisation"}
    for n in sorted(func):
        clone = ""
        if m.nodes[n]["kind"] == "species" and len({l1[r] for r in m.rx_of[n]}) >= CLONE_MIN_BOXES:
            targets = []
            for b in sorted({func_gid[r] for r in m.rx_of[n] if l1[r] != l1[n]}, key=natural):
                inner = f"{b}/{m.loc_key(n)}"                # a copy sits in the function box of the reaction,
                plain = f"{b}/{m.nodes[n]['loc_space']}"     # in the part for the species' own place if there is one
                targets.append(inner if inner in loc_groups else plain if plain in loc_groups else b)
            clone = "|".join(sorted(set(targets)))
        ids = chain[n] + [""] * (3 - len(chain[n]))
        rows.append(dict(node_id=n, L1=ids[0], L2=ids[1], L3=ids[2], rule=rule[n], clone_in=clone))

    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", sorted(rows, key=lambda r: r["node_id"]),
              ["node_id", "L1", "L2", "L3", "rule", "clone_in"])
    write_tsv(out / "groups_meta.tsv", [meta[g] for g in sorted(meta)],
              ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"])
    return rows, meta, m.log


def _ancestors(terms: Terms, tid):
    out = []
    while terms.parent.get(tid):
        tid = terms.parent[tid]
        out.append(tid)
    return out


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Ontology grouping strategy (see RULES.md)")
    ap.add_argument("--features", required=True, help="panel/features/<model>")
    ap.add_argument("--sbml", required=True, help="panel/models/<model>.xml (accepted; the feature tables suffice)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    a = ap.parse_args(argv)
    rows, meta, log = run(a.features, a.out)
    levels = collections.Counter(r["level"] for r in meta.values())
    print(f"{pathlib.Path(a.features).name}: {len(rows)} nodes, groups per level "
          f"{levels.get(1, 0)}/{levels.get(2, 0)}/{levels.get(3, 0)}, "
          f"{sum(1 for x in log if x[0] == 'R6 page full')} boxes dropped by the page rule -> {a.out}")


if __name__ == "__main__":
    main()
