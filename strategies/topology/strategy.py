"""Topology strategy: connectivity-driven modules of the species-reaction graph.

Run from the bundle root (the folder that holds panel/ and common/), with this folder unpacked as
strategies/topology/. One model, one algorithm:
    python strategies/topology/strategy.py --features panel/features/skm_ja --sbml panel/models/skm_ja.xml --out strategies/topology/topology_leiden/out/skm_ja --algorithm leiden
    python common/contract.py validate panel/features/skm_ja strategies/topology/topology_leiden/out/skm_ja
    python common/contract.py write-sbml panel/models/skm_ja.xml strategies/topology/topology_leiden/out/skm_ja --features panel/features/skm_ja --strategy topology_leiden -o strategies/topology/topology_leiden/out/skm_ja/skm_ja_grouped.xml

All six models and the three algorithms, PowerShell:
    foreach ($a in "leiden","infomap","walktrap") { foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") {
        $o = "strategies/topology/topology_$a/out/$m"
        python strategies/topology/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out $o --algorithm $a
        python common/contract.py validate panel/features/$m $o
        python common/contract.py write-sbml panel/models/$m.xml $o --features panel/features/$m --strategy topology_$a -o "$o/${m}_grouped.xml" } }

The same in bash:
    for a in leiden infomap walktrap; do for m in brenda_merged e_coli_core R-HSA-77289 skm_pss skm_ja skm_sa; do
        o=strategies/topology/topology_$a/out/$m
        python strategies/topology/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out $o --algorithm $a
        python common/contract.py validate panel/features/$m $o
        python common/contract.py write-sbml panel/models/$m.xml $o --features panel/features/$m --strategy topology_$a -o $o/${m}_grouped.xml
    done; done

Options beyond the contract's three:
    --algorithm leiden|infomap|walktrap   community detection used in rule R4 (default leiden, the recommended variant)
    --seed N                               seed of Leiden and Infomap (default 1; Walktrap has no random step)
    --stats FILE                           also write a JSON with the numbers quoted in NOTES.md

Rules (ids as in RULES.md; the three algorithms share everything except R4):
    R1  graph       species and reactions are the nodes; one link per (reaction, species) participant pair
    R2  candidates  species with a degree that no species would reach if links were spread at random (Poisson)
    R3  currency    a candidate is cut out of the graph unless most of its reactions would be stranded without it
    R4  modules     Leiden | Infomap | Walktrap on the graph without currency species
    R5  merge       a part with fewer than S_MIN nodes joins the part it is most strongly connected to
    R6  page        while there are more than K_MAX parts, the two connected parts with the largest modularity gain merge
    R7  next level  a leaf group with more than S_MAX nodes is partitioned by R4-R6 on its own subgraph (down to L3)
    R8  currency species: home = the leaf group holding most of its reactions, copies (clone_in) in the others
    R9  names       an annotation term (classification code in the notes, SBO, EC, GO) carried by most of the group
                    and found mostly in the group; else the hub member; else unnamed

Reads nodes.tsv and edges.tsv only. --sbml is accepted because the contract fixes the command line; the
file is not opened, everything needed is in the feature tables. Compartments are not used.
Needs python-igraph and leidenalg; nothing else outside the standard library.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import pathlib
import random
import re

import igraph as ig
import leidenalg

# ---- parameters (RULES.md gives the reason for each value; none is set per model) -----------------
CHANCE = 1.0         # R2: a degree is beyond chance when fewer than this many species would reach it at random
MAJORITY = 0.5       # R3, R9: "most of". R3: a reaction is stranded when it keeps less than this share of its
                     # connected network, and a candidate stays in the graph when at least this share of its
                     # reactions is stranded. R9: share of the group that carries a term, share of the term's
                     # carriers that lie in the group
S_MIN = 5            # R5: smallest module = two reactions in a row with their three species
K_MAX = 12           # R6: groups per page = classes of the largest qualitative colour palettes, a 4 x 3 grid
S_MAX = K_MAX * S_MIN  # R7: a leaf that could hold a full page of smallest modules gets its own level
MAX_LEVEL = 3        # contract: L1..L3
WALKTRAP_STEPS = 4   # R4: length of the random walks (value of the original paper, igraph default)
INFOMAP_TRIALS = 10  # R4: Infomap keeps the best of this many attempts (igraph default)
ALGORITHMS = ("leiden", "infomap", "walktrap")
BASIS = {"leiden": "topology", "walktrap": "topology", "infomap": "flow"}

GROUP_COLUMNS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
META_COLUMNS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]
CODE = re.compile(r"^(\d+(?:\.\d+)*)_(.+)$")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def chance_bound(degrees):
    """R2. The smallest degree d that fewer than CHANCE species would reach if the model's links were
    spread over its species at random: n * P(Poisson(mean degree) >= d) < CHANCE."""
    n = len(degrees)
    if n == 0:
        return 1
    lam = sum(degrees) / n
    term = math.exp(-lam)  # P(X = 0)
    below = 0.0            # P(X < d)
    d = 0
    while n * (1.0 - below) >= CHANCE and d < max(degrees) + 1:
        below += term
        d += 1
        term *= lam / d
    return d


# ---- R1: the graph -----------------------------------------------------------------------------------
class Net:
    """Species-reaction graph. Node index = position in the sorted id list."""

    def __init__(self, nodes, edges):
        self.rows = {n["node_id"]: n for n in nodes}
        self.ids = sorted(self.rows)
        self.idx = {v: i for i, v in enumerate(self.ids)}
        self.is_species = [self.rows[v]["kind"] == "species" for v in self.ids]
        links, arcs = set(), set()
        for e in edges:
            if e["reaction_id"] not in self.idx or e["species_id"] not in self.idx:
                continue
            r, s = self.idx[e["reaction_id"]], self.idx[e["species_id"]]
            links.add((r, s))
            if e["role"] == "product":
                arcs.add((r, s))
            else:  # reactant or modifier: the species acts on the reaction
                arcs.add((s, r))
            if e["role"] != "modifier" and self.rows[e["reaction_id"]]["reversible"] == "1":
                arcs.add((r, s))
                arcs.add((s, r))
        self.links = sorted(links)
        self.arcs = sorted(arcs)
        self.adj = [[] for _ in self.ids]
        for a, b in self.links:
            self.adj[a].append(b)
            self.adj[b].append(a)
        for lst in self.adj:
            lst.sort()
        self.currency = set()
        self.fence = None
        self.candidates = []

    def degree(self, i):
        return len(self.adj[i])

    # ---- R2 + R3: currency species ------------------------------------------------------------------
    def find_currency(self):
        species = [i for i in range(len(self.ids)) if self.is_species[i]]
        self.fence = chance_bound([self.degree(i) for i in species])
        cand = sorted((i for i in species if self.degree(i) >= self.fence), key=lambda i: (-self.degree(i), self.ids[i]))
        removed = set()
        whole = self._component_sizes(set())
        for s in cand:
            comp_size = self._component_sizes(removed | {s})
            stranded = sum(1 for r in self.adj[s] if comp_size[r] < MAJORITY * whole[r])
            share = stranded / len(self.adj[s])
            cut = share < MAJORITY
            if cut:
                removed.add(s)
            self.candidates.append(dict(node_id=self.ids[s], name=self.rows[self.ids[s]]["name"], degree=self.degree(s),
                                        stranded_share=round(share, 3), currency=cut))
        self.currency = removed

    def _component_sizes(self, without):
        size = [0] * len(self.ids)
        seen = [False] * len(self.ids)
        for start in range(len(self.ids)):
            if seen[start] or start in without:
                continue
            stack, comp = [start], []
            seen[start] = True
            while stack:
                v = stack.pop()
                comp.append(v)
                for w in self.adj[v]:
                    if not seen[w] and w not in without:
                        seen[w] = True
                        stack.append(w)
            for v in comp:
                size[v] = len(comp)
        return size


# ---- R4: the three algorithms ------------------------------------------------------------------------
def detect(net, members, algorithm, seed):
    """Communities of the subgraph induced by `members` (sorted node indices). Returns a list of lists."""
    pos = {v: i for i, v in enumerate(members)}
    und = [(pos[a], pos[b]) for a, b in net.links if a in pos and b in pos]
    if not und:
        return [[v] for v in members]
    if algorithm == "infomap":
        arcs = [(pos[a], pos[b]) for a, b in net.arcs if a in pos and b in pos]
        g = ig.Graph(n=len(members), edges=arcs, directed=True)
        ig.set_random_number_generator(random.Random(seed))
        membership = g.community_infomap(trials=INFOMAP_TRIALS).membership
    else:
        g = ig.Graph(n=len(members), edges=und, directed=False)
        if algorithm == "leiden":
            membership = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition, seed=seed,
                                                  n_iterations=-1).membership
        elif algorithm == "walktrap":
            membership = g.community_walktrap(steps=WALKTRAP_STEPS).as_clustering().membership
        else:
            raise ValueError(algorithm)
    parts = collections.defaultdict(list)
    for v, c in zip(members, membership):
        parts[c].append(v)
    return sorted((sorted(p) for p in parts.values()), key=lambda p: p[0])


def modularity(net, members, parts):
    """Newman modularity of `parts` on the undirected subgraph induced by `members`."""
    where = {v: k for k, p in enumerate(parts) for v in p}
    inside = set(members)
    m = 0
    e_in = collections.Counter()
    deg = collections.Counter()
    for a, b in net.links:
        if a in inside and b in inside:
            m += 1
            deg[where[a]] += 1
            deg[where[b]] += 1
            if where[a] == where[b]:
                e_in[where[a]] += 1
    if m == 0:
        return 0.0
    return sum(e_in[k] / m - (deg[k] / (2 * m)) ** 2 for k in range(len(parts)))


# ---- R5 + R6: size repair shared by the three algorithms ----------------------------------------------
def link_counts(net, parts):
    """Links between parts in the graph without currency species: {(i, j): n}, i < j."""
    where = {v: k for k, p in enumerate(parts) for v in p}
    w = collections.Counter()
    for a, b in net.links:
        if a in where and b in where and where[a] != where[b]:
            w[tuple(sorted((where[a], where[b])))] += 1
    return w


def currency_links(net, parts):
    """Links between parts that run through currency species: {(i, j): weight}, i <= j.

    One-mode projection of each currency species c onto its d reactions inside these parts: every pair
    of its reactions gets weight 1 / (d - 1), so each reaction keeps the weight 1 its link to c had
    (Newman 2001). (i, i) holds the weight inside part i."""
    where = {v: k for k, p in enumerate(parts) for v in p}
    w = collections.defaultdict(float)
    for c in sorted(net.currency):
        n = collections.Counter(where[r] for r in net.adj[c] if r in where)
        d = sum(n.values())
        if d < 2:
            continue
        keys = sorted(n)
        for a, i in enumerate(keys):
            w[(i, i)] += n[i] * (n[i] - 1) / 2 / (d - 1)
            for j in keys[a + 1:]:
                w[(i, j)] += n[i] * n[j] / (d - 1)
    return w


def merge_small(net, parts):
    """R5. Returns (parts, nodes merged by direct links, nodes merged through currency species, island nodes)."""
    parts = [list(p) for p in parts]
    merged, stuck, via_currency = set(), set(), set()
    while True:
        small = sorted((k for k, p in enumerate(parts) if len(p) < S_MIN and k not in stuck),
                       key=lambda k: (len(parts[k]), parts[k][0]))
        if not small or len(parts) == 1:
            break
        k = small[0]
        w = link_counts(net, parts)
        cand = {(j if i == k else i): n for (i, j), n in w.items() if k in (i, j)}
        if not cand:
            cand = {(j if i == k else i): n for (i, j), n in currency_links(net, parts).items()
                    if k in (i, j) and i != j and n > 0}
            via_currency.update(parts[k] if cand else ())
        if not cand:
            stuck.add(k)
            continue
        best = sorted(cand, key=lambda j: (-cand[j], len(parts[j]), parts[j][0]))[0]
        merged.update(parts[k])
        parts[best] = sorted(parts[best] + parts[k])
        moved = {old: new for new, old in enumerate(i for i in range(len(parts)) if i != k)}
        stuck = {moved[s] for s in stuck}
        del parts[k]
    islands = {v for k in stuck for v in parts[k]}
    return sorted(parts, key=lambda p: p[0]), merged - via_currency, via_currency, islands


def merge_to_page(net, parts):
    """R6. Greedy merge of connected parts by modularity gain until at most K_MAX parts remain.

    Connection weight of two parts = direct links + links through currency species (currency_links)."""
    parts = [list(p) for p in parts]
    joined = set()
    while len(parts) > K_MAX:
        w = collections.defaultdict(float)
        for pair, n in link_counts(net, parts).items():
            w[pair] += n
        inner = collections.Counter()
        where = {v: k for k, p in enumerate(parts) for v in p}
        for a, b in net.links:
            if a in where and b in where and where[a] == where[b]:
                inner[where[a]] += 1
        for (i, j), n in currency_links(net, parts).items():
            if i == j:
                inner[i] += n
            else:
                w[(i, j)] += n
        w = {pair: n for pair, n in w.items() if n > 0}
        if not w:
            break
        m = sum(w.values()) + sum(inner.values())
        deg = [2 * inner[k] for k in range(len(parts))]
        for (i, j), n in w.items():
            deg[i] += n
            deg[j] += n

        def gain(pair):
            i, j = pair
            return w[pair] / m - deg[i] * deg[j] / (2 * m * m)

        i, j = sorted(w, key=lambda pr: (-gain(pr), len(parts[pr[0]]) + len(parts[pr[1]]), parts[pr[0]][0],
                                         parts[pr[1]][0]))[0]
        joined.update(parts[i] if len(parts[i]) <= len(parts[j]) else parts[j])
        parts[i] = sorted(parts[i] + parts[j])
        del parts[j]
    return sorted(parts, key=lambda p: p[0]), joined


def partition(net, members, algorithm, seed):
    """R4 -> R5 -> R6 on one node set. Returns (parts, {node: rule id}, raw parts of the algorithm)."""
    raw = detect(net, members, algorithm, seed)
    parts, merged, via_currency, islands = merge_small(net, raw)
    parts, joined = merge_to_page(net, parts)
    rule = {}
    for v in members:
        rule[v] = ("R5a" if v in merged else "R5b" if v in via_currency else "R5c" if v in islands
                   else "R6" if v in joined else "R4")
    return parts, rule, raw


# ---- R9: names ---------------------------------------------------------------------------------------
def node_terms(row):
    """Identifier-backed terms of one node: {(source, xref): (label, depth)}. Read from nodes.tsv columns
    notes (classification codes, any node), sbo_path / sbo_path_names (reactions only: a process type says
    what a group does, a material type does not), ec and go (any node). depth = level in the term's hierarchy."""
    out = {}
    if row["notes"]:
        try:
            notes = json.loads(row["notes"])
        except ValueError:
            notes = {}
        for key in sorted(notes):
            if key == "text":
                continue
            for piece in str(notes[key]).split(";"):
                mt = CODE.match(piece.strip())
                if not mt:
                    continue
                code, names = mt.group(1).split("."), mt.group(2).split(".")
                for d in range(1, len(code) + 1):
                    label = ".".join(names[:d]) if len(names) == len(code) else ""
                    out[("notes", f"{key}:{'.'.join(code[:d])}")] = (label.strip(), d)
    if row["sbo_path"] and row["kind"] == "reaction":
        ids, names = row["sbo_path"].split(">"), row["sbo_path_names"].split(" > ")
        if len(ids) == len(names):
            for d, (t, nm) in enumerate(zip(ids, names)):
                out[("sbo", t)] = (nm.strip(), d)
    for ec in filter(None, row["ec"].split("|")):
        digits = [p for p in ec.split(".") if p and p != "-"]
        for d in range(1, len(digits) + 1):
            out[("ec", "ec-code:" + ".".join(digits[:d]))] = ("EC " + ".".join(digits[:d]), d)
    for go in filter(None, row["go"].split("|")):
        out[("go", go)] = (go, 0)
    return out


class Namer:
    """R9a: a term names a group when at least half of the group's reactions (or of its species) carry it
    and the group holds more than half of all nodes of the model that carry it. R9b: else the hub: the
    species that takes part in most reactions of the group; when no species links two of its reactions, the
    reaction with most participants in the group; then any remaining species. R9c: else unnamed."""

    def __init__(self, net):
        self.net = net
        self.terms = [node_terms(net.rows[v]) for v in net.ids]
        self.total = collections.Counter(t for terms in self.terms for t in terms)
        labels = collections.defaultdict(collections.Counter)
        self.depth = {}
        for terms in self.terms:
            for t, (label, depth) in terms.items():
                self.depth[t] = depth
                if label:
                    labels[t][label] += 1
        self.label = {t: sorted(c, key=lambda x: (-c[x], x))[0] for t, c in labels.items()}
        self.used_labels = set()
        self.used_terms = set()

    def name(self, members):
        net = self.net
        n_sp = sum(1 for v in members if net.is_species[v])
        n_rx = len(members) - n_sp
        got_sp, got_rx = collections.Counter(), collections.Counter()
        for v in members:
            for t in self.terms[v]:
                (got_sp if net.is_species[v] else got_rx)[t] += 1
        cand = []
        for t in sorted(set(got_sp) | set(got_rx)):
            n = got_sp[t] + got_rx[t]
            cover = max(got_sp[t] / n_sp if n_sp else 0.0, got_rx[t] / n_rx if n_rx else 0.0)
            label = self.label.get(t, "")
            if not label and t[0] == "notes":
                label = t[1].replace(":", " ")
            if n >= 2 and cover >= MAJORITY and n / self.total[t] > MAJORITY and label \
                    and t not in self.used_terms and label not in self.used_labels:
                # the most specific term first: fewest carriers in the model, then the deepest in its hierarchy
                cand.append((self.total[t], -self.depth[t], t, label))
        if cand:
            _n, _d, t, label = sorted(cand)[0]
            self.used_terms.add(t)
            self.used_labels.add(label)
            return label, t[0], t[1]
        inside = set(members)
        hubs = []
        for v in members:
            if v in net.currency:
                continue
            within = sum(1 for x in net.adj[v] if x in inside)
            # a species that links at least two reactions of the group comes before any reaction
            rank = 0 if net.is_species[v] and within >= 2 else 1 if not net.is_species[v] else 2
            hubs.append((rank, -within, -within / max(net.degree(v), 1), net.ids[v], v))
        for _rank, _w, _share, _id, v in sorted(hubs):
            label = net.rows[net.ids[v]]["name"].strip()
            if label and label not in self.used_labels:
                self.used_labels.add(label)
                return label, "member", net.ids[v]
        return "", "none", ""


# ---- the strategy --------------------------------------------------------------------------------------
def run(features, algorithm="leiden", seed=1):
    features = pathlib.Path(features)
    net = Net(read_tsv(features / "nodes.tsv"), read_tsv(features / "edges.tsv"))
    net.find_currency()
    core = [i for i in range(len(net.ids)) if i not in net.currency]

    # R4-R7: tree of groups, each a dict(members, children, level, parent)
    groups = []
    rule = {}

    def grow(members, level, parent):
        parts, part_rule, raw = partition(net, members, algorithm, seed)
        if level > 1 and len(parts) < 2:
            return None
        made = []
        for p in sorted(parts, key=lambda p: (-len(p), p[0])):
            g = dict(members=p, children=[], level=level, parent=parent)
            groups.append(g)
            made.append(g)
            for v in p:
                rule[v] = part_rule[v] if level == 1 else ("R7" if part_rule[v] == "R4" else part_rule[v])
        for g in made:
            if len(g["members"]) > S_MAX and level < MAX_LEVEL:
                kids = grow(g["members"], level + 1, g)
                if kids:
                    g["children"] = kids
        return made, raw

    top, raw_top = grow(core, 1, None)

    # ids: M01, M01.02, M01.02.03 - larger groups first
    def number(siblings, prefix):
        width = max(2, len(str(len(siblings))))
        for k, g in enumerate(siblings, start=1):
            g["id"] = f"{prefix}{k:0{width}d}"
            if g["children"]:
                number(g["children"][0], g["id"] + ".")

    number(top, "M")
    for g in groups:
        g["children"] = g["children"][0] if g["children"] else []

    leaf_of = {}
    chain = {}
    for g in groups:
        if not g["children"]:
            path = [g]
            while path[-1]["parent"] is not None:
                path.append(path[-1]["parent"])
            ids = [x["id"] for x in reversed(path)]
            for v in g["members"]:
                leaf_of[v] = g
                chain[v] = ids

    # R8: currency species
    clones = {}
    for s in sorted(net.currency, key=lambda i: net.ids[i]):
        count = collections.Counter(leaf_of[r]["id"] for r in net.adj[s] if r in leaf_of)
        if count:
            home = sorted(count, key=lambda gid: (-count[gid], gid))[0]
        else:  # a currency species always has reactions; kept as a guard
            home = top[0]["id"]
        g = next(x for x in groups if x["id"] == home)
        path = [g]
        while path[-1]["parent"] is not None:
            path.append(path[-1]["parent"])
        chain[s] = [x["id"] for x in reversed(path)]
        for x in path:
            x.setdefault("extra", []).append(s)
        clones[s] = sorted(gid for gid in count if gid != home)
        rule[s] = "R8"

    # R9: names, level by level, larger groups first
    namer = Namer(net)
    meta = []
    for g in sorted(groups, key=lambda g: (g["level"], -len(g["members"]), g["id"])):
        members = sorted(g["members"] + g.get("extra", []))
        label, source, xref = namer.name(members)
        g["label"], g["name_source"], g["name_xref"], g["size"] = label, source, xref, len(members)
        meta.append(dict(group_id=g["id"], level=g["level"], parent=g["parent"]["id"] if g["parent"] else "",
                         label=label, basis=BASIS[algorithm], name_source=source, name_xref=xref, compartment_id=""))

    rows = []
    for v, node_id in enumerate(net.ids):
        ids = chain[v] + [""] * (3 - len(chain[v]))
        rows.append(dict(node_id=node_id, L1=ids[0], L2=ids[1], L3=ids[2], rule=rule[v],
                         clone_in="|".join(clones.get(v, []))))
    rows.sort(key=lambda r: r["node_id"])
    meta.sort(key=lambda r: r["group_id"])

    leaves = [g for g in groups if not g["children"]]
    stats = dict(
        algorithm=algorithm, seed=seed, n_nodes=len(net.ids), n_links=len(net.links),
        degree_bound=net.fence, candidates=net.candidates, n_currency=len(net.currency),
        n_raw_modules=len(raw_top), raw_sizes=sorted((len(p) for p in raw_top), reverse=True),
        modularity_raw=round(modularity(net, core, raw_top), 4),
        modularity_L1=round(modularity(net, core, [g["members"] for g in top]), 4),
        modularity_leaves=round(modularity(net, core, [g["members"] for g in leaves]), 4),
        n_groups={f"L{k}": sum(1 for g in groups if g["level"] == k) for k in (1, 2, 3)},
        sizes={f"L{k}": sorted((g["size"] for g in groups if g["level"] == k), reverse=True) for k in (1, 2, 3)},
        leaf_sizes=sorted((g["size"] for g in leaves), reverse=True),
        rules=dict(sorted(collections.Counter(rule.values()).items())),
        name_sources=dict(sorted(collections.Counter(g["name_source"] for g in groups).items())),
        n_clone_species=sum(1 for s in clones if clones[s]), n_clones=sum(len(c) for c in clones.values()),
    )
    return rows, meta, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description="Topology grouping strategy (see the module docstring).")
    ap.add_argument("--features", required=True, help="folder with nodes.tsv, edges.tsv, compartments.tsv")
    ap.add_argument("--sbml", required=True, help="the SBML model (accepted for the contract; not opened)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    ap.add_argument("--algorithm", choices=ALGORITHMS, default="leiden")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--stats", default="", help="optional JSON file with summary numbers")
    a = ap.parse_args(argv)
    rows, meta, stats = run(a.features, a.algorithm, a.seed)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", rows, GROUP_COLUMNS)
    write_tsv(out / "groups_meta.tsv", meta, META_COLUMNS)
    if a.stats:
        pathlib.Path(a.stats).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(a.stats).write_text(json.dumps(stats, indent=1, sort_keys=True), encoding="utf-8", newline="\n")
    n = stats["n_groups"]
    print(f"{a.algorithm}: {len(rows)} nodes, {stats['n_currency']} currency species, "
          f"{n['L1']} L1 / {n['L2']} L2 / {n['L3']} L3 groups -> {out}")


if __name__ == "__main__":
    main()
