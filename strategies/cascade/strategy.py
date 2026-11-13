"""Run (one model; paths relative to the unpacked bundle, this file at strategies/cascade/strategy.py):
    python strategies/cascade/strategy.py --features panel/features/skm_ja --sbml panel/models/skm_ja.xml --out strategies/cascade/out/skm_ja
    python common/contract.py validate panel/features/skm_ja strategies/cascade/out/skm_ja
    python common/contract.py write-sbml panel/models/skm_ja.xml strategies/cascade/out/skm_ja --features panel/features/skm_ja --strategy cascade -o strategies/cascade/out/skm_ja/skm_ja_grouped.xml
All six models (PowerShell):
    foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") { python strategies/cascade/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/cascade/out/$m }
All six models (bash):
    for m in brenda_merged e_coli_core R-HSA-77289 skm_pss skm_ja skm_sa; do python strategies/cascade/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/cascade/out/$m; done

Cascade: flow-anchored modules for an SBML model (BioHackathon Europe 2026, project 22, strategy "cascade").
Signalling content is grouped around anchors (sensed small molecules, perturbing agents): each anchor gets the
steps that produce it and the steps downstream of it. Metabolic content without an anchor is cut into routes
joined at weak branch points and separated at hubs. Nothing is read but the feature tables of the model.

Rules, in the order applied (full text in RULES.md)
    R1  class of every reaction: transport, genetic, binding, metabolic, state
    R2  currency species: carrier pairs and hubs, found from the data
    R3  anchors: sensed small molecules and perturbing agents
    R4  anchors that share targets form one family
    R5  production core of a family: metabolic and transport steps on a mass-flow path through the anchor
    R6  catalyst supply: steps that make, move or activate a modifier of the production core
    R7  downstream steps go to the family that sends them the largest share of flow
    R8  side inputs: steps that feed a placed step
    R9  a family with fewer than MIN_RX reactions is not an anchor; R4-R8 are repeated without it
    R10 routes: remaining metabolic and transport steps are linked through shared non-currency species,
        weakest branch points first, up to CAP reactions per route
    R11 a route below MIN_RX reactions joins the route it shares a branch point with
    R12 collateral steps join the nearest placed step of the undirected flow graph
    R13 the currency species and the steps between them form the currency module
    R14 steps tied to nothing: one group per loc_space
    R15 level 2 of a family: stage (production, molecular interaction, genetic production)
    R16 level 2 of a route: its cycle and its linear branches
    R17 species: where it is made, else where it is used (anchors and currency species have their own home)
    R18 clones: one copy of a species in every other top-level group that uses it
    R19 a top-level group holding more than half of the model is replaced by its level-2 groups
    R20 innermost level: loc_space, where a flow group spans several spaces
    R21 names: anchor (member), stage (SBO term), most characteristic species of a route (member), compartment

Only the Python standard library is used. --sbml is accepted for the common CLI; every value the rules read
is in nodes.tsv / edges.tsv / compartments.tsv, so the file itself is not opened.
"""
from __future__ import annotations

import argparse
import collections
import csv
import fractions
import heapq
import math
import pathlib

# ---------------------------------------------------------------- parameters (RULES.md, parameter table)
MIN_RX = 3            # a flow module, a stage or a route element needs at least 3 reactions
MIN_LOC = 5           # a localisation sub-group needs at least 5 nodes
PAIR_MIN = 3          # a carrier pair faces each other across at least 3 reactions
FAMILY_OVERLAP = 0.5  # anchors sharing at least half of the smaller target set are one family
CATCH_ALL = 0.5       # a top-level group above this share of all nodes is replaced by its children
CAP_FACTOR = 1.0      # route cap = ceil(CAP_FACTOR * sqrt(number of reactions))
# SBO terms the rules read
SBO_SMALL = "SBO:0000247"        # simple chemical
SBO_STIMULUS = "SBO:0000405"     # perturbing agent
SBO_TRANSPORT = ("SBO:0000655", "SBO:0000185", "SBO:0000627", "SBO:0000628", "SBO:0000632")
SBO_GENETIC = ("SBO:0000589", "SBO:0000205")
SBO_BINDING = ("SBO:0000177",)
SIGNALLING = ("binding", "state", "genetic")
# stage vocabulary: key -> (SBO id, SBO name), as in the SBO release in common/SBO_OBO.obo
STAGE = {"s1": ("SBO:0000393", "production"),
         "s2": ("SBO:0000344", "molecular interaction"),
         "s3": ("SBO:0000589", "genetic production")}
GROUP_COLUMNS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
META_COLUMNS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def top(counter):
    """Key with the highest count; ties go to the smallest key."""
    return sorted(counter, key=lambda k: (-counter[k], k))[0]


# ---------------------------------------------------------------- model
class Model:
    def __init__(self, features_dir):
        f = pathlib.Path(features_dir)
        self.node = {r["node_id"]: r for r in read_tsv(f / "nodes.tsv")}
        self.comp = {r["compartment_id"]: r for r in read_tsv(f / "compartments.tsv")}
        self.species = sorted(n for n, r in self.node.items() if r["kind"] == "species")
        self.reactions = sorted(n for n, r in self.node.items() if r["kind"] == "reaction")
        self.part = {r: {"reactant": [], "product": [], "modifier": []} for r in self.reactions}
        for e in read_tsv(f / "edges.tsv"):
            if e["reaction_id"] in self.part and e["species_id"] in self.node and e["role"] in self.part[e["reaction_id"]]:
                self.part[e["reaction_id"]][e["role"]].append(e["species_id"])
        for r in self.reactions:
            for role in self.part[r]:
                self.part[r][role] = sorted(set(self.part[r][role]))
        self.small = {s for s in self.species if self._is_small(s)}
        self.stimulus = {s for s in self.species if SBO_STIMULUS in self.node[s]["sbo_path"].split(">")}
        self.uses = collections.defaultdict(list)             # species -> reactions it takes part in
        for r in self.reactions:
            for s in sorted(set(self.part[r]["reactant"]) | set(self.part[r]["product"]) | set(self.part[r]["modifier"])):
                self.uses[s].append(r)

    def _is_small(self, s):
        path = self.node[s]["sbo_path"].split(">") if self.node[s]["sbo_path"] else []
        return not path or SBO_SMALL in path          # fallback: a species without SBO term counts as a small molecule

    def name(self, n):
        return self.node[n]["name"] or n

    def sides(self, r):
        return self.part[r]["reactant"], self.part[r]["product"]


# ---------------------------------------------------------------- R1 reaction class
def classify(m):
    """transport | genetic | binding | metabolic | state. binding, state and genetic are the signalling classes."""
    out = {}
    for r in m.reactions:
        row = m.node[r]
        path = row["sbo_path"].split(">") if row["sbo_path"] else []
        a, b = m.sides(r)
        if any(t in path for t in SBO_TRANSPORT):
            out[r] = "transport"
        elif any(t in path for t in SBO_GENETIC):
            out[r] = "genetic"
        elif any(t in path for t in SBO_BINDING):
            out[r] = "binding"
        elif (a or b) and all(s in m.small for s in a + b):
            out[r] = "metabolic"
        elif not path and row["loc_kind"] != "inside" and len(a) == 1 and len(b) == 1:
            out[r] = "transport"                      # fallback without SBO: one species moved between spaces
        else:
            out[r] = "state"
    return out


# ---------------------------------------------------------------- R2 currency
def carrier_pairs(m):
    """Carrier pairs: two species on opposite sides of at least PAIR_MIN reactions, in more than half of the
    reactions of each, and - taken in order of that count - such that more than half of those reactions keep
    a reactant and a product when the pair and the carriers accepted before it are left out."""
    deg = collections.Counter()
    opp = collections.defaultdict(list)
    for r in m.reactions:
        a, b = m.sides(r)
        for s in sorted(set(a) | set(b)):
            deg[s] += 1
        for x in a:
            for y in b:
                if x != y and x not in b and y not in a:
                    opp[(min(x, y), max(x, y))].append(r)
    cand = [(len(rs), pair) for pair, rs in opp.items()
            if len(rs) >= PAIR_MIN and 2 * len(rs) > deg[pair[0]] and 2 * len(rs) > deg[pair[1]]]
    carriers, accepted = set(), []
    for n, (x, y) in sorted(cand, key=lambda t: (-t[0], t[1])):
        ok = 0
        for r in opp[(x, y)]:
            a, b = m.sides(r)
            if (set(a) - carriers - {x, y}) and (set(b) - carriers - {x, y}):
                ok += 1
        if 2 * ok > n:
            carriers |= {x, y}
            accepted.append((x, y, n))
    return carriers, deg, accepted


# ---------------------------------------------------------------- flow graphs
def build_graphs(m, rclass, currency):
    """Directed graphs over species and reactions (currency species left out).

    mass:    reactant -> reaction -> product, metabolic and transport reactions only.
    supply:  reactant -> reaction -> product and modifier -> reaction, all reactions (what a step needs).
    flow:    supply, and for signalling reactions also reaction -> reactant (a controlled binding, state
             change or degradation changes the amount of its substrate, so the signal passes to it).
    """
    mass, supply, flow = (collections.defaultdict(set) for _ in range(3))
    for r in m.reactions:
        a, b = m.sides(r)
        for s in a:
            if s in currency:
                continue
            supply[s].add(r)
            flow[s].add(r)
            if rclass[r] in SIGNALLING:
                flow[r].add(s)
            if rclass[r] in ("metabolic", "transport"):
                mass[s].add(r)
        for s in b:
            if s in currency:
                continue
            supply[r].add(s)
            flow[r].add(s)
            if rclass[r] in ("metabolic", "transport"):
                mass[r].add(s)
        for s in m.part[r]["modifier"]:
            if s not in currency:
                supply[s].add(r)
                flow[s].add(r)
    return mass, supply, flow


def reverse(adj):
    out = collections.defaultdict(set)
    for u in adj:
        for v in adj[u]:
            out[v].add(u)
    return out


def undirected(adj):
    out = collections.defaultdict(set)
    for u in adj:
        for v in adj[u]:
            out[u].add(v)
            out[v].add(u)
    return out


def spread(adj, sources, order):
    """Flow-share search. sources: {node: label}; order: {label: rank}.

    Leaving a node with k onward edges costs log(k): flow splits evenly, so the cost of a path is -log of the
    share of the source's flow that arrives. Ties: fewer steps, then the label of better rank, then node id.
    Returns {node: (cost, steps, label)} and the predecessor of every node.
    """
    done, pred = {}, {}
    heap = [(0.0, 0, order[lab], n, lab, "") for n, lab in sources.items()]
    heapq.heapify(heap)
    while heap:
        cost, steps, rk, u, lab, p = heapq.heappop(heap)
        if u in done:
            continue
        done[u] = (cost, steps, lab)
        pred[u] = p
        outs = [v for v in sorted(adj.get(u, ())) if v != p]
        if not outs:
            continue
        w = math.log(len(outs))
        for v in outs:
            if v not in done:
                heapq.heappush(heap, (round(cost + w, 9), steps + 1, rk, v, lab, u))
    return done, pred


def hops(adj, sources):
    """Breadth-first distance from a set of nodes."""
    dist = {s: 0 for s in sources}
    queue = collections.deque(sorted(sources))
    while queue:
        u = queue.popleft()
        for v in sorted(adj.get(u, ())):
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    return dist


# ---------------------------------------------------------------- R3 anchors, R4 families
def find_anchors(m, rclass, currency):
    """{anchor: its sensing reactions}. An anchor is a small molecule that is a modifier of a signalling
    reaction (binding, state, genetic) or a reactant of a binding reaction, or a species with SBO
    'perturbing agent' that takes part in a reaction. Currency species are never anchors."""
    sense = collections.defaultdict(set)
    for r in m.reactions:
        if rclass[r] in SIGNALLING:
            for s in m.part[r]["modifier"] + (m.part[r]["reactant"] if rclass[r] == "binding" else []):
                sense[s].add(r)
    out = {}
    for s in m.species:
        if s in currency:
            continue
        if s in m.stimulus and m.uses.get(s):
            out[s] = set(sense.get(s, ()))
        elif s in m.small and sense.get(s):
            out[s] = set(sense[s])
    return out


def anchor_families(m, anchors):
    """Anchors whose sensing reactions share targets are one family (e.g. several ligands of one receptor).
    Returns {primary anchor: [anchors]}; the primary anchor has the most sensing reactions (ties: id)."""
    targets = {}
    for a, rs in anchors.items():
        t = set()
        for r in rs:
            t |= set(m.part[r]["reactant"]) | set(m.part[r]["product"])
        targets[a] = t - set(anchors)
    parent = {a: a for a in anchors}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    ids = sorted(anchors)
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            ta, tb = targets[a], targets[b]
            if ta and tb and len(ta & tb) >= FAMILY_OVERLAP * min(len(ta), len(tb)):
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    fam = collections.defaultdict(list)
    for a in ids:
        fam[find(a)].append(a)
    out = {}
    for members in fam.values():
        members = sorted(members, key=lambda a: (-len(anchors[a]), a))
        out[members[0]] = members
    return out


# ---------------------------------------------------------------- R5-R8 cascade
def default_stage(m, rclass, r):
    """Stage of a step placed by R8 or R12: by what the step does."""
    if rclass[r] == "genetic":
        return "s3"
    a, b = m.sides(r)
    if rclass[r] == "metabolic" or (rclass[r] == "transport" and all(s in m.small for s in a + b)):
        return "s1"
    return "s2"


def cascade(m, rclass, currency, anchors, families):
    """{reaction: (family, stage, rule)} for every reaction tied to an anchor family by directed flow."""
    fam_of = {a: f for f, members in families.items() for a in members}
    rank_list = sorted(families, key=lambda f: (-sum(len(anchors[a]) for a in families[f]), f))
    order = {f: i for i, f in enumerate(rank_list)}
    mass, supply, flow = build_graphs(m, rclass, currency)
    rmass, rsupply = reverse(mass), reverse(supply)
    rx = set(m.reactions)
    placed = {}

    # R5 production core: metabolic/transport steps on a mass-flow path into or out of the anchor
    ligands = {a: fam_of[a] for a in anchors if a in m.small}
    best = {}
    for a in sorted(ligands):
        for dist in (hops(rmass, [a]), hops(mass, [a])):
            for n, d in dist.items():
                if n in rx:
                    key = (d, order[ligands[a]], a)
                    if n not in best or key < best[n][0]:
                        best[n] = (key, ligands[a])
    for r in sorted(best):
        placed[r] = (best[r][1], "s1", "R5")

    # R6 catalyst supply: steps on a reactant -> product path into a modifier of the production core
    # (the enzyme is made, moved or activated). The path stops at anchors and at production cores.
    best = {}
    for f in rank_list:
        start = sorted({s for r, (g, _s, _rule) in placed.items() if g == f for s in m.part[r]["modifier"]
                        if s not in currency and s not in anchors})
        dist, queue = {s: 0 for s in start}, collections.deque(start)
        while queue:
            s = queue.popleft()
            for r in sorted(rsupply.get(s, ())):
                if r in placed or r in dist or s not in m.part[r]["product"]:
                    continue
                if any(x in anchors for x in m.part[r]["reactant"]):
                    continue
                dist[r] = dist[s] + 1
                for x in m.part[r]["reactant"]:
                    if x not in dist and x not in currency:
                        dist[x] = dist[r] + 1
                        queue.append(x)
        for n, d in dist.items():
            if n in rx:
                key = (d, order[f])
                if n not in best or key < best[n][0]:
                    best[n] = (key, f)
    for r in sorted(best):
        placed[r] = (best[r][1], "s1", "R6")

    # R7 downstream by flow share
    src = {a: fam_of[a] for a in anchors}
    fwd, fpred = spread(flow, src, order)
    genetic_seen = {}
    for n in sorted(fwd, key=lambda n: (fwd[n][0], fwd[n][1], n)):
        p = fpred[n]
        genetic_seen[n] = (n in rx and rclass[n] == "genetic") or (bool(p) and genetic_seen.get(p, False))
    for r in m.reactions:
        if r not in placed and r in fwd:
            placed[r] = (fwd[r][2], "s3" if genetic_seen[r] else "s2", "R7")

    # R8 side inputs: steps from which a placed step is reached; they take its family
    side, _ = spread(rsupply, {r: f for r, (f, _s, _rule) in placed.items()}, order)
    for r in m.reactions:
        if r not in placed and r in side:
            placed[r] = (side[r][2], default_stage(m, rclass, r), "R8")
    return placed, rank_list


# ---------------------------------------------------------------- R10-R11 routes
class UnionFind:
    def __init__(self, items):
        self.parent = {i: i for i in items}
        self.size = {i: 1 for i in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return a
        if b < a:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]
        return a


def link_species(m, rclass, r):
    """Species through which a reaction can be linked: reactants and products; for a signalling reaction
    also its modifiers."""
    a, b = m.sides(r)
    return sorted(set(a) | set(b) | (set(m.part[r]["modifier"]) if rclass[r] in SIGNALLING else set()))


def routes(m, rclass, currency, rset, cap):
    """Route modules among the reactions rset: {reaction: module root} and the rule that placed each."""
    inc = collections.defaultdict(list)
    for r in sorted(rset):
        for s in link_species(m, rclass, r):
            if s not in currency:
                inc[s].append(r)
    uf = UnionFind(sorted(rset))
    for s in sorted(inc, key=lambda s: (len(inc[s]), s)):      # R10: weakest branch points first
        roots = sorted({uf.find(r) for r in inc[s]})
        if len(roots) > 1 and sum(uf.size[x] for x in roots) <= cap:
            for x in roots[1:]:
                uf.union(roots[0], x)
    mod = {r: uf.find(r) for r in sorted(rset)}
    rule = {r: "R10" for r in mod}
    while True:                                                 # R11: a route below MIN_RX joins a neighbour
        size = collections.Counter(mod.values())
        moved = False
        for k in sorted(size, key=lambda k: (size[k], k)):
            if size[k] >= MIN_RX:
                break
            members = sorted(r for r in mod if mod[r] == k)
            shared = sorted({s for r in members for s in link_species(m, rclass, r) if s not in currency},
                            key=lambda s: (len(inc[s]), s))
            for s in shared:
                cnt = collections.Counter(mod[r] for r in inc[s] if mod[r] != k)
                if cnt:
                    for r in members:
                        mod[r], rule[r] = top(cnt), "R11"
                    moved = True
                    break
            if moved:
                break
        if not moved:
            break
    size = collections.Counter(mod.values())
    return {r: k for r, k in mod.items() if size[k] >= MIN_RX}, rule


# ---------------------------------------------------------------- R16 route elements
def cycle_reactions(m, currency, members):
    """Reactions of a module that lie on a directed cycle through at least three species (reactant -> product
    as written; a reaction with reversible=1 may be passed either way)."""
    edges = collections.defaultdict(set)          # species -> {(species, reaction)}
    for r in members:
        a, b = m.sides(r)
        for x in a:
            for y in b:
                if x == y or x in currency or y in currency:
                    continue
                edges[x].add((y, r))
                if m.node[r]["reversible"] == "1":
                    edges[y].add((x, r))
    on_cycle = set()
    for u in sorted(edges):
        for v, r in sorted(edges[u]):
            if r in on_cycle:
                continue
            # a path v -> ... -> u that avoids reaction r and has an inner species closes such a cycle
            seen, queue, found = {v}, collections.deque([v]), False
            while queue and not found:
                x = queue.popleft()
                for y, q in sorted(edges.get(x, ())):
                    if q == r or y in seen:
                        continue
                    if y == u:
                        if x != v:
                            found = True
                            break
                        continue
                    seen.add(y)
                    queue.append(y)
            if found:
                on_cycle.add(r)
    return on_cycle


def route_elements(m, rclass, currency, members):
    """Split one route module into its cycle and its linear branches: {reaction: element key or ''}.
    No split when the cycle or the rest has fewer than MIN_RX reactions."""
    core = cycle_reactions(m, currency, members)
    out = {r: "" for r in members}
    if len(core) < MIN_RX or len(members) - len(core) < MIN_RX:
        return out
    inc = collections.defaultdict(list)
    for r in members:
        if r not in core:
            for s in link_species(m, rclass, r):
                if s not in currency:
                    inc[s].append(r)
    core_species = {s for r in core for s in link_species(m, rclass, r)}
    uf = UnionFind(sorted(r for r in members if r not in core))
    for s in sorted(inc):                          # branches: pieces outside the cycle, joined by species off it
        if s not in core_species:
            for r in inc[s][1:]:
                uf.union(inc[s][0], r)
    pieces = collections.defaultdict(list)
    for r in sorted(uf.parent):
        pieces[uf.find(r)].append(r)
    for r in core:
        out[r] = "c"
    big = [p for p in pieces.values() if len(p) >= MIN_RX]
    for i, piece in enumerate(sorted(big, key=lambda p: (-len(p), p[0])), start=1):
        for r in piece:
            out[r] = f"b{i}"
    return out


# ---------------------------------------------------------------- R21 naming helpers
def characteristic(m, reactions, species):
    """Species of a group ordered by how much they characterise it: n^3 / u^2, where n is the number of the
    group's reactions the species is a reactant or product of and u the number of all reactions it takes
    part in (many reactions here, few elsewhere). Ties: id. Species with the same name are listed once."""
    species = set(species)
    inside = collections.Counter(s for r in reactions for s in sorted(set(m.sides(r)[0]) | set(m.sides(r)[1]))
                                 if s in species)
    ranked = sorted(inside, key=lambda s: (-fractions.Fraction(inside[s] ** 3, len(m.uses[s]) ** 2), s))
    out, names = [], set()
    for s in ranked:
        if m.name(s) not in names:
            out.append(s)
            names.add(m.name(s))
    return out


# ---------------------------------------------------------------- the procedure
def run(features_dir):
    m = Model(features_dir)
    rclass = classify(m)                                                    # R1
    carriers, degree, pairs = carrier_pairs(m)                              # R2
    cap = max(MIN_RX, math.ceil(CAP_FACTOR * math.sqrt(len(m.reactions)))) if m.reactions else MIN_RX
    hubs = {s for s in m.species if degree[s] > cap and s in m.small}
    currency = carriers | hubs
    while True:      # the copy of a currency species across a transport step is currency too
        extra = set()
        for r in m.reactions:
            if rclass[r] == "transport":
                left = [s for s in m.sides(r)[0] + m.sides(r)[1] if s not in currency]
                if len(left) == 1 and len(m.sides(r)[0] + m.sides(r)[1]) > 1:
                    extra.add(left[0])
        if not extra:
            break
        currency |= extra

    anchors = find_anchors(m, rclass, currency)                             # R3
    while True:                                                             # R4-R8, R9
        families = anchor_families(m, anchors)
        placed, rank_list = cascade(m, rclass, currency, anchors, families) if anchors else ({}, [])
        size = collections.Counter(f for f, _s, _r in placed.values())
        weak = [f for f in families if size[f] < MIN_RX]
        if not weak:
            break
        for f in weak:
            for a in families[f]:
                del anchors[a]
    fam_of = {a: f for f, members in families.items() for a in members}

    # reaction -> [module key, sub key, rule]
    where = {r: [("C", f), s, rule] for r, (f, s, rule) in placed.items()}

    # R10-R11 routes among the metabolic and transport reactions no anchor reached
    rest = [r for r in m.reactions if r not in where and rclass[r] in ("metabolic", "transport")]
    rmod, rrule = routes(m, rclass, currency, rest, cap)
    for r, k in rmod.items():
        where[r] = [("P", k), "", rrule[r]]

    # R12 collateral: nearest placed step in the undirected flow graph
    _mass, _supply, flow = build_graphs(m, rclass, currency)
    uflow = undirected(flow)

    def collateral():
        rest = [r for r in m.reactions if r not in where]
        if not rest or not where:
            return
        mods_now = sorted({w[0] for w in where.values()})
        near, _ = spread(uflow, {r: w[0] for r, w in where.items()}, {k: i for i, k in enumerate(mods_now)})
        for r in rest:
            if r in near:
                k = near[r][2]
                where[r] = [k, default_stage(m, rclass, r) if k[0] == "C" else "", "R12"]

    collateral()
    # R10-R11 again for steps of any class that are still unplaced (fallback for content without anchors),
    # then R12 again for what touches these routes
    rest = [r for r in m.reactions if r not in where]
    rmod, rrule = routes(m, rclass, currency, rest, cap)
    for r, k in rmod.items():
        where[r] = [("P", k), "", rrule[r]]
    collateral()
    # R13 currency module: the currency species and the steps whose reactants and products are all currency
    rest = [r for r in m.reactions if r not in where]
    only_currency = [r for r in rest if (m.sides(r)[0] or m.sides(r)[1]) and all(s in currency for s in m.sides(r)[0] + m.sides(r)[1])]
    has_k = len(currency) + len(only_currency) >= MIN_LOC
    if has_k:
        for r in only_currency:
            where[r] = [("K", ""), "", "R13"]
    # R14 unconnected steps
    for r in m.reactions:
        if r not in where:
            where[r] = [("U", m.node[r]["loc_space"]), "", "R14"]

    # R15 stages: a stage below MIN_RX reactions merges into the neighbouring stage
    for f in sorted({w[0] for w in where.values() if w[0][0] == "C"}):
        mine = [r for r in m.reactions if where[r][0] == f]
        while True:
            cnt = collections.Counter(where[r][1] for r in mine)
            small = sorted(s for s in cnt if cnt[s] < MIN_RX)
            if not small or len(cnt) == 1:
                break
            prefer = {"s1": ("s2", "s3"), "s2": ("s3", "s1"), "s3": ("s2", "s1")}[small[0]]
            target = next(t for t in prefer if t in cnt)
            for r in mine:
                if where[r][1] == small[0]:
                    where[r][1] = target
        if len({where[r][1] for r in mine}) == 1:
            for r in mine:
                where[r][1] = ""
    # R16 route elements
    for k in sorted({w[0] for w in where.values() if w[0][0] == "P"}):
        mine = sorted(r for r in m.reactions if where[r][0] == k)
        for r, e in route_elements(m, rclass, currency, mine).items():
            where[r][1] = e

    # R17 species
    made, used = collections.defaultdict(list), collections.defaultdict(list)
    for r in m.reactions:
        g = (where[r][0], where[r][1])
        for s in m.part[r]["product"]:
            made[s].append(g)
        for s in m.part[r]["reactant"] + m.part[r]["modifier"]:
            used[s].append(g)
    home, srule = {}, {}
    for s in m.species:
        both = made[s] + used[s]
        own = [g for g in both if s in fam_of and g[0] == ("C", fam_of[s])]
        if own:                                                 # R17a an anchor is homed in its own family
            home[s], srule[s] = top(collections.Counter([g for g in made[s] if g in own] or own)), "R17a"
        elif has_k and s in currency:                           # R17b a currency species sits in the currency module
            home[s], srule[s] = (("K", ""), ""), "R17b"
        elif made[s]:                                           # R17c where it is made
            home[s], srule[s] = top(collections.Counter(made[s])), "R17c"
        elif used[s]:                                           # R17d where it is used
            home[s], srule[s] = top(collections.Counter(used[s])), "R17d"
        else:                                                   # in no reaction
            home[s], srule[s] = (("U", m.node[s]["loc_space"]), ""), "R14"
    for r in m.reactions:
        home[r] = (where[r][0], where[r][1])
    for n in sorted(home):                                      # R14: a node there sits in its own loc_space
        if home[n][0][0] == "U":
            home[n] = (("U", m.node[n]["loc_space"]), "")

    # ------------------------------------------------------------ group ids and members
    fam_rank = {f: i for i, f in enumerate(rank_list)}
    msize = collections.Counter(h[0] for h in home.values())
    mods = sorted(msize)
    c_mods = sorted((k for k in mods if k[0] == "C"), key=lambda k: fam_rank[k[1]])
    p_mods = sorted((k for k in mods if k[0] == "P"), key=lambda k: (-msize[k], k[1]))
    k_mods = [k for k in mods if k[0] == "K"]
    u_mods = sorted(k for k in mods if k[0] == "U")
    gid = {}
    for i, k in enumerate(c_mods, start=1):
        gid[k] = f"C{i:02d}"
    for i, k in enumerate(p_mods, start=1):
        gid[k] = f"P{i:02d}"
    for k in k_mods:
        gid[k] = "K01"
    for k in u_mods:
        gid[k] = "U_" + k[1]
    path, members = {}, collections.defaultdict(list)
    for n in sorted(home):
        k, sub = home[n]
        path[n] = [gid[k]] + ([f"{gid[k]}.{sub}"] if sub else [])
        for g in path[n]:
            members[g].append(n)

    # ------------------------------------------------------------ R21 names of flow groups
    meta = {}
    for k in c_mods:
        g, fam = gid[k], families[k[1]]
        shown = [m.name(a) for a in fam[:3]] + ([f"+{len(fam) - 3}"] if len(fam) > 3 else [])
        meta[g] = dict(label=" / ".join(shown), basis="flow", name_source="member", name_xref=k[1])
        for sub in sorted({home[n][1] for n in members[g] if home[n][1]}):
            meta[f"{g}.{sub}"] = dict(label=STAGE[sub][1], basis="flow", name_source="sbo", name_xref=STAGE[sub][0],
                                      parent=g)
    for k in p_mods:
        g = gid[k]
        for x in [g] + sorted({f"{g}.{home[n][1]}" for n in members[g] if home[n][1]}):
            rx = [n for n in members[x] if m.node[n]["kind"] == "reaction"]
            sp = [n for n in members[x] if m.node[n]["kind"] == "species"]
            shown = characteristic(m, rx, [s for s in sp if s not in currency]) or characteristic(m, rx, sp)
            basis = "topology" if x.endswith(".c") else "flow"
            if shown and x.endswith(".c"):                      # a cycle: its most characteristic member species
                label = f"cycle through {m.name(shown[0])}"
            else:                                               # a route: its two most characteristic species
                label = " / ".join(m.name(s) for s in shown[:2])
            if shown:
                meta[x] = dict(label=label, basis=basis, name_source="member", name_xref=shown[0])
            else:                                               # fallback: no member species, unnamed
                meta[x] = dict(label="", basis=basis, name_source="none", name_xref="")
            if x != g:
                meta[x]["parent"] = g
    for k in k_mods:                                            # currency module: its busiest currency species
        g = gid[k]
        shown, names = [], []
        for s in sorted((n for n in members[g] if n in currency), key=lambda s: (-degree[s], s)):
            if m.name(s) not in names and len(names) < 3:
                shown.append(s)
                names.append(m.name(s))
        meta[g] = dict(label=" / ".join(names), basis="flow", name_source="member", name_xref=shown[0])
    for k in u_mods:
        meta[gid[k]] = dict(label=m.comp[k[1]]["name"], basis="localisation", name_source="compartment",
                            name_xref=k[1], compartment_id=k[1])

    # ------------------------------------------------------------ R19 catch-all
    total = len(home)
    for g in sorted(gid.values()):
        kids = sorted(x for x in meta if meta[x].get("parent") == g)
        if len(members[g]) > CATCH_ALL * total and len(kids) > 1:
            for x in kids:
                meta[x].pop("parent")
                meta[x]["label"] = f"{meta[g]['label']}: {meta[x]['label']}"
            for n in members[g]:
                if len(path[n]) > 1:
                    path[n] = path[n][1:]

    # ------------------------------------------------------------ R20 localisation level
    deepest = collections.defaultdict(list)
    for n in sorted(path):
        deepest[path[n][-1]].append(n)
    for g in sorted(deepest):
        if meta[g]["basis"] == "localisation":
            continue
        by_space = collections.defaultdict(list)
        for n in deepest[g]:
            by_space[m.node[n]["loc_space"]].append(n)
        keep = sorted(c for c, ns in by_space.items() if len(ns) >= MIN_LOC and c in m.comp)
        if len(keep) < 2:
            continue
        for c in keep:
            x = f"{g}.{c}"
            meta[x] = dict(label=m.comp[c]["name"], basis="localisation", name_source="compartment", name_xref=c,
                           compartment_id=c)
            for n in by_space[c]:
                path[n] = path[n] + [x]

    # ------------------------------------------------------------ R18 clones
    clones = {}
    for s in m.species:
        use = collections.defaultdict(collections.Counter)
        for r in m.uses[s]:
            if path[r][0] != path[s][0]:
                use[path[r][0]][path[r][-1]] += 1
        clones[s] = sorted(top(cnt) for cnt in use.values())

    # ------------------------------------------------------------ rows
    level, parent_of = {}, {}
    for n in sorted(path):
        for i, g in enumerate(path[n], start=1):
            level[g] = i
            parent_of[g] = path[n][i - 2] if i > 1 else ""
    rule_of = {r: where[r][2] for r in m.reactions}
    rule_of.update(srule)
    groups = []
    for n in sorted(path):
        chain = path[n] + [""] * (3 - len(path[n]))
        groups.append(dict(node_id=n, L1=chain[0], L2=chain[1], L3=chain[2], rule=rule_of[n],
                           clone_in="|".join(clones.get(n, []))))
    rows = [dict(group_id=g, level=level[g], parent=parent_of[g], label=meta[g]["label"], basis=meta[g]["basis"],
                 name_source=meta[g]["name_source"], name_xref=meta[g]["name_xref"],
                 compartment_id=meta[g].get("compartment_id", "")) for g in sorted(level)]
    info = dict(model=m, rclass=rclass, carriers=carriers, pairs=pairs, hubs=hubs, currency=currency, cap=cap,
                anchors=anchors, families=families, where=where, home=home, path=path)
    return groups, rows, info


def main(argv=None):
    ap = argparse.ArgumentParser(description="cascade: flow-anchored modules")
    ap.add_argument("--features", required=True, help="folder with nodes.tsv, edges.tsv, compartments.tsv")
    ap.add_argument("--sbml", default="", help="the SBML file (accepted for the common CLI, not opened)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    a = ap.parse_args(argv)
    groups, rows, _info = run(a.features)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", groups, GROUP_COLUMNS)
    write_tsv(out / "groups_meta.tsv", rows, META_COLUMNS)
    n1 = sum(r["level"] == 1 for r in rows)
    print(f"{pathlib.Path(a.features).name}: {len(groups)} nodes, {n1} top-level groups, {len(rows)} groups -> {out}")


if __name__ == "__main__":
    main()
