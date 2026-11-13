"""Run (from the bundle root, the folder that holds common/ and panel/):

    python hybrid/strategy.py --features panel/features/skm_ja --sbml panel/models/skm_ja.xml --out hybrid/out/skm_ja
    python common/contract.py validate panel/features/skm_ja hybrid/out/skm_ja
    python common/contract.py write-sbml panel/models/skm_ja.xml hybrid/out/skm_ja --features panel/features/skm_ja --strategy hybrid -o hybrid/out/skm_ja/skm_ja_grouped.xml

    # all six models (PowerShell):
    foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") {
      python hybrid/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out hybrid/out/$m
      python common/contract.py validate panel/features/$m hybrid/out/$m
      python common/contract.py write-sbml panel/models/$m.xml hybrid/out/$m --features panel/features/$m --strategy hybrid -o hybrid/out/$m/${m}_grouped.xml
    }

    # all six models (bash):
    for m in brenda_merged e_coli_core R-HSA-77289 skm_pss skm_ja skm_sa; do
      python hybrid/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out hybrid/out/$m
      python common/contract.py validate panel/features/$m hybrid/out/$m
      python common/contract.py write-sbml panel/models/$m.xml hybrid/out/$m --features panel/features/$m --strategy hybrid -o hybrid/out/$m/${m}_grouped.xml
    done

Strategy "hybrid": a size-bounded hierarchical grouping (rules R1..R15 of RULES.md).

    level 1   functional module    function terms of the SBML notes seed the modules where partners confirm them
                                   (R3); connectivity carries the term to the unannotated nodes (R4) and forms
                                   the modules where there is no annotation (R5); split (R6) and merge (R7) keep
                                   the number of modules and the share of the largest within bounds
    below     finer function (R9), reaction type from the SBO hierarchy (R10), localisation as given by
              nodes.tsv (R11), connectivity parts (R12, R13) - only where a group is above the leaf bound or
              spans several compartments; a level that gives one child or mostly singletons is skipped

Needs python-igraph and leidenalg (Leiden with a fixed seed). Reads only --features (nodes.tsv, edges.tsv,
compartments.tsv); --sbml is accepted because the contract fixes the command line, the file is not opened.
"""
from __future__ import annotations

import argparse
import collections
import csv
import html
import json
import math
import pathlib
import re
from fractions import Fraction

import igraph as ig
import leidenalg as la

# ---------------------------------------------------------------- parameters (RULES.md, parameter table)
LEAF_MAX = 40          # M: most nodes a leaf group may hold
GROUP_MIN = 3          # m: a group below this size is merged where a merge rule exists
K_MIN, K_MAX = 4, 12   # the number of level-1 groups aimed at is ceil(n / M), kept inside [K_MIN, K_MAX]
SHARE_CAP = Fraction(1, 3)   # no group above a third of its parent (level 1: of the model)
HUB_FENCE = 1.5        # Tukey fence: a species with degree > Q3 + 1.5 IQR is a hub and is cloned (R15)
SEED = 0               # Leiden random seed
RESOLUTION = 1.0       # Leiden resolution (plain modularity)
FUNCTION_NOTE_KEYS = ("mapman",)   # notes keys holding hierarchical function terms "<code>_<name path>"

BIN_RE = re.compile(r"(?:^|;)\s*(\d+(?:\.\d+)*)_")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in columns) + "\n")


def percentile(values, q):
    """Linear-interpolation percentile (the numpy default)."""
    v = sorted(values)
    if not v:
        return 0.0
    pos = (len(v) - 1) * q
    lo, hi = math.floor(pos), math.ceil(pos)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


class Model:
    """The tables of one model, the weighted species-reaction graph and the function terms."""

    def __init__(self, features_dir):
        f = pathlib.Path(features_dir)
        self.nodes = {r["node_id"]: r for r in read_tsv(f / "nodes.tsv")}
        self.ids = sorted(self.nodes)
        self.n = len(self.ids)
        self.species = [i for i in self.ids if self.nodes[i]["kind"] == "species"]
        self.comp = {r["compartment_id"]: r for r in read_tsv(f / "compartments.tsv")}
        self.adj = {i: set() for i in self.ids}
        self.roles = collections.defaultdict(set)          # (reaction, species) -> roles
        for e in read_tsv(f / "edges.tsv"):
            r, s = e["reaction_id"], e["species_id"]
            if r in self.adj and s in self.adj:
                self.adj[r].add(s)
                self.adj[s].add(r)
                self.roles[(r, s)].add(e["role"])
        # R2: an edge weighs 1 / (number of reactions of its species); exact fractions keep every sum
        # independent of the order of summation
        self.w = {}
        for s in self.species:
            for r in self.adj[s]:
                self.w[(s, r)] = self.w[(r, s)] = Fraction(1, len(self.adj[s]))
        deg = [len(self.adj[s]) for s in self.species]
        q1, q3 = percentile(deg, 0.25), percentile(deg, 0.75)
        self.hub_fence = q3 + HUB_FENCE * (q3 - q1)
        self.hubs = {s for s in self.species if len(self.adj[s]) > self.hub_fence}
        self.terms, self.term_name = self._function_terms()              # R1
        self.prefixes = {v: sorted({".".join(c.split(".")[:d + 1]) for c in codes for d in range(c.count(".") + 1)})
                         for v, codes in self.terms.items()}
        self.carriers = collections.defaultdict(set)       # term (any depth) -> nodes that carry it
        self.children = collections.defaultdict(set)       # term -> child terms that occur in the model
        for v, ps in self.prefixes.items():
            for t in ps:
                self.carriers[t].add(v)
                if "." in t:
                    self.children[t.rsplit(".", 1)[0]].add(t)

    # ------------------------------------------------------------ R1: function terms from the notes
    def _function_terms(self):
        """node -> sorted term codes; term code -> its name path, as far as the notes spell it out."""
        terms, names = {}, {}
        for i in self.ids:
            raw = self.nodes[i]["notes"]
            if not raw:
                continue
            try:
                notes = json.loads(raw)
            except ValueError:
                continue
            codes = []
            for key in FUNCTION_NOTE_KEYS:
                text = html.unescape(str(notes.get(key, "")))
                hits = list(BIN_RE.finditer(text))
                for j, hit in enumerate(hits):
                    end = hits[j + 1].start() if j + 1 < len(hits) else len(text)
                    code = hit.group(1)
                    path = [p.strip() for p in text[hit.end():end].strip().strip(";").split(".")]
                    parts = code.split(".")
                    codes.append(code)
                    if len(path) == len(parts):            # name path lines up with the code: all prefixes named
                        for d in range(len(parts)):
                            names.setdefault(".".join(parts[:d + 1]), path[:d + 1])
                    elif path and path[0]:
                        names.setdefault(parts[0], path[:1])
            if codes:
                terms[i] = sorted(set(codes))
        return terms, names

    # ------------------------------------------------------------ R2: connectivity communities
    def leiden(self, nodes, resolution=RESOLUTION):
        """Leiden communities (weighted modularity, fixed seed) of the subgraph on `nodes`: sorted lists."""
        nodes = sorted(nodes)
        ix = {v: i for i, v in enumerate(nodes)}
        edges, weights = [], []
        for v in nodes:
            if self.nodes[v]["kind"] != "species":
                continue
            for r in sorted(self.adj[v]):
                if r in ix:
                    edges.append((ix[v], ix[r]))
                    weights.append(float(self.w[(v, r)]))
        g = ig.Graph(n=len(nodes), edges=edges)
        g.es["weight"] = weights
        part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight",
                                 resolution_parameter=resolution, seed=SEED)
        out = [sorted(nodes[i] for i in members) for members in part]
        return sorted(out, key=lambda c: (-len(c), c[0]))

    # ------------------------------------------------------------ R12: connectivity parts within a size cap
    def parts(self, nodes, cap):
        """`nodes` cut into connectivity parts of at most `cap` nodes each, largest first: Leiden communities of
        the set (a community above the cap is cut again), merged by R7 down to the fewest parts the cap allows."""
        nodes = sorted(nodes)
        if len(nodes) <= cap:
            return [nodes]
        inside = set(nodes)
        lone = [v for v in nodes if not (self.adj[v] & inside)]          # no neighbour inside the set
        todo, done = [[v for v in nodes if self.adj[v] & inside]], []
        while todo:
            cur = todo.pop(0)
            if len(cur) <= cap:
                if cur:
                    done.append(cur)
                continue
            sub, res = self.leiden(cur), RESOLUTION
            while len(sub) < 2 and res < 64:             # one community: raise the resolution
                res *= 2
                sub = self.leiden(cur, res)
            if len(sub) < 2:                             # still one block: cut the sorted id list
                size = math.ceil(len(cur) / math.ceil(len(cur) / cap))
                sub = [cur[i:i + size] for i in range(0, len(cur), size)]
            todo.extend(sub)
        groups = {c[0]: set(c) for c in done}
        merge_groups(self, groups, inside, cap, target=math.ceil(len(nodes) / cap), term={})   # as few as fit
        if lone:                                         # unconnected nodes: join the smallest part with room,
            size = math.ceil(len(lone) / math.ceil(len(lone) / cap))     # else form parts in id order
            for i in range(0, len(lone), size):
                chunk = lone[i:i + size]
                room = sorted((g for g in groups if len(groups[g]) + len(chunk) <= cap),
                              key=lambda g: (len(groups[g]), g))
                if room and len(chunk) < GROUP_MIN:
                    groups[room[0]] |= set(chunk)
                else:
                    groups[chunk[0]] = set(chunk)
        return sorted((sorted(g) for g in groups.values()), key=lambda c: (-len(c), c[0]))


def related(a, b):
    """Number of leading code components two function terms share (0: unrelated or no term)."""
    if not a or not b:
        return 0
    pa, pb = a.split("."), b.split(".")
    d = 0
    while d < min(len(pa), len(pb)) and pa[d] == pb[d]:
        d += 1
    return d


def merge_groups(m: Model, groups, inside, cap, target, term, on_merge=None):
    """R7: while there are more than `target` groups or a group below GROUP_MIN, the smallest group joins a
    group it is connected to - the one with the most closely related function term, then the one with the
    highest modularity gain - as long as the union stays within `cap`. Without any connection: the highest
    gain alone. `groups` (key -> set of nodes) is changed in place; the larger group keeps its key."""
    where = {v: g for g, members in groups.items() for v in members}
    strength = {v: sum((m.w[(v, u)] for u in m.adj[v] if u in inside), Fraction(0)) for v in where}
    total = sum(strength.values(), Fraction(0)) or Fraction(1)
    mass = {g: sum((strength[v] for v in members), Fraction(0)) for g, members in groups.items()}
    blocked = set()
    while len(groups) > 1:
        order = sorted((g for g in groups if g not in blocked), key=lambda g: (len(groups[g]), g))
        too_many = target is not None and len(groups) > target
        pick = next((g for g in order if too_many or len(groups[g]) < GROUP_MIN), None)
        if pick is None:
            break
        link = collections.defaultdict(Fraction)
        for v in groups[pick]:
            for u in m.adj[v]:
                if u in where and where[u] != pick:
                    link[where[u]] += m.w[(v, u)]
        cands = [g for g in sorted(groups) if g != pick and len(groups[g]) + len(groups[pick]) <= cap]
        gain = {g: link.get(g, Fraction(0)) - mass[pick] * mass[g] / total for g in cands}
        connected = [g for g in cands if link.get(g, 0) > 0]
        if connected:
            best = min(connected, key=lambda g: (-related(term.get(pick), term.get(g)), -gain[g], g))
        elif cands:
            best = min(cands, key=lambda g: (-gain[g], g))
        else:
            blocked.add(pick)
            continue
        small, big = (pick, best) if (len(groups[pick]), best) <= (len(groups[best]), pick) else (best, pick)
        if on_merge:
            on_merge(small, big)
        for v in groups[small]:
            where[v] = big
        groups[big] |= groups.pop(small)
        mass[big] += mass.pop(small)
        blocked.discard(big)
    return groups


# ==================================================================== functional modules (R3 - R8)
def best_terms(m: Model, members, base=None):
    """R8: function terms that can name the node set, best first: [(F1, term, precision, recall, support)].

    A term (any depth of the term tree, below `base`) qualifies when at least GROUP_MIN members carry it and
    they are a majority of the annotated members. F1 compares the term's carriers in the model with the
    annotated members of the set.
    """
    annotated = [v for v in members if v in m.terms]
    count = collections.Counter()
    for v in annotated:
        for t in m.prefixes[v]:
            if base is None or t.startswith(base + "."):
                count[t] += 1
    out = []
    for t, g in count.items():
        precision, recall = Fraction(g, len(annotated)), Fraction(g, len(m.carriers[t]))
        if g >= GROUP_MIN and precision > Fraction(1, 2):
            out.append((2 * precision * recall / (precision + recall), t, precision, recall, g))
    return sorted(out, key=lambda x: (-x[0], -x[1].count("."), x[1]))


def function_split(m: Model, nodes, base, k_target, cap):
    """Functional modules of the node set `nodes`.

    base: the function term the set already carries (None at level 1); seeds use the terms below it.
    Returns (node -> group key, node -> rule id, group key -> naming record or None, info).
    """
    nodes = sorted(nodes)
    inside = set(nodes)
    depth0 = base.count(".") + 1 if base else 0
    annotated = [v for v in nodes if v in m.terms]
    communities = m.leiden(nodes)                                            # R2
    com_of = {v: c for c, members in enumerate(communities) for v in members}
    partners = {v: sorted(({u for r in m.adj[v] for u in m.adj[r]} & inside) - {v}) for v in annotated}
    expanded, topo_split = set(), set()

    def active_terms(v):
        """The terms of an annotated node at the current cut of the term tree (below `base`)."""
        out = set()
        for code in m.terms[v]:
            parts = code.split(".")
            if base is not None and (len(parts) <= depth0 or ".".join(parts[:depth0]) != base):
                continue
            d = depth0 + 1
            while d < len(parts) and ".".join(parts[:d]) in expanded:
                d += 1
            if ".".join(parts[:d]) not in expanded:
                out.add(".".join(parts[:d]))
        return out

    def assign():
        cand = {v: active_terms(v) for v in annotated}
        cand = {v: c for v, c in cand.items() if c}
        freq = collections.Counter(t for c in cand.values() for t in c)
        label, rule = {}, {}
        for v in sorted(cand):                                               # R3: seeds confirmed by partners
            mine = [u for u in partners[v] if u in cand]
            score = {}
            for t in cand[v]:
                observed = sum(1 for u in mine if t in cand[u])
                if observed:                                                 # observed minus expected partners
                    score[t] = observed - Fraction(len(mine) * freq[t], len(cand))
            if score:
                label[v], rule[v] = min(score, key=lambda t: (-score[t], freq[t], t)), "R3"
        frontier = True
        while frontier:                                                      # R4: spread inside a community
            frontier = {}
            for v in nodes:
                if v in label:
                    continue
                score = collections.defaultdict(Fraction)
                for u in m.adj[v]:
                    if u in label and com_of[u] == com_of[v]:
                        score[label[u]] += m.w[(v, u)]
                if score:
                    frontier[v] = min(score, key=lambda t: (-score[t], t))
            for v, t in frontier.items():
                label[v], rule[v] = t, "R4"
        for v in nodes:                                                      # R5: a community without a seed
            if v not in label:
                label[v], rule[v] = "~" + communities[com_of[v]][0], "R5"
        for key in sorted(topo_split):                                       # R6, second case
            members = [v for v in nodes if label[v] == key]
            if len(members) > cap:
                for part in m.parts(members, cap):
                    for v in part:
                        label[v], rule[v] = f"{key}~{part[0]}", "R6"
        return label, rule

    while True:                                                              # R6: split a group above the cap
        label, rule = assign()
        size = collections.Counter(label.values())
        step = False
        for g in sorted((g for g in size if size[g] > cap), key=lambda g: (-size[g], g)):
            if "~" not in g and m.children.get(g) and g not in expanded:
                expanded.add(g)                                              # first case: one level down the tree
                step = True
                break
            if g not in topo_split and "~" not in g[1:]:
                topo_split.add(g)
                step = True
                break
        if not step:
            break
    for v in nodes:                               # a seed that kept a term below a split term was placed by R6
        if rule[v] == "R3" and label[v].count(".") > depth0:
            rule[v] = "R6"
    groups = collections.defaultdict(set)
    for v, g in label.items():
        groups[g].add(v)
    groups = dict(groups)
    n_before = len(groups)
    origin = {g: {g} for g in groups}
    merges = []

    def on_merge(small, big):
        merges.append((small, big, len(groups[small])))
        for g in origin.pop(small):
            for v in nodes:
                if label[v] == g:
                    rule[v] = "R7"
        # the absorbed nodes keep `label` (their group before merging), so later merges still find them

    term = {g: (None if "~" in g else g) for g in groups}
    merge_groups(m, groups, inside, cap, k_target, term, on_merge)           # R7
    # R8: a term names one group only - the group it fits best
    options = sorted(((f1, t, p, r, g, key) for key in groups for f1, t, p, r, g in best_terms(m, groups[key], base)),
                     key=lambda x: (-x[0], -x[1].count("."), x[1], x[5]))
    naming, used = {}, set()
    for f1, t, p, r, g, key in options:
        if key not in naming and t not in used:
            naming[key] = dict(term=t, f1=round(float(f1), 3), precision=round(float(p), 3),
                               recall=round(float(r), 3), support=g)
            used.add(t)
    where = {v: g for g, members in groups.items() for v in members}
    info = dict(k_target=k_target, cap=cap, n_communities=len(communities), merges=merges, n_before_merge=n_before,
                expanded=sorted(expanded), n_seeds=sum(1 for v in rule if rule[v] in ("R3", "R6") and v in m.terms),
                n_annotated=len(annotated), naming={g: naming.get(g) for g in groups})
    return where, rule, {g: naming.get(g) for g in groups}, info


# ==================================================================== reaction type (R10), localisation (R11)
def home_by_reactions(m, species, key_of_reaction):
    """A species goes where most of its reactions are; ties: where it is a product, then the smaller key."""
    count, product = collections.Counter(), collections.Counter()
    for r in m.adj[species]:
        k = key_of_reaction.get(r)
        if k is None:
            continue
        count[k] += 1
        if "product" in m.roles[(r, species)]:
            product[k] += 1
    if not count:
        return None
    return min(count, key=lambda k: (-count[k], -product[k], k))


def type_split(m, members):
    """R10: split a group by the SBO term of its reactions. Returns {SBO term: [nodes]} or None (skipped)."""
    rx = sorted(v for v in members if m.nodes[v]["kind"] != "species")
    sp = sorted(v for v in members if m.nodes[v]["kind"] == "species")
    path = {r: [t for t in m.nodes[r]["sbo_path"].split(">") if t] for r in rx}
    cls = {r: path[r][-1] for r in rx if path[r]}
    if len(set(cls.values())) < 2:
        return None                                        # no SBO term, or one type only
    term_path = {}
    for p in path.values():
        for d in range(1, len(p) + 1):
            term_path.setdefault(p[d - 1], p[:d])

    def lca(a, b):
        pa, pb = term_path[a], term_path[b]
        d = 0
        while d < min(len(pa), len(pb)) and pa[d] == pb[d]:
            d += 1
        return pa[:d]

    for r in rx:                                           # fallback: a reaction without SBO term takes the
        if r not in cls:                                   # type most frequent among the reactions next to it
            votes = collections.Counter(cls[q] for s in m.adj[r] for q in m.adj[s] if q in cls and path[q])
            cls[r] = min(votes, key=lambda k: (-votes[k], k)) if votes else min(set(cls.values()))
    while True:
        classes = sorted(set(cls.values()))
        biggest = min(classes, key=lambda k: (-sum(1 for r in cls if cls[r] == k), k))
        home = {s: home_by_reactions(m, s, cls) or biggest for s in sp}
        size = collections.Counter(cls.values()) + collections.Counter(home.values())
        tiny = sorted((k for k in classes if size[k] < GROUP_MIN), key=lambda k: (size[k], k))
        if not tiny or len(classes) <= 2:
            break
        a = tiny[0]                                        # roll the smallest type up the SBO hierarchy
        b = min((k for k in classes if k != a), key=lambda k: (-len(lca(a, k)), -size[k], k))
        up = lca(a, b)[-1]
        for r in cls:
            if cls[r] in (a, b):
                cls[r] = up
    if tiny or len(classes) < 2:
        return None                                        # one type would be left with a crumb: skipped
    out = collections.defaultdict(list)
    for v, k in list(cls.items()) + list(home.items()):
        out[k].append(v)
    return {k: sorted(v) for k, v in out.items()}


def loc_key(m, v):
    nd = m.nodes[v]
    return nd["loc_space"] + ("@boundary" if nd["loc_kind"] == "boundary" else "")


def loc_split(m, members):
    """R11: split a group by the given localisation. Returns {compartment_id[@boundary]: [nodes]} or None."""
    key = {v: loc_key(m, v) for v in members}
    size = collections.Counter(key.values())
    for v in members:                                      # fewer than GROUP_MIN boundary nodes: they join the
        if key[v].endswith("@boundary") and size[key[v]] < GROUP_MIN:    # group of their space (contract rule 5)
            key[v] = m.nodes[v]["loc_space"]
    size = collections.Counter(key.values())
    if len(size) < 2:
        return None                                        # one child: skipped
    if sum(1 for k in size if size[k] == 1) * 2 > len(size):
        return None                                        # mostly singletons: skipped
    out = collections.defaultdict(list)
    for v in members:
        out[key[v]].append(v)
    return {k: sorted(v) for k, v in out.items()}


def plan_loc(m, members, levels_below):
    """R11-R13 for one group: a plan [(kind, key, nodes, sub-plan)].

    Localisation groups where the split is admissible. A group above LEAF_MAX is cut into connectivity parts:
    one level down when a level is left (R12), otherwise as parts of the same compartment (R13).
    """
    loc = loc_split(m, members)
    if loc is None:
        if len(members) <= LEAF_MAX:
            return []
        return [("part", None, part, []) for part in m.parts(members, LEAF_MAX)]
    plan = []
    for key in sorted(loc):
        nodes = loc[key]
        if len(nodes) <= LEAF_MAX:
            plan.append(("loc", key, nodes, []))
        elif levels_below > 0:
            plan.append(("loc", key, nodes, [("part", None, part, []) for part in m.parts(nodes, LEAF_MAX)]))
        else:
            plan.extend(("locpart", key, part, []) for part in m.parts(nodes, LEAF_MAX))
    return plan


def plan_type_first(m, members):
    """Reaction type (R10), then localisation (R11). Without an admissible type split: localisation alone."""
    types = type_split(m, members)
    if types is None:
        return plan_loc(m, members, 1)
    return [("sbo", term, types[term], plan_loc(m, types[term], 0)) for term in sorted(types)]


def plan_loc_first(m, members):
    """Localisation (R11), then reaction type (R10) inside a localisation group above LEAF_MAX."""
    loc = loc_split(m, members)
    if loc is None:
        return plan_type_first(m, members)
    plan = []
    for key in sorted(loc):
        nodes = loc[key]
        sub = []
        if len(nodes) > LEAF_MAX:
            types = type_split(m, nodes)
            if types is not None and max(len(v) for v in types.values()) <= LEAF_MAX:
                sub = [("sbo", term, types[term], []) for term in sorted(types)]
            else:
                sub = [("part", None, part, []) for part in m.parts(nodes, LEAF_MAX)]
        plan.append(("loc", key, nodes, sub))
    return plan


def tiny_leaves(plan):
    """Number of leaf groups of a plan with fewer than GROUP_MIN nodes."""
    return sum(tiny_leaves(sub) if sub else int(len(nodes) < GROUP_MIN) for _kind, _key, nodes, sub in plan)


def representative(m, members):
    """R14: the member a connectivity group is named after: the species with the highest k * (k / d) ** 2,
    k = its reactions inside the group, d = all its reactions (a currency metabolite reacts mostly elsewhere
    and scores low). Ties: the smaller id. A group without species is named after its first reaction."""
    inside = set(members)
    pool = sorted(v for v in members if m.nodes[v]["kind"] == "species") or sorted(members)

    def score(v):
        k = len(m.adj[v] & inside)
        return Fraction(k ** 3, len(m.adj[v]) ** 2) if m.adj[v] else Fraction(0)

    return min(pool, key=lambda v: (-score(v), v))


def node_label(m, v):
    return m.nodes[v]["name"].strip() or v


# ==================================================================== the whole grouping
def run(features_dir, out_dir):
    m = Model(features_dir)
    k_target = min(K_MAX, max(K_MIN, math.ceil(m.n / LEAF_MAX)))
    cap = max(LEAF_MAX, math.floor(m.n * SHARE_CAP))
    l1, rule, naming, info = function_split(m, m.ids, None, k_target, cap)
    meta, rows = {}, {v: {"node_id": v, "L1": "", "L2": "", "L3": "", "rule": rule[v], "clone_in": ""} for v in m.ids}
    info["sub"], info["order"] = {}, {}

    def add_meta(gid, level, parent, label, basis, source, xref, comp=""):
        meta[gid] = dict(group_id=gid, level=level, parent=parent, label=label, basis=basis, name_source=source,
                         name_xref=xref, compartment_id=comp)

    def topo_meta(gid, level, parent, members):            # R14
        rep = representative(m, members)
        add_meta(gid, level, parent, node_label(m, rep), "topology", "member", rep)

    def term_meta(gid, level, parent, t, below=None):      # R8
        names = m.term_name.get(t)
        if names and len(names) == len(t.split(".")):      # the name path as written in the notes;
            skip = below.count(".") + 1 if below else 0    # a child drops the part its parent already shows
            label = ": ".join(names[skip:])
        else:
            label = f"{FUNCTION_NOTE_KEYS[0]} {t}"
        add_meta(gid, level, parent, label, "function", "notes", f"{FUNCTION_NOTE_KEYS[0]}:{t}")

    def sbo_meta(gid, level, parent, term, members):
        name = ""
        for v in members:
            p = m.nodes[v]["sbo_path"].split(">")
            if term in p:
                name = m.nodes[v]["sbo_path_names"].split(" > ")[p.index(term)]
                break
        add_meta(gid, level, parent, name or term, "reaction_type", "sbo", term)

    def loc_meta(gid, level, parent, key, suffix=""):
        cid, _, edge = key.partition("@")
        name = m.comp.get(cid, {}).get("name") or cid
        add_meta(gid, level, parent, name + (" (boundary)" if edge else "") + suffix, "localisation", "compartment",
                 cid, key)

    def place(members, level, gid):
        for v in members:
            rows[v][f"L{level}"] = gid

    def emit(plan, level, parent):
        """Write a plan [(kind, key, nodes, sub-plan)] below `parent`."""
        counter = collections.Counter()
        for kind, key, nodes, sub in plan:
            if kind == "loc":
                gid = f"{parent}.{key.replace('@boundary', '_b')}"
                loc_meta(gid, level, parent, key)
            elif kind == "locpart":                        # R13: named by compartment and one member
                counter[key] += 1
                gid = f"{parent}.{key.replace('@boundary', '_b')}.p{counter[key]}"
                loc_meta(gid, level, parent, key, f" [{node_label(m, representative(m, nodes))}]")
            elif kind == "sbo":
                gid = f"{parent}.{key.replace('SBO:', 'sbo')}"
                sbo_meta(gid, level, parent, key, nodes)
            else:
                counter["p"] += 1
                gid = f"{parent}.p{counter['p']}"
                topo_meta(gid, level, parent, nodes)
            place(nodes, level, gid)
            emit(sub, level + 1, gid)

    order = sorted(set(l1.values()), key=lambda g: (-sum(1 for v in l1 if l1[v] == g), g))
    for i, g in enumerate(order, start=1):
        gid = f"M{i:02d}"
        members = sorted(v for v in m.ids if l1[v] == g)
        if naming[g]:
            term_meta(gid, 1, "", naming[g]["term"])
        else:
            topo_meta(gid, 1, "", members)
        info.setdefault("l1_naming", {})[gid] = naming[g]
        place(members, 1, gid)
        if len(members) <= LEAF_MAX:                       # within the leaf bound: localisation only (R11)
            emit(plan_loc(m, members, 1), 2, gid)
            continue
        if naming[g]:                                      # R9: finer function terms, when they discriminate
            k2 = min(K_MAX, max(2, math.ceil(len(members) / LEAF_MAX)))
            cap2 = max(LEAF_MAX, math.floor(len(members) * SHARE_CAP))
            l2, _rule2, naming2, info2 = function_split(m, members, naming[g]["term"], k2, cap2)
            if sum(1 for k in naming2 if naming2[k]) >= 2:
                info["sub"][gid] = info2
                keys = sorted(set(l2.values()), key=lambda k: (-sum(1 for v in l2 if l2[v] == k), k))
                for j, k in enumerate(keys, start=1):
                    nodes = sorted(v for v in members if l2[v] == k)
                    sid = f"{gid}.f{j}"
                    if naming2[k]:
                        term_meta(sid, 2, gid, naming2[k]["term"], naming[g]["term"])
                    else:
                        topo_meta(sid, 2, gid, nodes)
                    info.setdefault("l2_naming", {})[sid] = naming2[k]
                    place(nodes, 2, sid)
                    emit(plan_loc(m, nodes, 0), 3, sid)
                continue
        # R10 and R11 in the order that leaves fewer groups below GROUP_MIN (ties: type first)
        plan_a, plan_b = plan_type_first(m, members), plan_loc_first(m, members)
        first = "type" if tiny_leaves(plan_a) <= tiny_leaves(plan_b) else "localisation"
        info["order"][gid] = (first, tiny_leaves(plan_a), tiny_leaves(plan_b))
        emit(plan_a if first == "type" else plan_b, 2, gid)

    # R15: clones of hub species
    def chain(v):
        return [rows[v][lv] for lv in ("L3", "L2", "L1") if rows[v][lv]]

    for s in sorted(m.hubs):
        own, targets = set(chain(s)), set()
        for r in sorted(m.adj[s]):
            for gid in chain(r):                           # deepest group of the reaction that can hold s
                c = meta[gid]["compartment_id"]
                if c and c != loc_key(m, s):
                    continue
                if not c:                                  # stepped over a foreign compartment: use the sibling
                    kids = sorted(k for k in meta if meta[k]["parent"] == gid          # group of s's own place
                                  and meta[k]["compartment_id"] == loc_key(m, s))
                    gid = kids[0] if kids else gid
                targets.add(gid)
                break
        rows[s]["clone_in"] = "|".join(sorted(t for t in targets if t not in own))
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", ["node_id", "L1", "L2", "L3", "rule", "clone_in"], [rows[v] for v in m.ids])
    write_tsv(out / "groups_meta.tsv", ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref",
                                        "compartment_id"], [meta[g] for g in sorted(meta)])
    return m, rows, meta, info


def main(argv=None):
    ap = argparse.ArgumentParser(description="hybrid grouping strategy (see RULES.md)")
    ap.add_argument("--features", required=True, help="folder with nodes.tsv, edges.tsv, compartments.tsv")
    ap.add_argument("--sbml", required=False, help="the SBML model (part of the common command line; not opened)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    a = ap.parse_args(argv)
    m, rows, meta, info = run(a.features, a.out)
    n1 = sum(1 for g in meta.values() if g["level"] == 1)
    print(f"{m.n} nodes -> {len(meta)} groups, {n1} at level 1 (aim {info['k_target']}, level-1 cap {info['cap']} "
          f"nodes) -> {a.out}")


if __name__ == "__main__":
    main()
