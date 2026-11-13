"""Readability metrics for a nested grouping of an SBML model (BioHackathon 2026, project 22).

Run:
    python metrics_readability.py panel/features/MODEL strategies/NAME/out/MODEL
    python metrics_readability.py panel/features/MODEL strategies/NAME/out/MODEL --reference REF.tsv --json
    python -m pytest test_metrics_readability.py -q

The first two lines print every metric, the verdict of each rule (pass / warn / fail / n/a) and the
combined readability score for one grouping. REF.tsv is the optional curated table (node_id, label).

Use from Python:
    import pandas as pd, metrics_readability as mr
    rd = lambda p: pd.read_csv(p, sep="\\t", dtype=str, keep_default_na=False)
    m = mr.metrics(rd("nodes.tsv"), rd("edges.tsv"), rd("compartments.tsv"), rd("groups.tsv"), rd("groups_meta.tsv"))
    mr.verdicts(m); mr.readability_score(m)

Viewpoint: what a reader can take in from a map drawn with one box per group. Rules and thresholds are
documented in rules_readability.md; the constants below are the single place where thresholds live.

Vocabulary used in the code
    path(v)   the non-empty group ids L1..L3 of node v; home(v) = the last of them (its deepest group)
    size(g)   number of nodes with g on their path;  load(g) = number of nodes whose home is g
    box       a group as drawn. A space group `c` and its sibling `c@boundary` (basis localisation, same
              parent) are ONE box when boxes are counted: a membrane is the border of the same compartment.
    NaN       the rule does not apply to this MODEL (e.g. one compartment, no annotation, 20 nodes or
              fewer) or the grouping has no nesting / no clone to measure. When the grouping itself leaves
              a denominator empty in a way that hides a defect (no group is named), the metric takes its
              worst value instead of NaN, so that omission is never rewarded.

The grouping is expected to pass `common/contract.py validate` first. Nodes of nodes.tsv missing from
groups.tsv are ignored and reported in `placed_share`.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import pathlib
import sys

import pandas as pd

NAN = float("nan")

# ----------------------------------------------------------------------------- thresholds (see .md)
SMALL_MODEL = 20          # nodes; a model this small reads as one box (Ghoniem et al. 2004)
BOX_MAX = 50              # nodes in one box before it stops being readable (Yoghourdjian et al. 2020)
BOX_HARD = 100            # reported only (box_over100_node_share)
LONG_SPAN = 100           # an edge is long when no box of at most this many nodes holds both its ends
FAN_MIN, FAN_MAX = 2, 12  # child boxes per parent box
TOP_MIN, TOP_MAX = 3, 12  # boxes at the top level
COMFORT_FAN, COMFORT_BOX = 7, 25   # used only to derive the depth a model of N nodes needs
RECALL_MIN = 0.5          # a name is informative if the box holds at least half of what it names
BACKED_SOURCES = {"sbo", "ec", "go", "xref", "compartment"}

# rule id -> (metric, direction, pass, warn, zero, weight). direction "low": smaller is better.
# score = 1 at/inside `pass`, 0.5 at `warn`, 0 at/after `zero`, linear in between.
RULES = collections.OrderedDict([
    ("RD01_TOP",       ("top_groups",           "range", (TOP_MIN, TOP_MAX), (2, 20), (1, 40), 3.0)),
    ("RD02_BOXSIZE",   ("oversize_node_share",  "low",  0.10, 0.30, 0.70, 3.0)),
    ("RD03_SINGLETON", ("singleton_share",      "low",  0.10, 0.25, 0.50, 3.0)),
    ("RD04_NAMED",     ("named_share",          "high", 0.90, 0.60, 0.20, 3.0)),
    ("RD05_LONGEDGE",  ("long_edge_share",      "low",  0.35, 0.60, 0.90, 2.0)),
    ("RD06_CLONE",     ("cloned_species_share", "low",  0.10, 0.25, 0.50, 2.0)),
    ("RD07_DISTINCT",  ("distinct_share",       "high", 0.95, 0.80, 0.50, 2.0)),
    ("RD08_FAN",       ("fan_ok_share",         "high", 0.90, 0.70, 0.40, 2.0)),
    ("RD09_BACKED",    ("backed_share",         "high", 0.80, 0.50, 0.20, 1.5)),
    ("RD10_COMPDRAW",  ("comp_redraw",          "low",  1.5,  3.0,  6.0,  1.5)),
    ("RD11_BALANCE",   ("size_balance",         "high", 0.50, 0.25, 0.10, 1.0)),
    ("RD12_DEPTH",     ("depth_excess",         "low",  0.0,  1.0,  2.0,  1.0)),
])
CORE_RULES = ("RD01_TOP", "RD02_BOXSIZE", "RD03_SINGLETON")   # gate of the combined score


# ----------------------------------------------------------------------------- helpers
def _norm(label: str) -> str:
    return " ".join(str(label).split()).casefold()


def _col(df: pd.DataFrame, name: str) -> list:
    return df[name].astype(str).tolist() if name in df.columns else [""] * len(df)


def _eff_number(counts) -> float:
    """exp(Shannon entropy) of a list of positive counts: the effective number of parts."""
    tot = float(sum(counts))
    if tot <= 0:
        return NAN
    h = -sum((c / tot) * math.log(c / tot) for c in counts if c > 0)
    return math.exp(h)


def depth_needed(n_nodes: int) -> int:
    """Smallest number of levels d >= 1 with COMFORT_FAN**d * COMFORT_BOX >= n_nodes."""
    d = 1
    while COMFORT_FAN ** d * COMFORT_BOX < n_nodes:
        d += 1
    return d


def _ec_class(x: str) -> str:
    x = x.split(":", 1)[1] if ":" in x else x
    parts = [p for p in x.strip().strip(".").split(".") if p not in ("", "-")]
    return ".".join(parts)


class _Model:
    """Indexes built once from the contract tables."""

    def __init__(self, nodes, edges, comps, groups, meta):
        self.kind = dict(zip(_col(nodes, "node_id"), _col(nodes, "kind")))
        ids = _col(nodes, "node_id")
        loc = _col(nodes, "loc_space")
        comp = _col(nodes, "compartment")
        self.loc = {i: (l or c) for i, l, c in zip(ids, loc, comp)}
        self.loc_kind = dict(zip(ids, _col(nodes, "loc_kind")))
        self.sbo_path = {i: set(filter(None, s.split(">"))) for i, s in zip(ids, _col(nodes, "sbo_path"))}
        self.ec = {i: [e for e in s.split("|") if e] for i, s in zip(ids, _col(nodes, "ec"))}
        self.go = {i: set(filter(None, s.split("|"))) for i, s in zip(ids, _col(nodes, "go"))}
        self.xrefs = {i: [x.split("=", 1)[-1].casefold() for x in s.split("|") if x]
                      for i, s in zip(ids, _col(nodes, "xrefs"))}
        self.annotated = any(self.sbo_path[i] or self.ec[i] or self.go[i] or self.xrefs[i] for i in ids)
        self.comp_ids = set(_col(comps, "compartment_id"))
        self.n_all = len(ids)

        self.meta = {}
        for row in meta.to_dict("records"):
            self.meta[str(row.get("group_id", ""))] = {k: str(v) for k, v in row.items()}

        self.path = {}
        self.clones = {}
        for row in groups.to_dict("records"):
            nid = str(row.get("node_id", ""))
            if nid not in self.kind or nid in self.path:
                continue
            p = []
            for lv in ("L1", "L2", "L3"):
                g = str(row.get(lv, ""))
                if not g:
                    break
                p.append(g)
            if not p:
                continue
            self.path[nid] = tuple(p)
            self.clones[nid] = [c for c in str(row.get("clone_in", "")).split("|") if c]
        self.nodes = sorted(self.path)
        self.n = len(self.nodes)

        self.parent, self.members, self.load = {}, collections.defaultdict(list), collections.Counter()
        self.children = collections.defaultdict(set)
        for nid in self.nodes:
            p = self.path[nid]
            self.load[p[-1]] += 1
            up = ""
            for g in p:
                self.members[g].append(nid)
                self.parent.setdefault(g, up)
                self.children[up].add(g)
                up = g
        self.groups = sorted(self.members)
        self.size = {g: len(self.members[g]) for g in self.groups}
        self.gpath = {}
        for g in self.groups:
            chain = [g]
            while self.parent.get(chain[-1], ""):
                chain.append(self.parent[chain[-1]])
            self.gpath[g] = tuple(reversed(chain))

        self.edges = [(r, s) for r, s in zip(_col(edges, "reaction_id"), _col(edges, "species_id"))
                      if r in self.path and s in self.path]

    # -- group attributes ---------------------------------------------------------------------------
    def m(self, g, key):
        return self.meta.get(g, {}).get(key, "")

    def is_loc(self, g):
        return self.m(g, "basis") == "localisation"

    def label(self, g):
        return _norm(self.m(g, "label"))

    def box_key(self, g):
        """Space group and its sibling boundary group share one key (they are one drawn compartment)."""
        cid = self.m(g, "compartment_id")
        if self.is_loc(g) and cid:
            return ("loc", cid.partition("@")[0])
        return ("grp", g)

    def n_boxes(self, parent):
        kids = self.children.get(parent, ())
        return len({self.box_key(k) for k in kids}) + (1 if parent and kids and self.load.get(parent, 0) else 0)

    def group_path(self, g):
        """Path of any group id (also one that only appears in clone_in), or None if unknown."""
        if g in self.gpath:
            return self.gpath[g]
        if g in self.meta:
            chain, seen = [g], {g}
            while self.m(chain[-1], "parent") and self.m(chain[-1], "parent") not in seen:
                chain.append(self.m(chain[-1], "parent"))
                seen.add(chain[-1])
            return tuple(reversed(chain))
        return None

    # -- which nodes does an identifier apply to ----------------------------------------------------
    def carriers(self, g):
        """Nodes of the model to which the identifier behind g's name applies, or None if it cannot
        be checked from the tables (name_source notes / none, or an unparsable name_xref)."""
        src, x = self.m(g, "name_source"), self.m(g, "name_xref").strip()
        if not x or src in ("none", "", "notes"):
            return None
        if src == "member":
            return {x} if x in self.kind else set()
        if src == "compartment" or (self.is_loc(g) and x.partition("@")[0] in self.comp_ids):
            cid, _, edge = (self.m(g, "compartment_id") or x).partition("@")
            return {i for i in self.kind if self.loc[i] == cid and (not edge or self.loc_kind[i] == "boundary")}
        if x.upper().startswith("SBO:"):
            return {i for i in self.kind if x.upper() in self.sbo_path[i]}
        if x.upper().startswith("GO:"):
            xl = x.casefold()
            return {i for i in self.kind if x.upper() in self.go[i] or any(r.endswith(xl) for r in self.xrefs[i])}
        if src == "ec" or x.casefold().startswith("ec-code:") or x.casefold().startswith("ec:"):
            cls = _ec_class(x)
            if not cls:
                return set()
            return {i for i in self.kind if any(e == cls or e.startswith(cls + ".") for e in self.ec[i])}
        xl = x.casefold()
        return {i for i in self.kind if any(r == xl or r.startswith(xl + ".") for r in self.xrefs[i])}


# ----------------------------------------------------------------------------- the metrics
def metrics(nodes, edges, comps, groups, meta, reference=None) -> dict:
    """All readability metrics of one grouping. Inputs are the contract tables as pandas DataFrames
    read with dtype=str, keep_default_na=False. `reference` is the optional curated table
    (node_id, label; several rows per node allowed). Returns {metric name: float}; NaN = not applicable.

    Rule metrics (thresholds in RULES): top_groups, oversize_node_share, singleton_share, named_share,
    long_edge_share, cloned_species_share, distinct_share, fan_ok_share, backed_share, comp_redraw,
    size_balance, depth_excess. Every other key is a descriptive or counter-rule value quoted in rules_readability.md.
    """
    M = _Model(nodes, edges, comps, groups, meta)
    out = collections.OrderedDict()
    N, G = M.n, len(M.groups)
    out["n_nodes"] = float(N)
    out["n_groups"] = float(G)
    out["placed_share"] = N / M.n_all if M.n_all else NAN
    if N == 0 or G == 0:
        for _rid, spec in RULES.items():
            out[spec[0]] = NAN
        return dict(out)
    for lv in (1, 2, 3):
        out[f"n_groups_L{lv}"] = float(sum(1 for g in M.groups if len(M.gpath[g]) == lv))

    # RD01 top-level boxes -----------------------------------------------------------------------
    top = M.n_boxes("")
    out["top_boxes"] = float(top)
    out["top_all_localisation"] = float(all(M.is_loc(g) for g in M.children[""]))
    out["top_groups"] = float(top) if N > SMALL_MODEL else NAN

    # RD02 box size ------------------------------------------------------------------------------
    loads = [M.load[g] for g in M.groups if M.load[g] > 0]
    out["oversize_node_share"] = sum(c for c in loads if c > BOX_MAX) / N
    out["box_over100_node_share"] = sum(c for c in loads if c > BOX_HARD) / N
    out["box_size_max"] = float(max(loads))
    out["box_size_median"] = float(pd.Series(loads).median())

    # RD03 singletons ----------------------------------------------------------------------------
    cid_count = collections.Counter(M.m(g, "compartment_id") for g in M.groups if M.is_loc(g))
    def forced_single(g):   # the only box of a compartment that holds one node: biology, not the strategy
        return M.is_loc(g) and M.m(g, "compartment_id") and cid_count[M.m(g, "compartment_id")] == 1
    single = [g for g in M.groups if M.size[g] == 1 and not forced_single(g)]
    out["singleton_share"] = len(single) / G
    out["singleton_node_share"] = sum(1 for v in M.nodes if M.size[M.path[v][-1]] == 1) / N

    # RD04 names ---------------------------------------------------------------------------------
    named = [g for g in M.groups if M.label(g)]
    out["named_share"] = len(named) / G
    out["named_node_share"] = sum(1 for v in M.nodes if M.label(M.path[v][-1])) / N

    # RD09 identifier-backed names ---------------
    carriers = {g: M.carriers(g) for g in named}
    def backed(g):
        if M.m(g, "name_source") not in BACKED_SOURCES or not M.m(g, "name_xref").strip():
            return False
        c = carriers[g]
        return bool(c) and any(v in c for v in M.members[g])
    multi_comp = len({M.loc[v] for v in M.nodes}) > 1
    if not M.annotated and not multi_comp:
        out["backed_share"] = NAN            # the model carries no identifier a name could come from
    else:
        out["backed_share"] = (sum(backed(g) for g in named) / len(named)) if named else 0.0
    out["member_named_share"] = (sum(M.m(g, "name_source") == "member" for g in named) / len(named)) if named else NAN

    # RD07 distinct and informative names --------------------------------------------------------
    sib = collections.Counter((M.parent[g], M.label(g)) for g in named)
    # a space box and its own boundary box may share the compartment's name: one drawn compartment
    sib_boxes = collections.defaultdict(set)
    for g in named:
        sib_boxes[(M.parent[g], M.label(g))].add(M.box_key(g))
    glob = collections.Counter(M.label(g) for g in named)
    dup, same_anc, low_recall, ok = 0, 0, 0, 0
    recalls = []
    for g in named:
        key = (M.parent[g], M.label(g))
        is_dup = sib[key] > 1 and len(sib_boxes[key]) > 1
        anc = any(M.label(a) == M.label(g) for a in M.gpath[g][:-1])
        c = carriers[g]
        rec = NAN
        if c is not None and M.m(g, "name_source") != "member":
            scope = set(M.members[M.parent[g]]) if M.parent[g] else set(M.nodes)
            in_scope = c & scope
            rec = (len(c & set(M.members[g])) / len(in_scope)) if in_scope else 0.0
            recalls.append(rec)
        low = (rec == rec) and rec < RECALL_MIN
        dup += is_dup
        same_anc += anc
        low_recall += low
        ok += not (is_dup or anc or low)
    out["distinct_share"] = ok / len(named) if named else 0.0
    out["sibling_dup_share"] = dup / len(named) if named else NAN
    out["same_as_ancestor_share"] = same_anc / len(named) if named else NAN
    out["name_low_recall_share"] = low_recall / len(named) if named else NAN
    out["name_recall_mean"] = (sum(recalls) / len(recalls)) if recalls else NAN
    out["global_dup_share"] = (sum(1 for g in named if glob[M.label(g)] > 1) / len(named)) if named else NAN
    out["label_len_median"] = float(pd.Series([len(M.m(g, "label")) for g in named]).median()) if named else NAN

    # RD08 children per parent -------------------------------------------------------------------
    parents = [g for g in M.groups if M.children.get(g)]
    def fan_ok(g):
        k = M.n_boxes(g)
        if k > FAN_MAX:
            return False
        return k >= FAN_MIN or all(M.is_loc(c) for c in M.children[g])
    out["fan_ok_share"] = (sum(fan_ok(g) for g in parents) / len(parents)) if parents else NAN
    out["fan_max"] = float(max(M.n_boxes(g) for g in parents)) if parents else NAN
    out["only_child_share"] = (sum(len(M.children[g]) == 1 and M.load.get(g, 0) == 0 for g in parents)
                               / len(parents)) if parents else NAN

    # RD05 long edges, RD06 clones ---------------------------------------------------------------
    n_species = sum(1 for v in M.nodes if M.kind[v] == "species")
    copies, n_clone, useless = {}, 0, 0
    adj = collections.defaultdict(set)
    for r, s in M.edges:
        adj[s].add(r)
    for v in M.nodes:
        cps = [M.path[v]]
        for c in M.clones[v]:
            p = M.group_path(c)
            if p is None:
                continue
            n_clone += 1
            cps.append(p)
            if not any(M.path[r][:len(p)] == p for r in adj[v]):
                useless += 1
        copies[v] = cps
    def crossings(pr, pc):
        k = 0
        while k < min(len(pr), len(pc)) and pr[k] == pc[k]:
            k += 1
        return len(pr) + len(pc) - 2 * k

    def is_long(pr, pc):
        """No common box of readable span: the deepest group holding both ends has more than
        LONG_SPAN nodes, or the two ends share no group at all."""
        k = 0
        while k < min(len(pr), len(pc)) and pr[k] == pc[k]:
            k += 1
        return k == 0 or M.size.get(pr[k - 1], 0) > LONG_SPAN

    if M.edges:
        raw = [is_long(M.path[r], M.path[s]) for r, s in M.edges]
        long_cl = sum(all(is_long(M.path[r], p) for p in copies[s]) for r, s in M.edges)
        d = [min(crossings(M.path[r], p) for p in copies[s]) for r, s in M.edges]
        out["long_edge_share"] = long_cl / len(M.edges)
        out["long_edge_share_raw"] = sum(raw) / len(M.edges)
        out["l1_cross_share"] = sum(all(M.path[r][0] != p[0] for p in copies[s]) for r, s in M.edges) / len(M.edges)
        out["cut_edge_share"] = sum(x > 0 for x in d) / len(M.edges)
        out["boundaries_crossed_mean"] = sum(d) / len(d)
        # Newman-Girvan modularity of the L1 partition on the species-reaction graph (clones ignored)
        deg, inside = collections.Counter(), collections.Counter()
        for r, s in M.edges:
            deg[M.path[r][0]] += 1
            deg[M.path[s][0]] += 1
            if M.path[r][0] == M.path[s][0]:
                inside[M.path[r][0]] += 1
        m2 = 2.0 * len(M.edges)
        out["modularity_L1"] = sum(inside[g] / len(M.edges) - (deg[g] / m2) ** 2 for g in deg)
    else:
        raw = []
        for k in ("long_edge_share", "long_edge_share_raw", "l1_cross_share", "cut_edge_share",
                  "boundaries_crossed_mean", "modularity_L1"):
            out[k] = NAN
    if n_species:
        # clones that would shorten every long edge a clone can shorten: one copy of the species in the
        # home box of each reaction it is long-connected to
        demand = {(sp, M.path[r][-1]) for (r, sp), lg in zip(M.edges, raw)
                  if lg and M.path[r][-1] not in M.path[sp]}
        out["clone_load"] = n_clone / n_species
        out["clone_demand"] = len(demand) / n_species
        out["cloned_species_share"] = sum(1 for v in M.nodes if len(copies[v]) > 1) / n_species
    else:
        out["clone_load"] = out["clone_demand"] = out["cloned_species_share"] = NAN
    out["clone_useless_share"] = (useless / n_clone) if n_clone else NAN

    # RD10 how often a compartment has to be drawn -----------------------------------------------
    spaces = collections.defaultdict(set)
    for g in M.groups:
        spaces[g] = {M.loc[v] for v in M.members[g]}
    if multi_comp:
        pieces = collections.defaultdict(set)
        for v in M.nodes:
            box = ""                       # the map itself
            for g in reversed(M.path[v]):
                if len(spaces[g]) > 1:
                    box = g
                    break
            pieces[M.loc[v]].add(box)
        cnt = [len(b) for b in pieces.values()]
        out["comp_redraw"] = sum(cnt) / len(cnt)
        out["comp_redraw_max"] = float(max(cnt))
    else:
        out["comp_redraw"] = out["comp_redraw_max"] = NAN
    out["mixed_box_node_share"] = sum(1 for v in M.nodes if len(spaces[M.path[v][-1]]) > 1) / N if multi_comp else NAN
    # the literal count: localisation boxes declared for the same compartment under different parents
    loc_parents = collections.defaultdict(set)
    for g in M.groups:
        if M.is_loc(g) and M.m(g, "compartment_id"):
            loc_parents[M.m(g, "compartment_id").partition("@")[0]].add(M.parent[g])
    out["loc_box_repeat"] = (sum(len(x) for x in loc_parents.values()) / len(loc_parents)) if loc_parents else NAN

    # RD11 balance -------------------------------------------------------------------------------
    out["size_balance"] = (_eff_number(loads) / len(loads)) if len(loads) > 1 else NAN
    out["largest_box_share"] = max(loads) / N

    # RD12 depth ---------------------------------------------------------------------------------
    own = [sum(1 for g in M.path[v] if not M.is_loc(g)) for v in M.nodes]
    out["depth_mean"] = sum(own) / N
    out["depth_levels_mean"] = sum(len(M.path[v]) for v in M.nodes) / N
    out["depth_max"] = float(max(len(M.path[v]) for v in M.nodes))
    out["depth_needed"] = float(depth_needed(N))
    out["depth_excess"] = max(0.0, out["depth_mean"] - out["depth_needed"])

    # optional curated reference: can a reader find a curated pathway in one place? ----------------
    out["ref_scatter"] = out["ref_mix"] = out["ref_labelled_share"] = NAN
    if reference is not None and len(reference):
        lab = collections.defaultdict(set)
        for nid, lb in zip(_col(reference, "node_id"), _col(reference, "label")):
            if nid in M.path and lb:
                lab[nid].add(lb)
        if lab:
            out["ref_labelled_share"] = len(lab) / N
            by_label, by_box = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
            for nid, lbs in lab.items():
                for lb in lbs:
                    by_label[lb][M.path[nid][0]] += 1
                    by_box[M.path[nid][-1]][lb] += 1
            big = {lb: c for lb, c in by_label.items() if sum(c.values()) >= 3}
            if big:
                w = {lb: sum(c.values()) for lb, c in big.items()}
                out["ref_scatter"] = sum(w[lb] * _eff_number(list(c.values())) for lb, c in big.items()) / sum(w.values())
            wb = {b: sum(c.values()) for b, c in by_box.items()}
            out["ref_mix"] = sum(wb[b] * _eff_number(list(c.values())) for b, c in by_box.items()) / sum(wb.values())
    return dict(out)


# ----------------------------------------------------------------------------- verdicts and score
def _score_one(x, direction, p, w, z):
    if x != x:
        return NAN
    if direction == "range":
        (plo, phi), (wlo, whi), (zlo, zhi) = p, w, z
        if plo <= x <= phi:
            return 1.0
        if x < plo:
            return _score_one(x, "high", plo, wlo, zlo)
        return _score_one(x, "low", phi, whi, zhi)
    if direction == "high":
        x, p, w, z = -x, -p, -w, -z
    if x <= p:
        return 1.0
    if x <= w:
        return 1.0 - 0.5 * (x - p) / (w - p)
    if x < z:
        return 0.5 - 0.5 * (x - w) / (z - w)
    return 0.0


def _waive_top_lower(m, rid, x):
    """Top level made only of compartments: their number is given by the model, so no lower bound."""
    return rid == "RD01_TOP" and m.get("top_all_localisation", 0.0) == 1.0 and x == x and x < TOP_MIN


def rule_scores(m: dict) -> dict:
    """Rule id -> score in [0, 1] (1 = inside the pass zone, 0.5 = at the warn limit), NaN = n/a."""
    out = collections.OrderedDict()
    for rid, (metric, direction, p, w, z, _wt) in RULES.items():
        x = m.get(metric, NAN)
        out[rid] = 1.0 if _waive_top_lower(m, rid, x) else _score_one(x, direction, p, w, z)
    return out


def verdicts(m: dict) -> dict:
    """Rule id -> 'pass' | 'warn' | 'fail' | 'na'."""
    out = collections.OrderedDict()
    for rid, s in rule_scores(m).items():
        out[rid] = "na" if s != s else "pass" if s >= 1.0 else "warn" if s >= 0.5 else "fail"
    return out


def readability_score(m: dict) -> float:
    """One number in [0, 1]: weighted mean of the applicable rule scores, multiplied by a gate
    0.5 + 0.5 * min(score of RD01, RD02, RD03) so that a map that fails a core rule cannot rank first."""
    s = rule_scores(m)
    num = sum(RULES[r][5] * v for r, v in s.items() if v == v)
    den = sum(RULES[r][5] for r, v in s.items() if v == v)
    if den == 0:
        return NAN
    core = [s[r] for r in CORE_RULES if s[r] == s[r]]
    gate = 0.5 + 0.5 * min(core) if core else 1.0
    return (num / den) * gate


# ----------------------------------------------------------------------------- command line
def read_tsv(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def metrics_from_dirs(features_dir, out_dir, reference=None) -> dict:
    f, o = pathlib.Path(features_dir), pathlib.Path(out_dir)
    ref = read_tsv(reference) if reference else None
    return metrics(read_tsv(f / "nodes.tsv"), read_tsv(f / "edges.tsv"), read_tsv(f / "compartments.tsv"),
                   read_tsv(o / "groups.tsv"), read_tsv(o / "groups_meta.tsv"), ref)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Readability metrics of one grouping (contract format).")
    ap.add_argument("features_dir", help="panel/features/MODEL")
    ap.add_argument("out_dir", help="strategies/NAME/out/MODEL (groups.tsv, groups_meta.tsv)")
    ap.add_argument("--reference", help="optional curated table (node_id, label), tab-separated")
    ap.add_argument("--json", action="store_true", help="print one JSON object instead of the table")
    a = ap.parse_args(argv)
    m = metrics_from_dirs(a.features_dir, a.out_dir, a.reference)
    v, s = verdicts(m), rule_scores(m)
    if a.json:
        clean = lambda d: {k: (None if isinstance(x, float) and x != x else x) for k, x in d.items()}
        print(json.dumps({"metrics": clean(m), "verdicts": v, "rule_scores": clean(s),
                          "readability_score": readability_score(m)}, indent=1))
        return
    for rid, (metric, *_rest) in RULES.items():
        x = m.get(metric, NAN)
        print(f"{rid:<15}{metric:<22}{x:>9.3f}  {v[rid]}")
    print(f"{'readability_score':<37}{readability_score(m):>9.3f}")
    rule_metrics = {spec[0] for spec in RULES.values()}
    for k, x in m.items():
        if k not in rule_metrics:
            print(f"  {k:<35}{x:>9.3f}")


if __name__ == "__main__":
    main()
