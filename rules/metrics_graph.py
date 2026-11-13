"""Graph-theoretic metrics for a nested grouping of an SBML model (BioHackathon 2026, project 22).

Run (PowerShell or bash, from the folder that holds this file):
    python metrics_graph.py score --features panel/features/skm_ja --grouping strategies/NAME/out/skm_ja
    python metrics_graph.py score --features panel/features/skm_ja --grouping strategies/NAME/out/skm_ja --reference reference/skm_ja.tsv --json skm_ja_graph.json
    python metrics_graph.py demo  --features panel/features/skm_ja
    python -m pytest test_metrics_graph.py -q

`score` prints the headline metric and the pass / warn / fail verdict of rules G1..G12 for one grouping
in the contract format (groups.tsv, groups_meta.tsv); --json writes every metric. `demo` builds the
synthetic groupings (one group, all singletons, one group per compartment, random with realistic sizes,
modularity optimiser, and adversarial cases) for a model and prints their verdicts side by side.
Needs pandas and numpy only. Deterministic: no random numbers in `metrics`.

Python:
    import pandas as pd, metrics_graph as mg
    rd = lambda p: pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    m = mg.metrics(rd("nodes.tsv"), rd("edges.tsv"), rd("compartments.tsv"), rd("groups.tsv"), rd("groups_meta.tsv"))
    verdicts = mg.judge(m)          # {"G1": "pass" | "warn" | "fail" | "na", ...}
    total = mg.score(m)             # 0..1, the combination of rules_graph.md

The graph
    Bipartite: species on one side, reactions on the other, one edge per distinct (reaction, species)
    pair of edges.tsv (reactant, product and modifier alike; a species with two roles in one reaction
    is one edge). d_s = number of reactions of species s.

Hub treatment (data-derived, no list of molecules)
    1. Species-normalised weight: edge (r, s) weighs w = 1 / d_s, so every species has a total weight
       of 1 ("one vote") whether it takes part in 1 reaction or in 65. All "_w" metrics use it; the
       "_raw" twins weigh every edge 1 and are reported as diagnostics.
    2. Hub set H: species with significantly more reactions than uniform attachment would give:
       d_s > c, c the smallest integer with P(Poisson(lambda) > c) <= 0.05 / |S|, lambda the mean
       degree of the species not in H (iterated to a fixed point). Used for the clone rule and for the
       "nohub" connectedness diagnostic.

Levels
    P_k(v) = the group of node v at level k, or its deepest group above k when v has no group at
    level k (a node placed directly in a parent that also has child groups forms a residual group
    with the parent's id; `residual_node_share` reports how many nodes that concerns).
    "leaf" = the deepest level used anywhere. Every per-level metric is suffixed _L1, _L2, _L3, _leaf;
    a level that is not used is NaN. A level with fewer than two groups separates nothing: cut,
    modularity, conductance and home-plurality are NaN there, and `judge` fails (does not skip)
    rules G1, G5, G6 when that level is the leaf.
    Combination over levels: G1, G2, G3, G4, G6, G8 are judged at the leaf; modularity (G5) takes
    the best level; the edge span (G7, Dasgupta's cost) integrates all levels in one number.

Clones
    A species with `clone_in` is drawn again in those groups. Only the cut (G1) and the edge span
    (G7) follow the drawing, and only for hub species: an edge (r, s) does not cross at level k when
    s is a hub and a copy of s lies in P_k(r). Copies of non-hub species never make an edge local.
    Every other metric uses the home position only, so cloning cannot improve it. G11 judges the
    clones themselves.

Returned keys (all floats)
    size      n_groups_*, largest_share_*, singleton_group_share_*, singleton_node_share_*,
              n_eff_groups_*, median_group_size_*, readable_node_share_*
    cut       cut_w_* (headline), cut_raw_*, cut_w_home_*, cut_raw_home_*, cut_w_mean, cut_raw_mean
    modular   Q_w_*, Q_raw_*, Q_w_best, Q_raw_best, Q_ref, n_groups_ref, Q_eff
    groups    connected_node_share_*, connected_nohub_node_share_*, disconnected_group_share_*,
              weak_node_share_*, conductance_mean_*, home_plurality_share_*
    tree      depth, edge_span_w, edge_span_null, edge_span_ratio, max_siblings, max_children,
              single_child_share, residual_node_share
    hubs      hub_degree_cutoff, n_hubs, hub_edge_share_raw
    clones    n_clones, clones_per_species, clone_unknown, clone_idle_share, clone_nonhub_share,
              clone_misuse_share, clone_absorbed_w
    reference ref_n_labelled, ref_coverage, ref_bcubed_{precision,recall,f1}_{L1,leaf},
              ref_bcubed_f1_best, ref_cut_w, ref_cut_raw      (NaN without a reference)
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np
import pandas as pd

LEVELS = ("L1", "L2", "L3")
NAN = float("nan")

# thresholds of rules_graph.md, in one place: (pass bound, warn bound). "lower" = lower is better.
READABLE_MIN, READABLE_MAX = 3, 50
SIBLINGS_PASS, SIBLINGS_WARN = 9, 50
THRESHOLDS = {
    "G1": ("cut_w_leaf", "lower", 0.25, 0.50),
    "G2": ("largest_share_leaf", "lower", 1 / 3, 1 / 2),
    "G3": ("singleton_group_share_leaf", "lower", 0.05, 0.20),
    "G4": ("connected_node_share_leaf", "higher", 0.95, 0.80),
    "G5": ("Q_eff", "higher", 0.70, 0.40),
    "G6": ("weak_node_share_leaf", "higher", 0.80, 0.50),
    "G7": ("edge_span_ratio", "lower", 1 / 3, 2 / 3),
    "G8": ("readable_node_share_leaf", "higher", 0.80, 0.50),
    "G9": ("max_siblings", "lower", SIBLINGS_PASS, SIBLINGS_WARN),
    "G10": ("single_child_share", "lower", 0.0, 0.20),
    "G11": ("clone_misuse_share", "lower", 0.10, 0.30),
    "G12": ("ref_bcubed_f1_best", "higher", 0.60, 0.40),
}
WEIGHTS = {"G1": 3, "G2": 2, "G3": 2, "G4": 2, "G5": 2, "G6": 1, "G7": 2, "G8": 1, "G9": 1, "G10": 1, "G11": 1}
GATES = ("G2", "G3")          # a grouping failing one of these cannot score above GATE_CAP
GATE_CAP = 0.40
# rules that are failed, not skipped, when the leaf level has fewer than two groups
FAIL_IF_ONE_GROUP = ("G1", "G5", "G6")


# ----------------------------------------------------------------------------------------------- graph
def hub_cutoff(deg, alpha: float = 0.05) -> int:
    """Degree above which a species is a hub: Poisson tail, family-wise alpha, lambda from non-hubs."""
    deg = np.asarray(deg, dtype=float)
    deg = deg[deg > 0]
    if deg.size == 0:
        return 0
    n, keep, cut = deg.size, deg, 0
    for _ in range(100):
        lam = float(keep.mean())
        a, p = alpha / n, math.exp(-lam)
        cdf, cut = p, 0
        while 1.0 - cdf > a and cut < 10000:
            cut += 1
            p *= lam / cut
            cdf += p
        new = deg[deg <= cut]
        if new.size == keep.size or new.size == 0:
            break
        keep = new
    return int(cut)


class _Graph:
    """Node index, distinct (reaction, species) pairs, degrees, weights and the hub set."""

    def __init__(self, nodes: pd.DataFrame, edges: pd.DataFrame):
        self.ids = list(nodes["node_id"])
        self.idx = {v: i for i, v in enumerate(self.ids)}
        if len(self.idx) != len(self.ids):
            raise ValueError("nodes.tsv: duplicate node_id")
        self.n = len(self.ids)
        self.is_species = (nodes["kind"].to_numpy() == "species")
        pairs = edges[["reaction_id", "species_id"]].drop_duplicates()
        pairs = pairs[pairs["reaction_id"].isin(self.idx) & pairs["species_id"].isin(self.idx)]
        self.r = pairs["reaction_id"].map(self.idx).to_numpy(dtype=np.int64)
        self.s = pairs["species_id"].map(self.idx).to_numpy(dtype=np.int64)
        self.m = len(self.r)
        self.deg = np.bincount(self.s, minlength=self.n) + np.bincount(self.r, minlength=self.n)
        self.w = 1.0 / self.deg[self.s] if self.m else np.zeros(0)
        self.raw = np.ones(self.m)
        sdeg = np.where(self.is_species, self.deg, 0)
        self.hub_cut = hub_cutoff(sdeg[self.is_species]) if self.is_species.any() else 0
        self.is_hub = self.is_species & (sdeg > self.hub_cut) if self.hub_cut > 0 else np.zeros(self.n, bool)


def _codes(labels) -> tuple[np.ndarray, int]:
    codes, uniq = pd.factorize(np.asarray(labels, dtype=object), sort=True)
    return codes.astype(np.int64), len(uniq)


def _barber_q(g: _Graph, w: np.ndarray, part: np.ndarray, k: int) -> float:
    """Barber's bipartite modularity: Q = sum_c [ I_c / W - KR_c * KS_c / W^2 ]."""
    W = float(w.sum())
    if W <= 0:
        return NAN
    same = part[g.r] == part[g.s]
    inner = np.bincount(part[g.r][same], weights=w[same], minlength=k)
    kr = np.bincount(part[g.r], weights=w, minlength=k)
    ks = np.bincount(part[g.s], weights=w, minlength=k)
    return float((inner / W - kr * ks / W ** 2).sum())


def barber_louvain(g: _Graph, w: np.ndarray) -> np.ndarray:
    """Deterministic Louvain-type optimiser of Barber's modularity (the reference value Q_ref).

    Nodes are visited in table order, ties keep the current community, then the lowest community id.
    Returns a community index per node.
    """
    W = float(w.sum())
    n = g.n
    member = np.arange(n)
    if W <= 0:
        return member
    kr = np.bincount(g.r, weights=w, minlength=n).astype(float)   # strength as a reaction
    ks = np.bincount(g.s, weights=w, minlength=n).astype(float)   # strength as a species
    adj = [dict() for _ in range(n)]
    for a, b, x in zip(g.r.tolist(), g.s.tolist(), w.tolist()):
        adj[a][b] = adj[a].get(b, 0.0) + x
        adj[b][a] = adj[b].get(a, 0.0) + x
    while True:
        nn = len(adj)
        comm = list(range(nn))
        KR, KS = kr.copy(), ks.copy()
        moved_any = False
        for _ in range(100):
            moved = False
            for i in range(nn):
                c0 = comm[i]
                KR[c0] -= kr[i]
                KS[c0] -= ks[i]
                to = {}
                for j, x in adj[i].items():
                    if j != i:
                        to[comm[j]] = to.get(comm[j], 0.0) + x
                best_c = c0
                best = to.get(c0, 0.0) - (kr[i] * KS[c0] + ks[i] * KR[c0]) / W
                for c in sorted(to):
                    gain = to[c] - (kr[i] * KS[c] + ks[i] * KR[c]) / W
                    if gain > best + 1e-12:
                        best, best_c = gain, c
                KR[best_c] += kr[i]
                KS[best_c] += ks[i]
                if best_c != c0:
                    comm[i] = best_c
                    moved = moved_any = True
            if not moved:
                break
        if not moved_any:
            break
        relabel = {c: k for k, c in enumerate(sorted(set(comm)))}
        comm = [relabel[c] for c in comm]
        k = len(relabel)
        member = np.asarray(comm)[member]
        kr2, ks2 = np.zeros(k), np.zeros(k)
        adj2 = [dict() for _ in range(k)]
        for i in range(nn):
            kr2[comm[i]] += kr[i]
            ks2[comm[i]] += ks[i]
            for j, x in adj[i].items():
                adj2[comm[i]][comm[j]] = adj2[comm[i]].get(comm[j], 0.0) + x
        kr, ks, adj = kr2, ks2, adj2
    return member


def _components(n: int, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Connected-component root of every node (union-find) for the edge list (a, b)."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for x, y in zip(a.tolist(), b.tolist()):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)
    return np.asarray([find(i) for i in range(n)], dtype=np.int64)


# ------------------------------------------------------------------------------------------- grouping
class _Grouping:
    """groups.tsv aligned to the node index: P_k, chains of group ids, clones."""

    def __init__(self, g: _Graph, groups: pd.DataFrame):
        seen = groups["node_id"]
        if seen.duplicated().any():
            raise ValueError("groups.tsv: a node appears more than once")
        missing = set(g.ids) - set(seen)
        if missing:
            raise ValueError(f"groups.tsv: {len(missing)} nodes of nodes.tsv are not grouped, e.g. {sorted(missing)[:3]}")
        gr = groups.set_index("node_id").reindex(g.ids)
        lv = {c: (gr[c].to_numpy(dtype=object) if c in gr.columns else np.full(g.n, "", dtype=object)) for c in LEVELS}
        if (lv["L1"] == "").any():
            raise ValueError("groups.tsv: a node has no L1")
        self.depth = max(i for i, c in enumerate(LEVELS, start=1) if (lv[c] != "").any())
        self.P = {}
        prev = lv["L1"]
        for i, c in enumerate(LEVELS, start=1):
            cur = np.where(lv[c] != "", lv[c], prev) if i > 1 else lv[c]
            self.P[i] = cur
            prev = cur
        # chain root..group for every group id that has members
        self.chain = {}
        for a, b, c in zip(lv["L1"], lv["L2"], lv["L3"]):
            ch = [x for x in (a, b, c) if x]
            for i in range(len(ch)):
                self.chain.setdefault(ch[i], tuple(ch[: i + 1]))
        self.lv = lv
        # clones: (species index, group id) with a known group
        self.clones, self.clone_unknown = [], 0
        col = gr["clone_in"].to_numpy(dtype=object) if "clone_in" in gr.columns else np.full(g.n, "", dtype=object)
        for i, txt in enumerate(col):
            if not txt:
                continue
            for c in sorted(set(filter(None, str(txt).split("|")))):
                if c in self.chain and g.is_species[i] and c not in self.chain_of_node(i):
                    self.clones.append((i, c))
                else:
                    self.clone_unknown += 1

    def chain_of_node(self, i: int) -> tuple:
        return tuple(x for x in (self.lv["L1"][i], self.lv["L2"][i], self.lv["L3"][i]) if x)

    def clone_location(self, c: str, k: int) -> str:
        """Group that holds a copy placed in group c, seen at level k."""
        ch = self.chain[c]
        return ch[min(k, len(ch)) - 1]


def _absorbed(g: _Graph, gp: _Grouping, k: int) -> np.ndarray:
    """Edges (r, s) that a copy of the hub species s keeps inside P_k(r)."""
    out = np.zeros(g.m, dtype=bool)
    if not gp.clones:
        return out
    where = {}
    for i, c in gp.clones:
        if g.is_hub[i]:                       # a copy of a non-hub species never makes an edge local
            where.setdefault(i, set()).add(gp.clone_location(c, k))
    P = gp.P[k]
    for e in range(g.m):
        loc = where.get(int(g.s[e]))
        if loc and P[g.r[e]] in loc:
            out[e] = True
    return out


def _level_metrics(g: _Graph, gp: _Grouping, k: int) -> dict:
    """All per-level metrics for level k (suffix added by the caller)."""
    part, ng = _codes(gp.P[k])
    sizes = np.bincount(part, minlength=ng)
    out = {
        "n_groups": float(ng),
        "largest_share": float(sizes.max() / g.n),
        "singleton_group_share": float((sizes == 1).sum() / ng),
        "singleton_node_share": float((sizes == 1).sum() / g.n),
        "n_eff_groups": float(g.n ** 2 / (sizes.astype(float) ** 2).sum()),
        "readable_node_share": float(sizes[(sizes >= READABLE_MIN) & (sizes <= READABLE_MAX)].sum() / g.n),
        "median_group_size": float(np.median(sizes)),
    }
    # connectedness of every group through its own home members
    same = part[g.r] == part[g.s]
    for tag, keep in (("connected", same), ("connected_nohub", same & ~g.is_hub[g.s])):
        root = _components(g.n, g.r[keep], g.s[keep])
        use = np.ones(g.n, bool) if tag == "connected" else ~g.is_hub
        big = sizes >= 2
        comp = pd.DataFrame({"g": part[use], "c": root[use]}).value_counts().reset_index(name="n")
        largest = comp.groupby("g")["n"].max().reindex(range(ng), fill_value=0).to_numpy()
        total = np.bincount(part[use], minlength=ng)
        out[f"{tag}_node_share"] = float(largest[big].sum() / total[big].sum()) if total[big].sum() else NAN
        if tag == "connected":
            out["disconnected_group_share"] = float((largest[big] < sizes[big]).sum() / big.sum()) if big.any() else NAN
    sep = ("cut_w", "cut_raw", "cut_w_home", "cut_raw_home", "Q_w", "Q_raw", "weak_node_share",
           "conductance_mean", "home_plurality_share")
    if ng < 2 or g.m == 0:
        out.update({key: NAN for key in sep})
        return out
    # cut: share of edge weight whose ends lie in different groups (clones respected / ignored)
    cross = ~same
    drawn = cross & ~_absorbed(g, gp, k)
    out["cut_w"] = float(g.w[drawn].sum() / g.w.sum())
    out["cut_raw"] = float(drawn.mean())
    out["cut_w_home"] = float(g.w[cross].sum() / g.w.sum())
    out["cut_raw_home"] = float(cross.mean())
    out["Q_w"] = _barber_q(g, g.w, part, ng)
    out["Q_raw"] = _barber_q(g, g.raw, part, ng)
    # conductance and the weak-community test per group (home position, species-normalised weight)
    W = float(g.w.sum())
    inner = np.bincount(part[g.r][same], weights=g.w[same], minlength=ng)
    vol = np.bincount(part[g.r], weights=g.w, minlength=ng) + np.bincount(part[g.s], weights=g.w, minlength=ng)
    cut = vol - 2 * inner
    den = np.minimum(vol, 2 * W - vol)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(den > 1e-12, cut / den, np.nan)
    weak = 2 * inner > cut + 1e-12
    out["weak_node_share"] = float(sizes[weak].sum() / g.n)
    ok = ~np.isnan(phi)
    out["conductance_mean"] = float((phi[ok] * sizes[ok]).sum() / sizes[ok].sum()) if ok.any() else NAN
    # a species sits in the group that holds most of its reactions (species with >= 2 reactions)
    multi = g.is_species & (g.deg >= 2)
    if multi.any():
        tab = pd.DataFrame({"s": g.s, "grp": part[g.r]}).value_counts().reset_index(name="n")
        top = tab.groupby("s")["n"].max()
        tab["home"] = tab["grp"].to_numpy() == part[tab["s"].to_numpy()]
        at_home = tab[tab["home"]].set_index("s")["n"]
        sp = np.flatnonzero(multi)
        good = at_home.reindex(sp, fill_value=0).to_numpy() >= top.reindex(sp, fill_value=0).to_numpy()
        out["home_plurality_share"] = float(good.mean())
    else:
        out["home_plurality_share"] = NAN
    return out


def _edge_span(g: _Graph, gp: _Grouping) -> tuple[float, float]:
    """Dasgupta's cost of the hierarchy, normalised: (span, span of a random filling of the same tree).

    span = sum_e w_e * |box(e)| / N / W, box(e) the smallest group that holds both ends of edge e
    (a copy of a hub species counts as an end), the whole model (N nodes) when no group does.
    The null replaces the edges by all species-reaction pairs, i.e. nodes dealt at random into the
    same boxes.
    """
    full, n_r, n_s = {}, {}, {}
    for i in range(g.n):
        for x in gp.chain_of_node(i):
            full[x] = full.get(x, 0) + 1
            d = n_s if g.is_species[i] else n_r
            d[x] = d.get(x, 0) + 1
    box = np.full(g.m, float(g.n))
    done = np.zeros(g.m, dtype=bool)
    for k in range(gp.depth, 0, -1):
        P = gp.P[k]
        together = ((P[g.r] == P[g.s]) | _absorbed(g, gp, k)) & ~done
        if together.any():
            box[together] = [full[x] for x in P[g.r[together]]]
            done |= together
    W = float(g.w.sum())
    span = float((g.w * box).sum() / (g.n * W)) if W > 0 else NAN
    # null: pairs (reaction, species) whose smallest common box is x = nR_x * nS_x - sum over children
    kids = {}
    for ch in {gp.chain_of_node(i) for i in range(g.n)}:
        for a, b in zip(ch, ch[1:]):
            kids.setdefault(a, set()).add(b)
    n_rx, n_sp = int((~g.is_species).sum()), int(g.is_species.sum())
    if n_rx == 0 or n_sp == 0:
        return span, NAN
    tot, inside_top = 0.0, 0.0
    for x in sorted(full):
        own = n_r.get(x, 0) * n_s.get(x, 0) - sum(n_r.get(c, 0) * n_s.get(c, 0) for c in kids.get(x, ()))
        tot += own * full[x] / g.n
        if len(gp.chain[x]) == 1:
            inside_top += n_r.get(x, 0) * n_s.get(x, 0)
    tot += (n_rx * n_sp - inside_top) * 1.0
    return span, float(tot / (n_rx * n_sp))


def _bcubed(part: np.ndarray, labels: list[set]) -> tuple[float, float, float]:
    """Extended BCubed precision, recall, F1 for a partition against overlapping reference labels."""
    idx = [i for i, l in enumerate(labels) if l]
    if len(idx) < 2:
        return NAN, NAN, NAN
    grp = {}
    for i in idx:
        grp.setdefault(int(part[i]), []).append(i)
    lab = {}
    for i in idx:
        for l in labels[i]:
            lab.setdefault(l, []).append(i)
    prec = []
    for i in idx:                      # items sharing my group: do they share a label with me?
        mates = grp[int(part[i])]
        prec.append(sum(1 for j in mates if labels[i] & labels[j]) / len(mates))
    rec = []
    for i in idx:                      # items sharing a label with me: are they in my group?
        mates = set()
        for l in labels[i]:
            mates.update(lab[l])
        rec.append(sum(1 for j in mates if part[j] == part[i]) / len(mates))
    p, r = float(np.mean(prec)), float(np.mean(rec))
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)


def _reference_metrics(g: _Graph, gp: _Grouping, reference: pd.DataFrame) -> dict:
    out = {}
    labels = [set() for _ in range(g.n)]
    for v, l in zip(reference["node_id"], reference["label"]):
        if v in g.idx and l != "":
            labels[g.idx[v]].add(l)
    n_lab = sum(1 for l in labels if l)
    out["ref_n_labelled"] = float(n_lab)
    out["ref_coverage"] = float(n_lab / g.n)
    best = NAN
    for name, k in (("L1", 1), ("leaf", gp.depth)):
        part, _ = _codes(gp.P[k])
        p, r, f = _bcubed(part, labels)
        out[f"ref_bcubed_precision_{name}"], out[f"ref_bcubed_recall_{name}"], out[f"ref_bcubed_f1_{name}"] = p, r, f
        if not math.isnan(f) and (math.isnan(best) or f > best):
            best = f
    out["ref_bcubed_f1_best"] = best
    # calibration: the cut share of the curated reference itself. An unlabelled reaction takes the
    # most frequent label of its labelled participants (ties: first in sorted order); an edge is
    # internal when both ends share a label. Edges with an unlabelled end are left out.
    lab2 = [set(l) for l in labels]
    votes = {}
    for a, b in zip(g.r.tolist(), g.s.tolist()):
        if not labels[a]:
            for l in labels[b]:
                votes.setdefault(a, {}).setdefault(l, 0)
                votes[a][l] += 1
    for a, v in votes.items():
        lab2[a] = {min(v, key=lambda l: (-v[l], l))}
    use = np.asarray([bool(lab2[a]) and bool(labels[b]) for a, b in zip(g.r.tolist(), g.s.tolist())], dtype=bool)
    if use.any():
        cross = np.asarray([not (lab2[a] & labels[b]) for a, b in zip(g.r.tolist(), g.s.tolist())], dtype=bool)
        out["ref_cut_w"] = float(g.w[use & cross].sum() / g.w[use].sum())
        out["ref_cut_raw"] = float(cross[use].mean())
    else:
        out["ref_cut_w"] = out["ref_cut_raw"] = NAN
    return out


# -------------------------------------------------------------------------------------------- metrics
def metrics(nodes, edges, comps, groups, meta, reference=None) -> dict[str, float]:
    """All graph metrics of one grouping. NaN where a metric does not apply.

    nodes, edges, comps: the feature tables of the model; groups, meta: groups.tsv and groups_meta.tsv
    of the strategy; all read with dtype=str, keep_default_na=False. `comps` and `meta` are accepted
    for a uniform signature; the graph metrics need neither (the hierarchy is read from groups.tsv).
    reference: optional table (node_id, label), several rows per node allowed, unlabelled nodes absent.
    """
    g = _Graph(nodes, edges)
    gp = _Grouping(g, groups)
    out = {
        "n_nodes": float(g.n), "n_species": float(g.is_species.sum()), "n_reactions": float((~g.is_species).sum()),
        "n_edges": float(g.m), "depth": float(gp.depth),
        "hub_degree_cutoff": float(g.hub_cut), "n_hubs": float(g.is_hub.sum()),
        "hub_edge_share_raw": float(g.is_hub[g.s].mean()) if g.m else NAN,
    }
    per = {}
    for k in (1, 2, 3):
        per[k] = _level_metrics(g, gp, k) if k <= gp.depth else None
    keys = list(per[1])
    for k in (1, 2, 3):
        for key in keys:
            out[f"{key}_L{k}"] = per[k][key] if per[k] else NAN
    for key in keys:
        out[f"{key}_leaf"] = per[gp.depth][key]
    # combination over levels
    used = [per[k] for k in range(1, gp.depth + 1)]
    for key in ("cut_w", "cut_raw", "cut_w_home"):          # diagnostics; G1 is judged at the leaf
        vals = [p[key] for p in used if not math.isnan(p[key])]
        out[f"{key}_mean"] = float(np.mean(vals)) if vals else NAN
    for key in ("Q_w", "Q_raw"):
        vals = [p[key] for p in used if not math.isnan(p[key])]
        out[f"{key}_best"] = float(max(vals)) if vals else NAN
    ref_part = barber_louvain(g, g.w)
    out["Q_ref"] = _barber_q(g, g.w, *_codes(ref_part)) if g.m else NAN
    out["n_groups_ref"] = float(len(set(ref_part.tolist())))
    span, span_null = _edge_span(g, gp)
    out["edge_span_w"], out["edge_span_null"] = span, span_null
    out["edge_span_ratio"] = float(span / span_null) if not math.isnan(span) and span_null and span_null > 0 else NAN
    q, qr = out["Q_w_best"], out["Q_ref"]
    out["Q_eff"] = float(q / qr) if not math.isnan(q) and not math.isnan(qr) and qr > 1e-9 else NAN
    # shape of the hierarchy
    kids = {}
    for ch in {gp.chain_of_node(i) for i in range(g.n)}:
        for a, b in zip(ch, ch[1:]):
            kids.setdefault(a, set()).add(b)
    direct = {}
    for i in range(g.n):
        direct[gp.chain_of_node(i)[-1]] = direct.get(gp.chain_of_node(i)[-1], 0) + 1
    sib = [len(set(gp.lv["L1"]))] + [len(c) + (1 if direct.get(p) else 0) for p, c in sorted(kids.items())]
    out["residual_node_share"] = float(sum(n for p, n in direct.items() if p in kids) / g.n)
    out["max_siblings"] = float(max(sib))
    out["max_children"] = float(max(sib[1:])) if len(sib) > 1 else NAN
    # a parent whose only content is one child group adds a box and no information; the (undrawn) root
    # is counted as such a parent when a nested grouping has a single top-level group
    single = [1.0 if (len(c) == 1 and not direct.get(p)) else 0.0 for p, c in sorted(kids.items())]
    if single and sib[0] == 1:
        single.append(1.0)
    out["single_child_share"] = float(np.mean(single)) if single else NAN
    # clones
    n_sp = max(int(g.is_species.sum()), 1)
    out["n_clones"] = float(len(gp.clones))
    out["clones_per_species"] = float(len(gp.clones) / n_sp)
    out["clone_unknown"] = float(gp.clone_unknown)
    if gp.clones:
        in_group = {}
        for e in range(g.m):
            in_group.setdefault(int(g.s[e]), set()).update(gp.chain_of_node(int(g.r[e])))
        idle = np.asarray([c not in in_group.get(i, ()) for i, c in gp.clones])
        nonhub = np.asarray([not g.is_hub[i] for i, _ in gp.clones])
        out["clone_idle_share"] = float(idle.mean())
        out["clone_nonhub_share"] = float(nonhub.mean())
        out["clone_misuse_share"] = float((idle | nonhub).mean())
    else:
        out["clone_idle_share"] = out["clone_nonhub_share"] = out["clone_misuse_share"] = NAN
    a, b = out["cut_w_home_mean"], out["cut_w_mean"]
    out["clone_absorbed_w"] = float(a - b) if not math.isnan(a) and not math.isnan(b) else NAN
    if reference is not None and len(reference):
        out.update(_reference_metrics(g, gp, reference))
    else:
        out.update({key: NAN for key in (
            "ref_n_labelled", "ref_coverage", "ref_bcubed_precision_L1", "ref_bcubed_recall_L1", "ref_bcubed_f1_L1",
            "ref_bcubed_precision_leaf", "ref_bcubed_recall_leaf", "ref_bcubed_f1_leaf", "ref_bcubed_f1_best",
            "ref_cut_w", "ref_cut_raw")})
    return out


# ---------------------------------------------------------------------------------------------- rules
def judge(m: dict) -> dict[str, str]:
    """Verdict of every rule: 'pass', 'warn', 'fail' or 'na' (rule does not apply)."""
    out = {}
    one_group = m["n_groups_leaf"] < 2
    for rule, (key, direction, p, w) in THRESHOLDS.items():
        v = m.get(key, NAN)
        if rule == "G9":      # both the widest sibling set and "more than one top-level group"
            if m["n_groups_L1"] < 2:
                out[rule] = "fail"
                continue
        if v is None or math.isnan(v):
            out[rule] = "fail" if (one_group and rule in FAIL_IF_ONE_GROUP) else "na"
            continue
        if direction == "lower":
            out[rule] = "pass" if v <= p + 1e-12 else "warn" if v <= w + 1e-12 else "fail"
        else:
            out[rule] = "pass" if v >= p - 1e-12 else "warn" if v >= w - 1e-12 else "fail"
    return out


def score(m: dict, verdicts: dict | None = None) -> float:
    """One number in [0, 1]: weighted mean of pass=1 / warn=0.5 / fail=0 over the applicable rules
    G1..G11, capped at GATE_CAP when a gate rule (largest group, singletons) fails."""
    v = verdicts or judge(m)
    pts = {"pass": 1.0, "warn": 0.5, "fail": 0.0}
    num = sum(WEIGHTS[r] * pts[v[r]] for r in WEIGHTS if v[r] != "na")
    den = sum(WEIGHTS[r] for r in WEIGHTS if v[r] != "na")
    total = num / den if den else NAN
    if any(v[r] == "fail" for r in GATES):
        total = min(total, GATE_CAP)
    return float(total)


# ------------------------------------------------------------------------------- synthetic groupings
def _frame(ids, l1, l2=None, clone=None) -> pd.DataFrame:
    n = len(ids)
    return pd.DataFrame({"node_id": ids, "L1": l1, "L2": l2 if l2 is not None else [""] * n, "L3": [""] * n,
                         "rule": ["R1"] * n, "clone_in": clone if clone is not None else [""] * n})


def _meta(groups: pd.DataFrame) -> pd.DataFrame:
    rows = {}
    for a, b, c in zip(groups["L1"], groups["L2"], groups["L3"]):
        up = ""
        for i, x in enumerate((a, b, c), start=1):
            if x:
                rows[x] = (x, str(i), up, "", "topology", "none", "", "")
                up = x
    return pd.DataFrame(sorted(rows.values()), columns=["group_id", "level", "parent", "label", "basis",
                                                        "name_source", "name_xref", "compartment_id"])


def _coarsen(g: _Graph, part: np.ndarray, target: int) -> np.ndarray:
    """Merge groups greedily (largest gain in Barber modularity, connected pairs only) down to `target`."""
    part = np.asarray(part).copy()
    W = float(g.w.sum())
    while len(set(part.tolist())) > target:
        codes, k = _codes(part)
        kr = np.bincount(codes[g.r], weights=g.w, minlength=k)
        ks = np.bincount(codes[g.s], weights=g.w, minlength=k)
        between = {}
        for a, b, x in zip(codes[g.r].tolist(), codes[g.s].tolist(), g.w.tolist()):
            if a != b:
                key = (min(a, b), max(a, b))
                between[key] = between.get(key, 0.0) + x
        if not between:
            break
        best = max(sorted(between), key=lambda p: between[p] / W - (kr[p[0]] * ks[p[1]] + kr[p[1]] * ks[p[0]]) / W ** 2)
        codes[codes == best[1]] = best[0]
        part = codes
    return part


def synthetic(nodes: pd.DataFrame, edges: pd.DataFrame, seed: int = 0) -> dict[str, pd.DataFrame]:
    """Synthetic groupings of a model, for testing the metrics: name -> groups table."""
    ids = list(nodes["node_id"])
    g = _Graph(nodes, edges)
    opt = barber_louvain(g, g.w)
    opt_l = [f"m{c:03d}" for c in opt]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(ids))
    rand_l = [None] * len(ids)
    for src, dst in enumerate(perm):          # same sizes as the optimiser's groups, members shuffled
        rand_l[dst] = "r" + opt_l[src][1:]
    comp = ["c_" + (x or "none") for x in nodes["loc_space"]]
    # adversarial cases
    first = {}                                   # species -> its first reaction in table order
    for a, b in zip(g.r.tolist(), g.s.tolist()):
        first.setdefault(b, a)
    per_rx = [f"x{(first.get(i, i) if g.is_species[i] else i):05d}" for i in range(g.n)]
    used = {}                                    # species -> groups (random grouping) of its reactions
    for a, b in zip(g.r.tolist(), g.s.tolist()):
        used.setdefault(b, set()).add(rand_l[a])
    clone_all = ["|".join(sorted(used.get(i, set()) - {rand_l[i]})) if g.is_species[i] else "" for i in range(g.n)]
    hub_only = ["|".join(sorted({opt_l[a] for a, b in zip(g.r.tolist(), g.s.tolist()) if b == i} - {opt_l[i]}))
                if g.is_hub[i] else "" for i in range(g.n)]
    split = [("sp_" + (c or "none")) if sp else ("rx_" + (o or "none"))
             for sp, c, o in zip(g.is_species, nodes["loc_space"], nodes["sbo"])]
    top = _coarsen(g, opt, SIBLINGS_PASS)
    top_l = [f"t{c:03d}" for c in top]
    n_kids = pd.Series(opt_l).groupby(pd.Series(top_l)).nunique()
    nested_l2 = [o if n_kids[t] > 1 else "" for t, o in zip(top_l, opt_l)]   # no box-in-box for a lone child
    return {
        "optimiser_nested": _frame(ids, top_l, nested_l2),
        "per_reaction": _frame(ids, per_rx),
        "random_clone_all": _frame(ids, rand_l, clone=clone_all),
        "optimiser_clone_hubs": _frame(ids, opt_l, clone=hub_only),
        "species_apart_from_reactions": _frame(ids, split),
        "one_group": _frame(ids, ["all"] * len(ids)),
        "singletons": _frame(ids, [f"s{i:05d}" for i in range(len(ids))]),
        "compartment": _frame(ids, comp),
        "random_sized": _frame(ids, rand_l),
        "optimiser": _frame(ids, opt_l),
        "compartment_optimiser": _frame(ids, comp, [f"{c}.{o}" for c, o in zip(comp, opt_l)]),
        "dummy_root_optimiser": _frame(ids, ["all"] * len(ids), opt_l),
    }


def _read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Graph-theoretic metrics of a grouping (rules_graph.md).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score")
    s.add_argument("--features", required=True)
    s.add_argument("--grouping", required=True)
    s.add_argument("--reference")
    s.add_argument("--json")
    d = sub.add_parser("demo")
    d.add_argument("--features", required=True)
    a = ap.parse_args(argv)
    f = pathlib.Path(a.features)
    nodes, edges, comps = _read(f / "nodes.tsv"), _read(f / "edges.tsv"), _read(f / "compartments.tsv")
    if a.cmd == "score":
        o = pathlib.Path(a.grouping)
        ref = _read(a.reference) if a.reference else None
        m = metrics(nodes, edges, comps, _read(o / "groups.tsv"), _read(o / "groups_meta.tsv"), ref)
        v = judge(m)
        for rule, (key, *_rest) in THRESHOLDS.items():
            print(f"{rule:4s} {v[rule]:5s} {key} = {m.get(key, NAN):.4g}")
        print(f"graph score = {score(m, v):.3f}")
        if a.json:
            pathlib.Path(a.json).write_text(json.dumps({"metrics": m, "verdicts": v, "score": score(m, v)}, indent=1),
                                            encoding="utf-8")
            print(f"all {len(m)} metrics -> {a.json}")
    else:
        rows = []
        for name, grp in synthetic(nodes, edges).items():
            m = metrics(nodes, edges, comps, grp, _meta(grp))
            v = judge(m)
            rows.append({"grouping": name, "score": round(score(m, v), 3), **v})
        print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
