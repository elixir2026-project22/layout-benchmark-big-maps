"""Run:
    python metrics_biology.py panel/features/MODEL strategies/NAME/out/MODEL
    python metrics_biology.py panel/features/MODEL strategies/NAME/out/MODEL --reference REF.tsv --rerun SECOND_RUN_DIR --json
    pytest -q test_metrics_biology.py

Biology and FAIR metrics for a nested grouping of an SBML model (BioHackathon Europe 2026, project 22).

The first argument is the folder with nodes.tsv, edges.tsv, compartments.tsv; the second the folder with
groups.tsv and groups_meta.tsv (CONTRACT.md). REF.tsv is the optional held-out curated reference with the
columns node_id, label (a node may have several rows, a node may be absent). SECOND_RUN_DIR is the output of
a second run of the same strategy on the same model; it is compared byte for byte (rule FAIR-1).

Python API:
    metrics(nodes, edges, comps, groups, meta, reference=None, rerun=None) -> dict[str, float]
    judge(m)  -> {rule id: "pass" | "warn" | "fail" | "na"}
    score(m)  -> {"score": float in [0, 1], "block_<name>": float, ...}
All tables are pandas DataFrames read with dtype=str, keep_default_na=False. A metric that does not apply
is NaN. Nothing here uses a random number, the network, or a file outside the arguments: the same tables give
the same numbers. Only numpy and pandas are needed.

Rule ids and the metric each rule reads are in RULES below; the rules are explained in rules_biology.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

NAN = float("nan")
LEVELS = ("L1", "L2", "L3")
# BioModels.net qualifiers that link an element to a class, a whole or a place - not to itself.
CLASS_QUALIFIERS = {"bqbiol:isVersionOf", "bqbiol:isPartOf", "bqbiol:occursIn", "bqbiol:hasProperty",
                    "bqbiol:isPropertyOf"}
CHANNELS = ("sbo", "ec", "go", "class")
_CODE = re.compile(r"^(\d+(?:\.\d+)*)_|^(\d+(?:\.\d+)+)(?=\s|$)")

# One row per rule, most important first. passes / warns are the thresholds on the metric; worst / best its
# end points; block and weight say where the rule enters score().
RULES = [
    dict(id="BIO-1", metric="ref_skill", direction="higher", passes=0.25, warns=0.05, worst=0.0, best=1.0,
         block="reference", weight=1.0),
    dict(id="LOC-1", metric="loc_expressed", direction="higher", passes=1.0, warns=0.90, worst=0.0, best=1.0,
         block="localisation", weight=0.25),
    dict(id="BIO-2", metric="hom", direction="higher", passes=0.21, warns=0.05, worst=0.0, best=1.0,
         block="homogeneity", weight=1.0),
    dict(id="LOC-2", metric="loc_repeat", direction="lower", passes=1.0, warns=2.0, worst=math.inf, best=1.0,
         block="localisation", weight=0.25),
    dict(id="LOC-3", metric="loc_boundary", direction="higher", passes=1.0, warns=0.80, worst=0.0, best=1.0,
         block="localisation", weight=0.20),
    dict(id="LOC-4", metric="loc_nesting", direction="higher", passes=1.0, warns=0.90, worst=0.0, best=1.0,
         block="localisation", weight=0.15),
    dict(id="FAIR-1", metric="fair_deterministic", direction="higher", passes=1.0, warns=1.0, worst=0.0, best=1.0,
         block="gate", weight=0.0),
    dict(id="FAIR-2", metric="name_lift", direction="higher", passes=0.61, warns=0.21, worst=0.0, best=1.0,
         block="fair", weight=0.40),
    dict(id="LOC-5", metric="loc_span_ok", direction="higher", passes=1.0, warns=0.90, worst=0.0, best=1.0,
         block="localisation", weight=0.15),
    dict(id="FAIR-3", metric="name_grounded", direction="higher", passes=0.80, warns=0.50, worst=0.0, best=1.0,
         block="fair", weight=0.30),
    dict(id="FAIR-4", metric="stranded_poor", direction="lower", passes=0.10, warns=0.30, worst=1.0, best=0.0,
         block="fair", weight=0.30),
]
BLOCK_WEIGHTS = {"reference": 0.35, "homogeneity": 0.25, "localisation": 0.25, "fair": 0.15}
# utility of a metric value in score(): 0 at worst, 0.25 at the warn threshold, 0.75 at the pass threshold, 1 at best
UTILITY_AT_WARN, UTILITY_AT_PASS = 0.25, 0.75


# ----------------------------------------------------------------------------------------------- helpers
def _isnan(x) -> bool:
    return isinstance(x, float) and math.isnan(x)


def _chains(groups: pd.DataFrame) -> dict[str, list[str]]:
    """node id -> its group ids from L1 down to the deepest level it has."""
    out = {}
    for row in groups.itertuples(index=False):
        out[row.node_id] = [g for g in (row.L1, row.L2, row.L3) if g]
    return out


def _clones(groups: pd.DataFrame) -> dict[str, list[str]]:
    out = {}
    if "clone_in" in groups.columns:
        for nid, c in zip(groups["node_id"], groups["clone_in"]):
            cl = [g for g in c.split("|") if g]
            if cl:
                out[nid] = cl
    return out


class _Tree:
    """The compartment containment tree of compartments.tsv (same reading as common/compartments.Tree)."""

    def __init__(self, comps: pd.DataFrame):
        self.rows = {r["compartment_id"]: r for r in comps.to_dict("records")}

    def space_of(self, cid: str) -> str:
        r = self.rows.get(cid)
        if r is not None and r.get("kind") == "membrane" and r.get("encloses_id"):
            return r["encloses_id"]
        return cid

    def ancestors(self, cid: str) -> list[str]:
        """Proper ancestors of a space: parent, ..., root."""
        out, cur = [], cid
        while cur in self.rows and self.rows[cur].get("parent_id") and self.rows[cur]["parent_id"] not in out \
                and self.rows[cur]["parent_id"] != cid:
            cur = self.rows[cur]["parent_id"]
            out.append(cur)
        return out


# ------------------------------------------------------------------------------ ontology terms of a node
def _ec_terms(ec: str) -> set[str]:
    out = set()
    for code in filter(None, ec.split("|")):
        parts = [p for p in code.strip().split(".") if p and p != "-"]
        for i in range(1, len(parts) + 1):
            out.add("ec:" + ".".join(parts[:i]))
    return out


def _note_code_terms(notes: str) -> set[str]:
    """Hierarchical codes in the SBML notes ("key: 14.5.7_text; 18.4_text") as key:14, key:14.5, key:14.5.7."""
    out = set()
    if not notes:
        return out
    try:
        obj = json.loads(notes)
    except ValueError:
        return out
    for key, value in obj.items():
        if key == "text" or not isinstance(value, str):
            continue
        for piece in re.split(r"[;|]", value):
            m = _CODE.match(piece.strip())
            if m:
                parts = (m.group(1) or m.group(2)).split(".")
                for i in range(1, len(parts) + 1):
                    out.add(f"notes:{key.strip().lower()}:" + ".".join(parts[:i]))
    return out


def _class_terms(xrefs: str, notes: str) -> set[str]:
    out = set()
    for item in filter(None, xrefs.split("|")):
        qual, _, res = item.partition("=")
        if qual in CLASS_QUALIFIERS and res:
            db = res.split(":", 1)[0].lower()
            if db not in ("ec-code", "ec", "go"):
                out.add(res)
    return out | _note_code_terms(notes)


def node_terms(nodes: pd.DataFrame) -> dict[str, dict[str, frozenset]]:
    """channel -> node id -> ancestor-closed set of ontology terms (only nodes that have a term).

    sbo    reactions only: the SBO term and its ancestors (sbo_path). The SBO term of a species is its
           material type, and a pathway group rightly mixes metabolites and proteins.
    ec     EC numbers as prefixes: 1.2.1.105 -> ec:1, ec:1.2, ec:1.2.1, ec:1.2.1.105
    go     GO ids, flat (the bundle has no GO hierarchy)
    class  cross-references with a class qualifier (isVersionOf, isPartOf, occursIn, hasProperty,
           isPropertyOf) other than EC and GO, flat; and hierarchical codes in the notes, as prefixes
    """
    out = {c: {} for c in CHANNELS}
    for r in nodes.itertuples(index=False):
        if r.kind == "reaction" and r.sbo_path:
            out["sbo"][r.node_id] = frozenset(t for t in r.sbo_path.split(">") if t)
        t = _ec_terms(r.ec)
        if t:
            out["ec"][r.node_id] = frozenset(t)
        t = {g for g in r.go.split("|") if g}
        if t:
            out["go"][r.node_id] = frozenset(t)
        t = _class_terms(r.xrefs, r.notes)
        if t:
            out["class"][r.node_id] = frozenset(t)
    return out


class _Similarity:
    """Resnik similarity between nodes, per channel and combined (maximum over channels).

    IC(t) = ln(N / n_t) / ln(N): N nodes eligible for the channel (reactions for sbo, all nodes otherwise),
    n_t of them carry t. sim(i, j) = max IC(t) over the terms the two nodes share; 0 when they share none.
    Nodes with the same term sets are pooled into one signature, so the cost is in signatures, not nodes.
    """

    def __init__(self, nodes: pd.DataFrame):
        terms = node_terms(nodes)
        n_all = len(nodes)
        n_rx = int((nodes["kind"] == "reaction").sum())
        self.channel = {}
        for ch in CHANNELS:
            eligible = n_rx if ch == "sbo" else n_all
            per_node = terms[ch]
            sets = sorted(set(per_node.values()), key=lambda s: sorted(s))
            index = {s: i for i, s in enumerate(sets)}
            count = Counter()
            for s in per_node.values():
                count.update(s)
            sim = np.zeros((len(sets), len(sets)))
            if eligible > 1:
                holders = defaultdict(list)
                for s, i in index.items():
                    for t in s:
                        holders[t].append(i)
                for t in sorted(holders):
                    ic = math.log(eligible / count[t]) / math.log(eligible)
                    if ic <= 0:
                        continue
                    idx = np.array(holders[t])
                    block = sim[np.ix_(idx, idx)]
                    np.maximum(block, ic, out=block)
                    sim[np.ix_(idx, idx)] = block
            self.channel[ch] = dict(of={n: index[s] for n, s in per_node.items()}, sim=sim, n=len(per_node))
        self.annotated = sorted({n for ch in CHANNELS for n in terms[ch]})

    def matrices(self, channels=CHANNELS):
        """(node id -> signature index, similarity matrix M, comparable matrix K) over the given channels."""
        sig_of, sigs = {}, {}
        for n in self.annotated:
            key = tuple(self.channel[ch]["of"].get(n, -1) for ch in channels)
            if all(k < 0 for k in key):
                continue
            sig_of[n] = sigs.setdefault(key, len(sigs))
        keys = np.array(sorted(sigs, key=lambda k: sigs[k]), dtype=int).reshape(len(sigs), len(channels))
        m = np.zeros((len(sigs), len(sigs)))
        k = np.zeros((len(sigs), len(sigs)), dtype=bool)
        for c, ch in enumerate(channels):
            col = keys[:, c]
            has = col >= 0
            if not has.any():
                continue
            idx = np.where(has)[0]
            sub = self.channel[ch]["sim"][np.ix_(col[idx], col[idx])]
            block = m[np.ix_(idx, idx)]
            np.maximum(block, sub, out=block)
            m[np.ix_(idx, idx)] = block
            k[np.ix_(idx, idx)] = True
        return sig_of, m, k


def _homogeneity(sig_of: dict, m: np.ndarray, k: np.ndarray, part: dict[str, str]) -> float:
    """(S_in - S_all) / (1 - S_all); see rule BIO-2."""
    ids = [n for n in sig_of if n in part]
    if len(ids) < 2:
        return NAN
    sig = np.array([sig_of[n] for n in ids])
    glabels = sorted({part[n] for n in ids})
    gindex = {g: i for i, g in enumerate(glabels)}
    grp = np.array([gindex[part[n]] for n in ids])
    n_sig = m.shape[0]
    counts = np.zeros((len(glabels), n_sig))
    np.add.at(counts, (grp, sig), 1.0)
    kf = k.astype(float)
    self_m, self_k = np.diag(m), np.diag(kf)
    total = counts.sum(axis=0)
    all_num = total @ m - self_m          # per signature, against every other annotated node
    all_den = total @ kf - self_k
    in_num = counts @ m - self_m[None, :]  # per group and signature, against the other members
    in_den = counts @ kf - self_k[None, :]
    a_num, a_den = all_num[sig], all_den[sig]
    use = a_den > 0
    if not use.any():
        return NAN
    base = np.where(use, a_num / np.where(a_den > 0, a_den, 1.0), 0.0)
    i_num, i_den = in_num[grp, sig], in_den[grp, sig]
    inside = np.where(i_den > 0, i_num / np.where(i_den > 0, i_den, 1.0), base)  # nobody to compare with: chance
    s_in, s_all = float(inside[use].mean()), float(base[use].mean())
    if s_all >= 1.0:
        return NAN
    return (s_in - s_all) / (1.0 - s_all)


# --------------------------------------------------------------------------------- reference (BCubed)
def _bcubed(c: np.ndarray, lab: np.ndarray, chunk: int = 512) -> tuple[float, float, float]:
    """Extended BCubed precision, recall, F1 (Amigo et al. 2009) for item-by-cluster and item-by-label
    0/1 matrices over the same items. Every item has at least one cluster and one label."""
    n = c.shape[0]
    c = c.astype(np.float32)
    lab = lab.astype(np.float32)
    p_sum = r_sum = 0.0
    for a in range(0, n, chunk):
        cc = c[a:a + chunk] @ c.T
        ll = lab[a:a + chunk] @ lab.T
        mn = np.minimum(cc, ll)
        with np.errstate(divide="ignore", invalid="ignore"):
            pm = np.where(cc > 0, mn / cc, 0.0)
            rm = np.where(ll > 0, mn / ll, 0.0)
        p_sum += float((pm.sum(axis=1) / (cc > 0).sum(axis=1)).sum())
        r_sum += float((rm.sum(axis=1) / (ll > 0).sum(axis=1)).sum())
    p, r = p_sum / n, r_sum / n
    return p, r, (2 * p * r / (p + r) if p + r > 0 else 0.0)


def _indicator(items: list[str], sets_of: dict[str, set]) -> np.ndarray:
    cols = {}
    for i in items:
        for s in sorted(sets_of[i]):
            cols.setdefault(s, len(cols))
    out = np.zeros((len(items), len(cols)), dtype=np.float32)
    for r, i in enumerate(items):
        for s in sets_of[i]:
            out[r, cols[s]] = 1.0
    return out


def _reference_metrics(chains, clones, level_of, parent_of, reference, node_ids, space_of_node) -> dict[str, float]:
    keys = ["ref_skill", "ref_f1", "ref_precision", "ref_recall", "ref_best_level", "ref_f1_L1", "ref_f1_L2",
            "ref_f1_L3", "ref_baseline_f1", "ref_f1_one_group", "ref_f1_singletons", "ref_f1_compartments",
            "ref_labelled_share", "ref_n_labels"]
    out = dict.fromkeys(keys, NAN)
    if reference is None or len(reference) == 0:
        return out
    labels = defaultdict(set)
    for nid, lab in zip(reference["node_id"], reference["label"]):
        if nid in node_ids and lab != "":
            labels[nid].add(lab)
    items = sorted(n for n in labels if n in chains)
    out["ref_labelled_share"] = len(items) / max(1, len(node_ids))
    out["ref_n_labels"] = float(len({x for i in items for x in labels[i]}))
    if len(items) < 2:
        return out
    lab = _indicator(items, labels)

    def up(gid, depth):  # a group as seen at a level: itself if not deeper, else its ancestor at that level
        while gid in level_of and level_of[gid] > depth and parent_of.get(gid):
            gid = parent_of[gid]
        return gid

    best = (-1.0, NAN, NAN, NAN)
    for depth in (1, 2, 3):
        cl = {}
        for i in items:
            ch = chains[i]
            home = ch[min(depth, len(ch)) - 1]
            cl[i] = {home} | {up(g, depth) for g in clones.get(i, [])}
        p, r, f = _bcubed(_indicator(items, cl), lab)
        out[f"ref_f1_L{depth}"] = f
        if f > best[0] + 1e-12:
            best = (f, p, r, float(depth))
    one = _bcubed(np.ones((len(items), 1), dtype=np.float32), lab)[2]
    single = _bcubed(np.eye(len(items), dtype=np.float32), lab)[2]
    place = max(_bcubed(_indicator(items, {i: {space_of_node.get(i, ("", ""))[:d]} for i in items}), lab)[2]
                for d in (1, 2))         # one box per compartment, without and with its boundary box
    base = max(one, single, place)
    out.update(ref_f1_one_group=one, ref_f1_singletons=single, ref_f1_compartments=place)
    out.update(ref_f1=best[0], ref_precision=best[1], ref_recall=best[2], ref_best_level=best[3],
               ref_baseline_f1=base, ref_skill=(best[0] - base) / (1 - base) if base < 1 else NAN)
    return out


# -------------------------------------------------------------------------------------- name support
def _carriers(nodes: pd.DataFrame, source: str, xref: str, tree: _Tree, id_strings):
    """Nodes of the model that carry the identifier behind a group name, or None if it cannot be checked."""
    x = xref.strip()
    if not x:
        return set()
    if source == "sbo":
        return {r.node_id for r in nodes.itertuples(index=False) if x in r.sbo_path.split(">")}
    if source == "ec":
        code = x.split(":", 1)[1] if ":" in x else x
        code = "ec:" + ".".join(p for p in code.split(".") if p and p != "-")
        return {r.node_id for r in nodes.itertuples(index=False) if code in _ec_terms(r.ec)}
    if source == "go":
        return {r.node_id for r in nodes.itertuples(index=False) if x in r.go.split("|")}
    if source == "compartment":
        cid = x.partition("@")[0]
        return {r.node_id for r in nodes.itertuples(index=False)
                if r.loc_space == cid or r.compartment == cid or tree.space_of(r.compartment) == cid}
    if source in ("xref", "notes"):
        low = x.lower()
        return {n for n, ids in id_strings().items()
                if any(i == low or (i.startswith(low) and i[len(low):len(low) + 1] in "._ ") for i in ids)}
    return None


def _id_strings(nodes: pd.DataFrame) -> dict[str, set]:
    """node id -> every identifier string the node carries (lower case): db:id of each cross-reference, EC
    prefixes, GO ids, SBO path, and for the notes key:value, the value, and key:code for hierarchical codes."""
    out = {}
    for r in nodes.itertuples(index=False):
        ids = {res.lower() for res in (i.partition("=")[2] for i in r.xrefs.split("|") if i)}
        ids |= {t.lower() for t in _ec_terms(r.ec)} | {g.lower() for g in r.go.split("|") if g}
        if r.sbo_path:
            ids |= {t.lower() for t in r.sbo_path.split(">")}
        if r.notes:
            try:
                obj = json.loads(r.notes)
            except ValueError:
                obj = {}
            for k, value in obj.items():
                if isinstance(value, str):
                    for piece in re.split(r"[;|]", value):
                        piece = piece.strip().lower()
                        ids.add(f"{k.strip().lower()}:{piece}")
                        ids.add(piece)
            ids |= {t[len("notes:"):] for t in (c.lower() for c in _note_code_terms(r.notes))}
        out[r.node_id] = ids
    return out


def _resolvable(source: str, xref: str, dbs: set[str], comp_go: dict[str, str]) -> bool:
    """Does the identifier behind a name resolve outside the model (ontology term or database record)?"""
    if source == "sbo":
        return bool(re.fullmatch(r"SBO:\d{7}", xref))
    if source == "go":
        return bool(re.fullmatch(r"GO:\d{7}", xref))
    if source == "ec":
        return bool(re.fullmatch(r"(ec-code:|ec:)?\d+(\.(\d+|-)){0,3}", xref))
    if source == "compartment":
        return bool(comp_go.get(xref.partition("@")[0]))
    if source == "xref":
        db = xref.split(":", 1)[0].lower() if ":" in xref else ""
        return db in dbs and db not in ("uri", "")
    return False  # notes, member, none: local to the model


# ----------------------------------------------------------------------------------------- main entry
def metrics(nodes, edges, comps, groups, meta, reference=None, rerun=None) -> dict[str, float]:
    """All biology and FAIR metrics of one grouping of one model.

    nodes, edges, comps   the model's nodes.tsv, edges.tsv, compartments.tsv
    groups, meta          the strategy's groups.tsv and groups_meta.tsv
    reference             optional DataFrame (node_id, label): held-out curated labels, several rows per node
                          allowed, unlabelled nodes simply absent
    rerun                 optional (groups, meta) DataFrames of a second run of the same strategy
    Returns metric name -> float; NaN where a metric does not apply. Keys read by a rule are in RULES.
    """
    out: dict[str, float] = {}
    node = {r.node_id: r for r in nodes.itertuples(index=False)}
    node_ids = set(node)
    chains = {n: c for n, c in _chains(groups).items() if n in node_ids}
    clones = {n: c for n, c in _clones(groups).items() if n in node_ids}
    mrow = {r["group_id"]: r for r in meta.to_dict("records")}
    level_of = {g: int(r["level"]) for g, r in mrow.items() if str(r.get("level", "")).isdigit()}
    parent_of = {g: r.get("parent", "") for g, r in mrow.items()}
    leaf = {n: c[-1] for n, c in chains.items() if c}
    top = {n: c[0] for n, c in chains.items() if c}
    members = defaultdict(set)            # group id -> all nodes under it
    for n, c in chains.items():
        for g in c:
            members[g].add(n)
    tree = _Tree(comps)
    n_nodes = len(node_ids)
    out["n_nodes"] = float(n_nodes)
    out["n_groups"] = float(len(members))
    out["n_leaf_groups"] = float(len(set(leaf.values())))
    out["coverage"] = len([n for n in node_ids if n in leaf]) / max(1, n_nodes)

    # ---------------------------------------------------------------- BIO-1 agreement with the reference
    out.update(_reference_metrics(chains, clones, level_of, parent_of, reference, node_ids,
                                  {n: (r.loc_space, r.loc_kind == "boundary") for n, r in node.items()}))

    # ---------------------------------------------------------------- BIO-2 homogeneity
    sim = _Similarity(nodes)
    sig_of, m, k = sim.matrices()
    out["hom"] = _homogeneity(sig_of, m, k, leaf)
    out["hom_L1"] = _homogeneity(sig_of, m, k, top)
    for ch in CHANNELS:
        s1, m1, k1 = sim.matrices((ch,))
        out[f"hom_{ch}"] = _homogeneity(s1, m1, k1, leaf)
    out["hom_compartments"] = _homogeneity(sig_of, m, k, {n: r.loc_space for n, r in node.items()})
    out["hom_gain"] = (out["hom"] - max(0.0, out["hom_compartments"])
                       if not (_isnan(out["hom"]) or _isnan(out["hom_compartments"])) else NAN)
    out["annotated_share"] = len(sim.annotated) / max(1, n_nodes)

    # ---------------------------------------------------------------- LOC localisation fidelity
    def loc_group(g):  # (space, is boundary group) of a localisation group, else None
        r = mrow.get(g)
        if r is None or r.get("basis") != "localisation":
            return None
        cid, _, edge = r.get("compartment_id", "").partition("@")
        return (cid, edge == "boundary") if cid else None

    spaces = sorted({r.loc_space for r in node.values() if r.loc_space})
    loc_keys = ["loc_expressed", "loc_pure_leaf", "loc_repeat", "loc_boxes_mean", "loc_boxes_max", "loc_boundary",
                "loc_boundary_in_interior", "loc_nesting", "loc_span_ok", "loc_span_share", "loc_span_connected"]
    if len(spaces) < 2:
        out.update(dict.fromkeys(loc_keys, NAN))
    else:
        placed, outer_group, inner_group = {}, {}, {}
        for n, c in chains.items():
            hits = [g for g in c if (loc_group(g) or ("", False))[0] == node[n].loc_space]
            placed[n] = bool(hits)
            if hits:
                outer_group[n], inner_group[n] = hits[0], hits[-1]
        out["loc_expressed"] = sum(placed.values()) / max(1, len(chains))
        space_sets = {g: {node[x].loc_space for x in mem} for g, mem in members.items()}
        out["loc_pure_leaf"] = sum(len(space_sets[leaf[n]]) == 1 for n in leaf) / max(1, len(leaf))

        # LOC-3 transport and membrane participants on the boundary they cross
        bnodes = [n for n in chains if node[n].loc_kind == "boundary"]
        if bnodes:
            on = sum(any(loc_group(g) == (node[n].loc_space, True) for g in chains[n]) for n in bnodes)
            inside = sum(any(loc_group(g) == (node[n].loc_space, False) for g in chains[n])
                         and not any(loc_group(g) == (node[n].loc_space, True) for g in chains[n]) for n in bnodes)
            out["loc_boundary"] = on / len(bnodes)
            out["loc_boundary_in_interior"] = inside / len(bnodes)
        else:
            out["loc_boundary"] = out["loc_boundary_in_interior"] = NAN

        # LOC-4 nesting: every populated enclosing space is drawn exactly once around the node, in its scope
        populated = set(spaces)
        scope_boxes = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # parent group -> space -> [interior, boundary]
        for g in members:
            lg = loc_group(g)
            if lg:
                scope_boxes[parent_of.get(g, "")][lg[0]][int(lg[1])] += 1
        total = drawn = 0
        for n in chains:
            anc = [a for a in tree.ancestors(node[n].loc_space) if a in populated]
            total += len(anc)
            if anc and n in outer_group:
                here = scope_boxes[parent_of.get(outer_group[n], "")]
                for a in anc:
                    interior, boundary = here[a] if a in here else (0, 0)
                    drawn += (interior == 1) or (interior == 0 and boundary == 1)
        out["loc_nesting"] = drawn / total if total else NAN

        # LOC-2 how often the same compartment (or boundary) is drawn: effective number of boxes
        box_size = Counter(inner_group.values())
        key_of = {g: mrow[g]["compartment_id"] for g in box_size}
        key_size, key_boxes = Counter(), Counter()
        for g, size in box_size.items():
            key_size[key_of[g]] += size
            key_boxes[key_of[g]] += 1
        if inner_group:
            together = sum(box_size[g] / key_size[key_of[g]] for g in inner_group.values()) / len(inner_group)
            out["loc_repeat"] = 1.0 / together
            out["loc_boxes_mean"] = sum(key_boxes.values()) / len(key_boxes)
            out["loc_boxes_max"] = float(max(key_boxes.values()))
        else:
            out["loc_repeat"] = out["loc_boxes_mean"] = out["loc_boxes_max"] = NAN

        # LOC-5 a group that spans compartments holds the transport that links them and is split by place
        def connected(g) -> bool:
            par = {s: s for s in space_sets[g]}

            def find(a):
                par.setdefault(a, a)
                while par[a] != a:
                    par[a] = par[par[a]]
                    a = par[a]
                return a

            for x in members[g]:
                r = node[x]
                if r.kind != "reaction":
                    continue
                touched = sorted({tree.space_of(c) for c in r.compartments.split("|") if c} | {r.loc_space})
                for a in touched[1:]:
                    par[find(a)] = find(touched[0])
            return len({find(s) for s in space_sets[g]}) == 1

        spanning = {g for g in members if len(space_sets[g]) > 1}
        conn = {g: connected(g) for g in sorted(spanning)}
        in_span = [n for n in chains if any(g in spanning for g in chains[n])]
        bad = sum(1 for n in in_span if not placed[n] or not all(conn[g] for g in chains[n] if g in spanning))
        out["loc_span_ok"] = 1.0 - bad / max(1, len(chains))
        out["loc_span_share"] = len(in_span) / max(1, len(chains))
        out["loc_span_connected"] = (sum(all(conn[g] for g in chains[n] if g in spanning) for n in in_span)
                                     / len(in_span)) if in_span else NAN

    # ---------------------------------------------------------------- FAIR-1 determinism
    if rerun is None:
        out["fair_deterministic"] = out["fair_rerun_same_rows"] = NAN
    else:
        g2, m2 = rerun
        same = (list(groups.columns) == list(g2.columns) and list(meta.columns) == list(m2.columns)
                and groups.reset_index(drop=True).equals(g2.reset_index(drop=True))
                and meta.reset_index(drop=True).equals(m2.reset_index(drop=True)))
        out["fair_deterministic"] = 1.0 if same else 0.0
        a = {r[0]: tuple(r[1:]) for r in groups[["node_id", *LEVELS, "rule"]].itertuples(index=False)}
        b = {r[0]: tuple(r[1:]) for r in g2[["node_id", *LEVELS, "rule"]].itertuples(index=False)}
        out["fair_rerun_same_rows"] = sum(a[n] == b.get(n) for n in a) / max(1, len(a))
    out["fair_sorted"] = float(list(groups["node_id"]) == sorted(groups["node_id"])
                               and list(meta["group_id"]) == sorted(meta["group_id"]))
    rule = dict(zip(groups["node_id"], groups["rule"]))
    out["fair_rule_share"] = sum(1 for n in chains if rule.get(n, "").strip()) / max(1, len(chains))
    out["fair_n_rules"] = float(len({r for r in rule.values() if r.strip()}))

    # ---------------------------------------------------------------- FAIR-2 / FAIR-3 names
    dbs = set()
    for x in nodes["xrefs"]:
        for item in x.split("|"):
            res = item.partition("=")[2]
            if ":" in res:
                dbs.add(res.split(":", 1)[0].lower())
    comp_go = {r["compartment_id"]: r.get("go_id", "") for r in comps.to_dict("records")}
    kind = {n: r.kind for n, r in node.items()}
    eligible_total = carrier_total = lift_num = lift_den = 0.0
    named = traceable = vacuous = 0
    good_name, local_name, any_name = set(), set(), set()
    cache, lazy = {}, {}

    def id_strings():
        if "ids" not in lazy:
            lazy["ids"] = _id_strings(nodes)
        return lazy["ids"]

    for g in sorted(members):
        r = mrow.get(g)
        if r is None or r.get("name_source", "none") in ("none", ""):
            continue
        any_name.add(g)
        if len(members[g]) < 2:          # a name on a single node says nothing about a group
            continue
        named += 1
        source, xref = r["name_source"], r.get("name_xref", "")
        if source == "member":
            ok = xref in members[g]
            traceable += ok
            if ok:
                local_name.add(g)
            continue
        key = (source, xref)
        if key not in cache:
            cache[key] = _carriers(nodes, source, xref, tree, id_strings)
        car = cache[key]
        kinds = {kind[c] for c in car} if car else set()
        elig = [x for x in members[g] if kind[x] in kinds]
        hit = sum(x in car for x in elig) if car else 0
        if not hit:                      # the identifier is on no member: the name does not trace to the group
            pool = elig or sorted(members[g])
            prev = len(car) / sum(1 for x in node_ids if kind[x] in kinds) if car else 0.0
            eligible_total += len(pool)
            lift_num += len(pool) * (0.0 - prev) / (1 - prev) if prev < 1 else 0.0
            lift_den += len(pool)
            continue
        traceable += 1
        prev = len(car) / sum(1 for x in node_ids if kind[x] in kinds)
        if prev >= 1.0:                  # carried by every such node of the model: true, and says nothing
            vacuous += 1
            continue
        (good_name if _resolvable(source, xref, dbs, comp_go) else local_name).add(g)
        support = hit / len(elig)
        eligible_total += len(elig)
        carrier_total += hit
        lift_num += len(elig) * (support - prev) / (1 - prev)
        lift_den += len(elig)
    out["name_support"] = carrier_total / eligible_total if eligible_total else NAN
    out["name_lift"] = lift_num / lift_den if lift_den else NAN
    out["name_traceable"] = traceable / named if named else NAN
    out["name_vacuous_share"] = vacuous / named if named else NAN
    n_chain = max(1, len(chains))
    out["name_coverage"] = sum(any(g in any_name for g in c) for c in chains.values()) / n_chain
    out["name_resolvable"] = sum(any(g in good_name for g in c) for c in chains.values()) / n_chain
    only_local = sum(any(g in local_name for g in c) and not any(g in good_name for g in c) for c in chains.values())
    out["name_grounded"] = out["name_resolvable"] + 0.5 * only_local / n_chain
    out["named_group_share"] = len(any_name) / max(1, len(members))

    # ---------------------------------------------------------------- FAIR-4 nodes without annotation
    nbrs = defaultdict(set)
    for rid, sid in zip(edges["reaction_id"], edges["species_id"]):
        if rid in node_ids and sid in node_ids:
            nbrs[rid].add(sid)
            nbrs[sid].add(rid)
    annotated = set(sim.annotated)
    chain_set = {n: set(c) for n, c in chains.items()}

    def together(a, b) -> bool:
        return (leaf.get(a) is not None and leaf.get(a) == leaf.get(b)
                or any(g in chain_set.get(b, ()) for g in clones.get(a, ()))
                or any(g in chain_set.get(a, ()) for g in clones.get(b, ())))

    def stranded(pool):
        pool = [n for n in pool if nbrs.get(n) and n in chains]
        if not pool:
            return NAN
        return sum(not any(together(n, x) for x in nbrs[n]) for n in pool) / len(pool)

    poor = sorted(node_ids - annotated)
    out["poor_share"] = len(poor) / max(1, n_nodes)
    out["stranded_poor"] = stranded(poor)
    out["stranded_annotated"] = stranded(sorted(annotated))
    return out


# --------------------------------------------------------------------------------- verdicts and score
def judge(m: dict[str, float]) -> dict[str, str]:
    """Rule id -> pass | warn | fail | na, from the thresholds in RULES."""
    out = {}
    for r in RULES:
        v = m.get(r["metric"], NAN)
        if v is None or _isnan(v):
            out[r["id"]] = "na"
        elif r["direction"] == "higher":
            out[r["id"]] = "pass" if v >= r["passes"] - 1e-12 else "warn" if v >= r["warns"] - 1e-12 else "fail"
        else:
            out[r["id"]] = "pass" if v <= r["passes"] + 1e-12 else "warn" if v <= r["warns"] + 1e-12 else "fail"
    return out


def unit(rule: dict, v: float) -> float:
    """Utility of a metric value on 0..1: 0 at the rule's worst value, UTILITY_AT_WARN at the warn threshold,
    UTILITY_AT_PASS at the pass threshold, 1 at its best value, linear in between (NaN stays NaN).
    loc_repeat (1 .. infinity, lower is better) is first turned into 1 / loc_repeat."""
    if v is None or _isnan(v):
        return NAN
    worst, warn, ok, best = rule["worst"], rule["warns"], rule["passes"], rule["best"]
    if rule["metric"] == "loc_repeat":
        v, worst, warn, ok, best = 1.0 / max(v, 1.0), 0.0, 1.0 / warn, 1.0 / ok, 1.0
    elif rule["direction"] == "lower":
        v, worst, warn, ok, best = -v, -worst, -warn, -ok, -best
    xs, ys = [worst], [0.0]
    for x, y in ((warn, UTILITY_AT_WARN), (ok, UTILITY_AT_PASS), (best, 1.0)):
        if x > xs[-1] + 1e-12:
            xs.append(x)
            ys.append(y)
        else:                            # a threshold that coincides with the previous knot takes its utility
            ys[-1] = max(ys[-1], y)
    return float(np.interp(min(max(v, worst), best), xs, ys))


def score(m: dict[str, float]) -> dict[str, float]:
    """One number per grouping and model: weighted mean of the blocks that apply, 0 if not deterministic.

    Inside a block, the rules that apply are averaged with their weights; a block with no applicable rule
    is left out and the block weights are renormalised. name_lift counts as 0 when it is undefined because
    no group carries an informative name (nothing named, or only names that hold for every node).
    """
    blocks = {}
    for b in BLOCK_WEIGHTS:
        num = den = 0.0
        for r in RULES:
            if r["block"] != b:
                continue
            v = m.get(r["metric"], NAN)
            if r["metric"] == "name_lift" and (v is None or _isnan(v)):
                v = 0.0
            u = unit(r, v)
            if not _isnan(u):
                num += r["weight"] * u
                den += r["weight"]
        blocks[b] = num / den if den else NAN
    num = sum(BLOCK_WEIGHTS[b] * v for b, v in blocks.items() if not _isnan(v))
    den = sum(BLOCK_WEIGHTS[b] for b, v in blocks.items() if not _isnan(v))
    total = num / den if den else NAN
    if m.get("fair_deterministic", NAN) == 0.0:
        total = 0.0
    out = {"score": total}
    out.update({f"block_{b}": v for b, v in blocks.items()})
    return out


# ------------------------------------------------------------------------------------------------ CLI
def read_tsv(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def same_bytes(dir_a, dir_b) -> bool:
    """Rule FAIR-1 as the contract states it: groups.tsv and groups_meta.tsv byte-identical in two runs."""
    def digest(p):
        return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
    return all(digest(pathlib.Path(dir_a) / f) == digest(pathlib.Path(dir_b) / f)
               for f in ("groups.tsv", "groups_meta.tsv"))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Biology and FAIR metrics of one grouping of one SBML model.")
    ap.add_argument("features_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--reference", help="TSV with columns node_id, label")
    ap.add_argument("--rerun", help="output folder of a second run of the same strategy")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    f, o = pathlib.Path(a.features_dir), pathlib.Path(a.out_dir)
    ref = read_tsv(a.reference) if a.reference else None
    rerun = None
    if a.rerun:
        rerun = (read_tsv(pathlib.Path(a.rerun) / "groups.tsv"), read_tsv(pathlib.Path(a.rerun) / "groups_meta.tsv"))
    m = metrics(read_tsv(f / "nodes.tsv"), read_tsv(f / "edges.tsv"), read_tsv(f / "compartments.tsv"),
                read_tsv(o / "groups.tsv"), read_tsv(o / "groups_meta.tsv"), reference=ref, rerun=rerun)
    if a.rerun:
        m["fair_deterministic"] = float(same_bytes(o, a.rerun))
    verdict, total = judge(m), score(m)
    if a.json:
        print(json.dumps({"metrics": m, "verdict": verdict, "score": total}, indent=1, sort_keys=True))
        return
    for r in RULES:
        v = m.get(r["metric"], NAN)
        print(f"{r['id']:7s} {r['metric']:20s} {v:8.3f}  {verdict[r['id']]}")
    print("score   " + "  ".join(f"{k}={v:.3f}" for k, v in total.items()))


if __name__ == "__main__":
    main()
