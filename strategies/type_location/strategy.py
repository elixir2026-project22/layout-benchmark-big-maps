"""Run (one model; the CLI every strategy shares):
    python strategy.py --features panel/features/MODEL --sbml panel/models/MODEL.xml --out strategies/type_location/out/MODEL

Run (whole panel, PowerShell, from the folder that holds panel/, common/ and strategies/type_location/):
    foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") {
      python strategies/type_location/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/type_location/out/$m
      python common/contract.py validate panel/features/$m strategies/type_location/out/$m
      python common/contract.py write-sbml panel/models/$m.xml strategies/type_location/out/$m --features panel/features/$m --strategy type_location -o strategies/type_location/out/$m/${m}_grouped.xml
    }

Run (whole panel, bash):
    for m in brenda_merged e_coli_core R-HSA-77289 skm_pss skm_ja skm_sa; do
      python strategies/type_location/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/type_location/out/$m
      python common/contract.py validate panel/features/$m strategies/type_location/out/$m
      python common/contract.py write-sbml panel/models/$m.xml strategies/type_location/out/$m --features panel/features/$m --strategy type_location -o strategies/type_location/out/$m/${m}_grouped.xml
    done

Optional parameters (the defaults are the documented values, the same for every model):
    --theta 0.3333   --n-min 3   --k-max 11

type_location: top level (L1) by reaction type - the reaction's SBO term at a depth chosen by rule, refined by
EC class where SBO does not divide - and second level (L2) by cell localisation, one group per space inside
each type and an "@boundary" group for transport and membrane steps. Rules R1..R16 are in RULES.md; the
comments below carry the same numbers.

Reads nodes.tsv, edges.tsv, compartments.tsv of the features folder; writes groups.tsv and groups_meta.tsv.
The --sbml file is accepted for the common CLI but not read: everything the rules need is in the feature
tables. Standard library only, no random numbers, every iteration runs over a sorted sequence.
"""
from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import re

THETA = 1.0 / 3.0   # a divisible type holding more than this share of the reactions is split (R9)
N_MIN = 3           # a subtype with fewer reactions is not a group of its own (R10)
K_MAX = 11          # a split that would give more top-level types than this is not carried out (R9)

GROUP_COLUMNS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
META_COLUMNS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]

# SBO terms the fallback rules R2-R7 may assign to a reaction that carries no SBO term.
# (id, name) pairs along the SBO is_a path; checked against common/SBO_OBO.obo when this file was written.
_P = [("SBO:0000000", "systems biology representation"), ("SBO:0000231", "occurring entity representation"),
      ("SBO:0000375", "process")]
_BT = _P + [("SBO:0000167", "biochemical or transport reaction")]
_BR = _BT + [("SBO:0000176", "biochemical reaction")]
FALLBACK_PATHS = {
    "process": _P,
    "pseudoreaction": _P + [("SBO:0000631", "pseudoreaction")],
    "translocation": _BT + [("SBO:0000185", "translocation reaction")],
    "biochemical": _BR,
    "binding": _BR + [("SBO:0000177", "non-covalent binding")],
    "dissociation": _BR + [("SBO:0000180", "dissociation")],
}
SBO_COMPLEX = "SBO:0000253"  # non-covalent complex; protein complex and macromolecular complex are below it

# Enzyme Commission nomenclature (IUBMB): names of the seven classes and of their subclasses, keyed by the
# EC identifier. A constant of the nomenclature, like the SBO names; a subclass missing here is labelled "EC x.y".
EC_NAMES = {
    "1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases", "4": "Lyases", "5": "Isomerases",
    "6": "Ligases", "7": "Translocases",
    "1.1": "Acting on the CH-OH group of donors", "1.2": "Acting on the aldehyde or oxo group of donors",
    "1.3": "Acting on the CH-CH group of donors", "1.4": "Acting on the CH-NH2 group of donors",
    "1.5": "Acting on the CH-NH group of donors", "1.6": "Acting on NADH or NADPH",
    "1.7": "Acting on other nitrogenous compounds as donors", "1.8": "Acting on a sulfur group of donors",
    "1.9": "Acting on a heme group of donors", "1.10": "Acting on diphenols and related substances as donors",
    "1.11": "Acting on a peroxide as acceptor", "1.12": "Acting on hydrogen as donor",
    "1.13": "Acting on single donors with incorporation of molecular oxygen",
    "1.14": "Acting on paired donors, with incorporation or reduction of molecular oxygen",
    "1.15": "Acting on superoxide as acceptor", "1.16": "Oxidizing metal ions",
    "1.17": "Acting on CH or CH2 groups", "1.18": "Acting on iron-sulfur proteins as donors",
    "1.19": "Acting on reduced flavodoxin as donor", "1.20": "Acting on phosphorus or arsenic in donors",
    "1.21": "Catalysing the reaction X-H + Y-H = X-Y", "1.22": "Acting on halogen in donors",
    "1.23": "Reducing C-O-C group as acceptor", "1.97": "Other oxidoreductases",
    "1.98": "Enzymes using H2 as reductant", "1.99": "Other enzymes using O2 as oxidant",
    "2.1": "Transferring one-carbon groups", "2.2": "Transferring aldehyde or ketonic groups",
    "2.3": "Acyltransferases", "2.4": "Glycosyltransferases",
    "2.5": "Transferring alkyl or aryl groups, other than methyl groups",
    "2.6": "Transferring nitrogenous groups", "2.7": "Transferring phosphorus-containing groups",
    "2.8": "Transferring sulfur-containing groups", "2.9": "Transferring selenium-containing groups",
    "2.10": "Transferring molybdenum- or tungsten-containing groups",
    "3.1": "Acting on ester bonds", "3.2": "Glycosylases", "3.3": "Acting on ether bonds",
    "3.4": "Acting on peptide bonds (peptidases)",
    "3.5": "Acting on carbon-nitrogen bonds, other than peptide bonds", "3.6": "Acting on acid anhydrides",
    "3.7": "Acting on carbon-carbon bonds", "3.8": "Acting on halide bonds",
    "3.9": "Acting on phosphorus-nitrogen bonds", "3.10": "Acting on sulfur-nitrogen bonds",
    "3.11": "Acting on carbon-phosphorus bonds", "3.12": "Acting on sulfur-sulfur bonds",
    "3.13": "Acting on carbon-sulfur bonds",
    "4.1": "Carbon-carbon lyases", "4.2": "Carbon-oxygen lyases", "4.3": "Carbon-nitrogen lyases",
    "4.4": "Carbon-sulfur lyases", "4.5": "Carbon-halide lyases", "4.6": "Phosphorus-oxygen lyases",
    "4.7": "Carbon-phosphorus lyases", "4.8": "Nitrogen-oxygen lyases",
    "4.98": "ATP-independent chelatases", "4.99": "Other lyases",
    "5.1": "Racemases and epimerases", "5.2": "cis-trans-Isomerases", "5.3": "Intramolecular oxidoreductases",
    "5.4": "Intramolecular transferases", "5.5": "Intramolecular lyases",
    "5.6": "Isomerases altering macromolecular conformation", "5.99": "Other isomerases",
    "6.1": "Forming carbon-oxygen bonds", "6.2": "Forming carbon-sulfur bonds",
    "6.3": "Forming carbon-nitrogen bonds", "6.4": "Forming carbon-carbon bonds",
    "6.5": "Forming phosphoric-ester bonds", "6.6": "Forming nitrogen-metal bonds",
    "6.7": "Forming nitrogen-nitrogen bonds",
    "7.1": "Catalysing the translocation of hydrons", "7.2": "Catalysing the translocation of inorganic cations",
    "7.3": "Catalysing the translocation of inorganic anions and their chelates",
    "7.4": "Catalysing the translocation of amino acids and peptides",
    "7.5": "Catalysing the translocation of carbohydrates and their derivatives",
    "7.6": "Catalysing the translocation of other compounds",
}


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def safe(text):
    return re.sub(r"[^A-Za-z0-9]", "_", text)


def split(value):
    return [v for v in value.split("|") if v]


# ---------------------------------------------------------------------------------------------- type path
def ec_steps(ec_field):
    """R8: the EC class and subclass of a reaction as two further steps of its type path.

    Several EC numbers: the class most of them share, then the subclass most of those share; a tie goes to
    the lower number. Returns [] when the reaction has no EC number.
    """
    parsed = []
    for ec in split(ec_field):
        parts = ec.strip().split(".")
        if parts and parts[0].isdigit():
            parsed.append((int(parts[0]), int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else None))
    if not parsed:
        return []
    n_class = collections.Counter(c for c, _ in parsed)
    cls = min(n_class, key=lambda c: (-n_class[c], c))
    steps = [(f"EC:{cls}", f"EC {cls} {EC_NAMES.get(str(cls), '')}".strip(), "ec", f"ec-code:{cls}")]
    n_sub = collections.Counter(s for c, s in parsed if c == cls and s is not None)
    if n_sub:
        sub = min(n_sub, key=lambda s: (-n_sub[s], s))
        key = f"{cls}.{sub}"
        steps.append((f"EC:{key}", f"EC {key} {EC_NAMES.get(key, '')}".strip(), "ec", f"ec-code:{key}"))
    return steps


def type_path(rx, parts, species, space_of):
    """(steps, rule id): the reaction's type from general to specific. A step is (key, label, source, xref)."""
    if rx["sbo"]:                                                                       # R1
        ids = rx["sbo_path"].split(">")
        names = rx["sbo_path_names"].split(" > ")
        if len(names) != len(ids):
            names = ids
        pairs, rule = list(zip(ids, names)), "R1"
    else:
        reactants = [s for s, role in parts if role == "reactant"]
        products = [s for s, role in parts if role == "product"]
        modifiers = [s for s, role in parts if role == "modifier"]
        spaces = {space_of(c) for c in split(rx["compartments"])}

        def is_complex(s):
            return SBO_COMPLEX in species.get(s, {}).get("sbo_path", "").split(">")

        if reactants and products and len(spaces) >= 2:                                 # R2
            pairs, rule = FALLBACK_PATHS["translocation"], "R2"
        elif (not reactants or not products) and not modifiers:                         # R3
            pairs, rule = FALLBACK_PATHS["pseudoreaction"], "R3"
        elif len(set(reactants)) >= 2 and len(set(products)) == 1 and is_complex(products[0]):     # R4
            pairs, rule = FALLBACK_PATHS["binding"], "R4"
        elif len(set(reactants)) == 1 and len(set(products)) >= 2 and is_complex(reactants[0]):    # R5
            pairs, rule = FALLBACK_PATHS["dissociation"], "R5"
        elif reactants and products:                                                    # R6
            pairs, rule = FALLBACK_PATHS["biochemical"], "R6"
        else:                                                                           # R7
            pairs, rule = FALLBACK_PATHS["process"], "R7"
    steps = [(i, n, "sbo", i) for i, n in pairs]
    return steps + ec_steps(rx["ec"]), rule                                             # R8


# ---------------------------------------------------------------------------------------------- depth cut
class Type:
    """One top-level type: the reactions under `node` (a tuple of step keys).

    stage "sbo": the node may still be divided by its SBO subtypes; "ec": only by EC class (what is left of
    a term after its SBO subtypes were taken out); None: not divisible. `pooled` = reactions whose own
    subtype was too small to stand alone (R10).
    """

    def __init__(self, node, rx, stage="sbo", residual=False, pooled=()):
        self.node, self.rx, self.stage, self.residual, self.pooled = node, sorted(rx), stage, residual, set(pooled)


def children(node, rx, keys):
    """(direct reactions, {child node: reactions}) of a node of the type tree."""
    direct, kids = [], collections.defaultdict(list)
    for r in rx:
        if len(keys[r]) == len(node):
            direct.append(r)
        else:
            kids[keys[r][:len(node) + 1]].append(r)
    return direct, dict(kids)


def divide(node, rx, stage, keys, n_min):
    """One division step of a type: (subtypes that stand alone, the rest, pooled reactions, stage of the rest).

    SBO subtypes are used first (R9). Where a term has none of at least n_min reactions, or for what is left
    of a term after its SBO subtypes were taken out, the EC classes are used (R8). Subtypes with fewer than
    n_min reactions stay with the rest (R10).
    """
    direct, kids = children(node, rx, keys)
    for use in ("sbo", "ec"):
        if use == "sbo" and stage != "sbo":
            continue
        mine = {k: v for k, v in kids.items() if k[-1].startswith("EC:") == (use == "ec")}
        big = {k: v for k, v in mine.items() if len(v) >= n_min}
        if big:
            pooled = [r for k, v in kids.items() if k not in big and (use == "ec" or k in mine) for r in v]
            rest = sorted(direct + [r for k, v in kids.items() if k not in big for r in v])
            return big, rest, pooled, ("ec" if use == "sbo" else None)
    return {}, sorted(rx), [], None


def cut(keys, theta, n_min, k_max):
    """R9 + R10: choose the depth of the type tree. Returns the list of top-level types."""
    n_total = len(keys)
    types = [Type((), list(keys))]
    while True:
        # descend through steps that do not divide: no reaction stops here and all go one way
        moved = False
        for t in types:
            if t.residual:
                continue
            direct, kids = children(t.node, t.rx, keys)
            if not direct and len(kids) == 1:
                t.node, t.stage = next(iter(kids)), "sbo"
                moved = True
        if moved:
            continue
        # split the largest type that holds more than theta of the reactions and can be divided
        cand = sorted((t for t in types if t.stage and len(t.rx) > theta * n_total),
                      key=lambda t: (-len(t.rx), t.node, t.residual))
        if not cand:
            break
        t = cand[0]
        big, rest, pooled, rest_stage = divide(t.node, t.rx, t.stage, keys, n_min)
        new = [Type(k, v) for k, v in sorted(big.items())]
        if rest:
            new.append(Type(t.node, rest, stage=rest_stage, residual=True, pooled=set(pooled) | (t.pooled & set(rest))))
        if not big or len(types) - 1 + len(new) > k_max:
            t.stage = None
            continue
        types = [x for x in types if x is not t] + new
    return sorted(types, key=lambda t: (t.node, t.residual))


def common_prefix(paths):
    out = []
    for column in zip(*paths):
        if len(set(column)) != 1:
            break
        out.append(column[0])
    return tuple(out)


def name_type(node, rx, residual, keys, step):
    """(id stem, label, name_source, name_xref) of a group of the type tree.

    A plain group is named by its term. A pooled group (R10) is named by the most specific term all its
    reactions share; when that is the divided term itself, the label says "other <term>".
    """
    if residual:
        lca = common_prefix([keys[r] for r in rx])
        if len(lca) > len(node):
            node, residual = lca, False
    if not node:
        return "T_untyped", "", "none", ""
    _key, label, source, xref = step[node[-1]]
    stem = "T_" + safe(node[-1])
    return (stem + "_other", "other " + label, source, xref) if residual else (stem, label, source, xref)


# ---------------------------------------------------------------------------------------------- strategy
def run(features, out, theta=THETA, n_min=N_MIN, k_max=K_MAX):
    features, out = pathlib.Path(features), pathlib.Path(out)
    nodes = sorted(read_tsv(features / "nodes.tsv"), key=lambda n: n["node_id"])
    edges = read_tsv(features / "edges.tsv")
    comps = {c["compartment_id"]: c for c in read_tsv(features / "compartments.tsv")}
    node = {n["node_id"]: n for n in nodes}
    reactions = [n["node_id"] for n in nodes if n["kind"] == "reaction"]
    species = {n["node_id"]: n for n in nodes if n["kind"] == "species"}

    def space_of(cid):
        c = comps.get(cid)
        return c["encloses_id"] if c and c["kind"] == "membrane" and c["encloses_id"] else cid

    parts = collections.defaultdict(set)       # reaction -> {(species, role)}
    for e in edges:
        parts[e["reaction_id"]].add((e["species_id"], e["role"]))
    parts = {r: sorted(v) for r, v in parts.items()}

    # R1-R8: type path of every reaction
    keys, step, type_rule = {}, {}, {}
    for r in reactions:
        steps, rule = type_path(node[r], parts.get(r, []), species, space_of)
        keys[r] = tuple(s[0] for s in steps)
        type_rule[r] = rule
        for s in steps:
            step.setdefault(s[0], s)

    # R9, R10: top-level types
    meta, l1_of, placed_by = {}, {}, {}
    type_node = {}                              # L1 id -> its Type
    for t in (cut(keys, theta, n_min, k_max) if reactions else []):
        stem, label, source, xref = name_type(t.node, t.rx, t.residual, keys, step)
        meta[stem] = dict(group_id=stem, level=1, parent="", label=label, basis="reaction_type",
                          name_source=source, name_xref=xref, compartment_id="")
        type_node[stem] = t
        for r in t.rx:
            l1_of[r] = stem
            if source == "ec":
                placed_by[r] = "R8"             # the type that names its group is its EC class
            elif r in t.pooled and stem.endswith("_other"):
                placed_by[r] = "R10"            # its own subtype was too small to stand alone
            else:
                placed_by[r] = type_rule[r]     # R1 (own SBO term) or R2-R7 (fallback)

    # R11: home type of a species = the type holding most of its reactions
    sp_rx = collections.defaultdict(set)        # species -> reactions
    sp_made = collections.defaultdict(set)      # species -> reactions that list it as a product
    for r, ps in parts.items():
        for s, role in ps:
            sp_rx[s].add(r)
            if role == "product":
                sp_made[s].add(r)

    def cell(n):
        return (n["loc_space"], n["loc_kind"] == "boundary")

    sp_types = {}
    for s in sorted(species):
        score = collections.defaultdict(lambda: [0, 0, 0])
        for r in sorted(sp_rx.get(s, ())):
            sc = score[l1_of[r]]
            sc[0] += 1                                              # reactions of this type it takes part in
            sc[1] += node[r]["loc_space"] == species[s]["loc_space"]  # ... that lie in or on its own space
            sc[2] += r in sp_made[s]                                # ... that produce it
        sp_types[s] = score
        if score:
            l1_of[s] = min(score, key=lambda g: (-score[g][0], -score[g][1], -score[g][2], g))
            placed_by[s] = "R11"
        else:                                                       # R12: a species in no reaction
            l1_of[s] = "S_unconnected"
            placed_by[s] = "R12"
            meta.setdefault("S_unconnected", dict(group_id="S_unconnected", level=1, parent="", label="",
                                                  basis="mixed", name_source="none", name_xref="",
                                                  compartment_id=""))

    # second level
    cells = {cell(n) for n in nodes}
    l2_of = {}
    if len(cells) > 1:
        # R13, R14: one localisation group per space inside each type, boundary nodes in "@boundary" groups
        membrane_of = {}
        for cid in sorted(comps):
            c = comps[cid]
            if c["kind"] == "membrane" and c["encloses_id"]:
                membrane_of.setdefault(c["encloses_id"], cid)
        for n in nodes:
            nid = n["node_id"]
            space, on_boundary = cell(n)
            gid = f"{l1_of[nid]}__{safe(space)}" + ("__boundary" if on_boundary else "")
            if gid not in meta:
                name = comps.get(space, {}).get("name") or space
                if on_boundary and space in membrane_of:
                    label, xref = comps[membrane_of[space]]["name"] or membrane_of[space], membrane_of[space]
                elif on_boundary:
                    label, xref = name + " boundary", space
                else:
                    label, xref = name, space
                meta[gid] = dict(group_id=gid, level=2, parent=l1_of[nid], label=label, basis="localisation",
                                 name_source="compartment", name_xref=xref,
                                 compartment_id=space + ("@boundary" if on_boundary else ""))
            l2_of[nid] = gid
    else:
        # R15: one localisation cell in the whole model - the second level is the type one step deeper
        sub_of_rx = {}
        for stem in sorted(type_node):
            t = type_node[stem]
            big, rest, _pooled, _stage = divide(t.node, t.rx, t.stage or "ec", keys, n_min)
            if not big or (t.residual and t.stage is None):
                continue
            for k, v in sorted(big.items()) + ([(t.node, rest)] if rest else []):
                s_stem, label, source, xref = name_type(k, v, k == t.node, keys, step)
                gid = f"{stem}__{s_stem[2:]}"
                meta[gid] = dict(group_id=gid, level=2, parent=stem, label=label, basis="reaction_type",
                                 name_source=source, name_xref=xref, compartment_id="")
                for r in v:
                    sub_of_rx[r] = gid
        l2_of.update(sub_of_rx)
        for s in sorted(species):
            n_in = collections.Counter(sub_of_rx[r] for r in sp_rx.get(s, ())
                                       if l1_of[r] == l1_of[s] and r in sub_of_rx)
            if n_in:
                l2_of[s] = min(n_in, key=lambda g: (-n_in[g], g))

    # R16: a copy of the species in every other type it takes part in
    by_parent = collections.defaultdict(set)
    for nid, gid in l2_of.items():
        by_parent[l1_of[nid]].add(gid)
    clone_in = {}
    for s in sorted(species):
        targets = []
        for g in sorted(sp_types[s]):
            if g == l1_of[s]:
                continue
            if len(cells) > 1:
                space, on_boundary = cell(species[s])
                want = f"{g}__{safe(space)}" + ("__boundary" if on_boundary else "")
                targets.append(want if want in by_parent[g] else g)
            else:
                n_in = collections.Counter(l2_of[r] for r in sp_rx[s] if l1_of[r] == g and r in l2_of)
                targets.append(min(n_in, key=lambda x: (-n_in[x], x)) if n_in else g)
        clone_in[s] = "|".join(targets)

    rows = [dict(node_id=n["node_id"], L1=l1_of[n["node_id"]], L2=l2_of.get(n["node_id"], ""), L3="",
                 rule=placed_by[n["node_id"]], clone_in=clone_in.get(n["node_id"], "")) for n in nodes]
    used = {r["L1"] for r in rows} | {r["L2"] for r in rows if r["L2"]}
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", sorted(rows, key=lambda r: r["node_id"]), GROUP_COLUMNS)
    write_tsv(out / "groups_meta.tsv", [meta[g] for g in sorted(meta) if g in used], META_COLUMNS)
    return rows, [meta[g] for g in sorted(meta) if g in used]


def main(argv=None):
    ap = argparse.ArgumentParser(description="type_location grouping: reaction type (L1), cell localisation (L2)")
    ap.add_argument("--features", required=True, help="folder with nodes.tsv, edges.tsv, compartments.tsv")
    ap.add_argument("--sbml", required=True, help="the SBML model (accepted for the common CLI, not read)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    ap.add_argument("--theta", type=float, default=THETA)
    ap.add_argument("--n-min", type=int, default=N_MIN)
    ap.add_argument("--k-max", type=int, default=K_MAX)
    a = ap.parse_args(argv)
    rows, meta = run(a.features, a.out, a.theta, a.n_min, a.k_max)
    levels = collections.Counter(m["level"] for m in meta)
    print(f"{pathlib.Path(a.features).name}: {len(rows)} nodes -> {levels[1]} L1, {levels[2]} L2 groups -> {a.out}")


if __name__ == "__main__":
    main()
