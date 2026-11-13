"""location_type: cell localisation on top, reaction type inside each compartment.

Run (from the bundle root, the folder that holds common/ and panel/):
    python strategies/location_type/strategy.py --features panel/features/skm_pss --sbml panel/models/skm_pss.xml --out strategies/location_type/out/skm_pss
    python common/contract.py validate panel/features/skm_pss strategies/location_type/out/skm_pss
    python common/contract.py write-sbml panel/models/skm_pss.xml strategies/location_type/out/skm_pss --features panel/features/skm_pss --strategy location_type -o strategies/location_type/out/skm_pss/skm_pss_grouped.xml

All six panel models in one go (PowerShell; in bash: for m in brenda_merged e_coli_core R-HSA-77289 skm_pss skm_ja skm_sa; do ...; done):
    foreach ($m in "brenda_merged","e_coli_core","R-HSA-77289","skm_pss","skm_ja","skm_sa") { python strategies/location_type/strategy.py --features panel/features/$m --sbml panel/models/$m.xml --out strategies/location_type/out/$m }

Variant that keeps the compartment as level 1 even when the model has one compartment only:
    python strategies/location_type/strategy.py --features panel/features/brenda_merged --sbml panel/models/brenda_merged.xml --out strategies/location_type/out_keep_single/brenda_merged --keep-single-compartment

What it does (the numbered rules R1-R12 are in RULES.md; the `rule` column of groups.tsv uses the same ids):
    level 1   one group per cell space (nodes.tsv loc_space), plus "<space>@boundary" for nodes drawn
              on that space's boundary (loc_kind = boundary)                                  R1, R2
    level 2   inside each level-1 group, reactions by type: their SBO process term, or their EC class
              where the SBO term does not divide the group (no SBO term, or one term for all)  R4, R5
    level 3   only inside a level-2 box of more than MAX_BOX nodes: the next key of the type ladder
              that divides it (SBO term > EC class > EC sub-class > EC sub-sub-class > SBO role of
              the modifiers)                                                                   R11, R6, R7
    species   home = the type box of its own compartment holding most of its reactions (R9, R10);
              copies (clone_in) where it is only a modifier, and for a species shared by three or
              more boxes (R12)
    no box    a reaction with no value for the key, a value held by fewer than MIN_RX reactions, or
              a box of fewer than MIN_BOX nodes: it stays directly in the enclosing group        R8
    one compartment only: the localisation level says nothing and is skipped; type starts at
              level 1 (R3). --keep-single-compartment switches R3 off.

Reads only nodes.tsv, edges.tsv and compartments.tsv of --features. --sbml is accepted because the
contract's command line has it; nothing is taken from the file. Standard library only, no random
numbers, every iteration over a sorted sequence: two runs give byte-identical files.
"""
from __future__ import annotations

import argparse
import collections
import csv
import pathlib

# ---- parameters (RULES.md, parameter table) --------------------------------------------------------
MAX_BOX = 50      # a box with more nodes than this gets one more level (if a key divides it)
MIN_BOX = 5       # fewest nodes (reactions + species at home) for a type box (= MAX_BOX / 10)
MIN_RX = 3        # fewest reactions for a type box
HUB_BOXES = 3     # a species taking part in this many boxes of its compartment is copied into each
SBO_PROCESS = "SBO:0000375"   # "process", the root of all reaction types: neither it nor its direct children name a box
SBO_MODIFIER = "SBO:0000019"  # "modifier", the root of all modifier roles: too general to name a box
# The seven top-level classes of the Enzyme Commission nomenclature (the only vocabulary in this file).
EC_CLASS = {"1": "oxidoreductases", "2": "transferases", "3": "hydrolases", "4": "lyases", "5": "isomerases",
            "6": "ligases", "7": "translocases"}

KEYS = ["sbo", "ec1", "ec2", "ec3", "role"]                       # the type ladder, in order
FIRST_KEYS = ["sbo", "ec1"]                                       # keys that may open the first type level
KEY_RULE = {"sbo": "R4", "ec1": "R5", "ec2": "R6", "ec3": "R6", "role": "R7"}   # rule id per ladder key
GROUP_COLUMNS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
META_COLUMNS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def _num(part):
    return (0, int(part)) if part.isdigit() else (1, part)


def ec_key(code):
    return tuple(_num(p) for p in code.split("."))


def ec_prefix(ec_field, depth, within=""):
    """Most frequent EC prefix of `depth` fields among a reaction's EC numbers (below prefix `within`).

    Ties go to the numerically smallest prefix. "" when the reaction has no EC number that deep."""
    count = collections.Counter()
    for code in filter(None, ec_field.split("|")):
        parts = [p for p in code.strip().split(".") if p not in ("", "-")]
        if len(parts) < depth or not all(p.isdigit() for p in parts[:depth]):
            continue
        pre = ".".join(parts[:depth])
        if within and not (pre == within or pre.startswith(within + ".")):
            continue
        count[pre] += 1
    if not count:
        return ""
    top = max(count.values())
    return sorted((p for p, n in count.items() if n == top), key=ec_key)[0]


class Model:
    def __init__(self, features):
        f = pathlib.Path(features)
        self.nodes = sorted(read_tsv(f / "nodes.tsv"), key=lambda n: n["node_id"])
        self.edges = read_tsv(f / "edges.tsv")
        self.comps = {c["compartment_id"]: c for c in read_tsv(f / "compartments.tsv")}
        self.node = {n["node_id"]: n for n in self.nodes}
        self.reactions = [n["node_id"] for n in self.nodes if n["kind"] == "reaction"]
        self.species = [n["node_id"] for n in self.nodes if n["kind"] != "reaction"]
        self.sbo_name, self.sbo_parent, self.sbo_depth = {}, {}, {}
        for n in self.nodes:
            ids = [x.strip() for x in n["sbo_path"].split(">") if x.strip()]
            names = [x.strip() for x in n["sbo_path_names"].split(">")]
            for i, term in enumerate(ids):
                self.sbo_depth[term] = i
                if i < len(names) and len(names) == len(ids):
                    self.sbo_name.setdefault(term, names[i])
                if i:
                    self.sbo_parent.setdefault(term, ids[i - 1])
            if n["sbo"]:
                self.sbo_name[n["sbo"]] = n["sbo_name"] or self.sbo_name.get(n["sbo"], "")
        self.role = {}            # reaction -> SBO role of its modifiers, generic "modifier" ignored (most frequent, ties: lowest id)
        roles = collections.defaultdict(collections.Counter)
        self.by_species = collections.defaultdict(list)   # species -> [(reaction, role)]
        for e in self.edges:
            if e["reaction_id"] in self.node and e["species_id"] in self.node:
                self.by_species[e["species_id"]].append((e["reaction_id"], e["role"]))
            if e["role"] == "modifier" and e["sbo"] and e["sbo"] != SBO_MODIFIER:
                roles[e["reaction_id"]][e["sbo"]] += 1
                self.sbo_name.setdefault(e["sbo"], e["sbo_name"])
        for r, c in roles.items():
            top = max(c.values())
            self.role[r] = sorted(t for t, n in c.items() if n == top)[0]

    # ---- rules R1, R2: where a node is drawn -------------------------------------------------------
    def loc_group(self, nid):
        n = self.node[nid]
        return n["loc_space"] + ("@boundary" if n["loc_kind"] == "boundary" else "")

    # ---- the type ladder ----------------------------------------------------------------------------
    def can_roll_up(self, term):
        """A rare SBO term may be replaced by its parent unless the parent is "process" or a direct
        child of it (e.g. "biochemical or transport reaction"): those say too little to name a box."""
        parent = self.sbo_parent.get(term)
        return bool(parent) and parent != SBO_PROCESS and self.sbo_parent.get(parent) not in (None, SBO_PROCESS)

    def split(self, rxns, key, ec_within="", banned=()):
        """Divide reactions by one ladder key. Returns ({value: [reactions]}, [left over]) or None.

        A key divides a group when at least two of its values hold MIN_RX reactions each. Reactions
        with no value, or with a value held by fewer than MIN_RX reactions, are left over."""
        value = {}
        if key == "sbo":
            value = {r: self.node[r]["sbo"] for r in rxns if self.node[r]["sbo"]}
            while True:  # roll rare terms up to their SBO parent, deepest first, never up to "process"
                count = collections.Counter(value.values())
                rare = sorted((t for t, n in count.items() if n < MIN_RX and self.can_roll_up(t)),
                              key=lambda t: (-self.sbo_depth.get(t, 0), t))
                if not rare:
                    break
                value = {r: (self.sbo_parent[v] if v == rare[0] else v) for r, v in value.items()}
        elif key in ("ec1", "ec2", "ec3"):
            depth = int(key[2])
            for r in rxns:
                v = ec_prefix(self.node[r]["ec"], depth, ec_within)
                if v:
                    value[r] = v
        elif key == "role":
            value = {r: self.role[r] for r in rxns if r in self.role}
        boxes = collections.defaultdict(list)
        for r in rxns:
            if r in value:
                boxes[value[r]].append(r)
        boxes = {v: rs for v, rs in boxes.items() if len(rs) >= MIN_RX and v not in banned}
        if len(boxes) < 2:
            return None
        placed = {r for rs in boxes.values() for r in rs}
        return boxes, [r for r in rxns if r not in placed]

    def describe(self, key, value):
        """(group id part, label, name_source, name_xref) of a type box."""
        if key == "sbo":
            return "sbo_" + value.split(":")[-1], self.sbo_name.get(value, value), "sbo", value
        if key == "role":
            return "role_" + value.split(":")[-1], "with " + self.sbo_name.get(value, value), "sbo", value
        label = "EC " + value + (" " + EC_CLASS[value] if value in EC_CLASS else "")
        return "ec_" + value, label, "ec", "ec-code:" + value


def build(features, keep_single=False):
    m = Model(features)
    meta = {}                                  # group id -> meta row
    path = {}                                  # node id -> [group ids]
    rule = {}                                  # node id -> rule id
    todo = {}                                  # group id -> (remaining keys, ec prefix): may get one more level

    def add_group(gid, level, parent, label, basis, source, xref, comp=""):
        meta[gid] = dict(group_id=gid, level=level, parent=parent, label=label, basis=basis, name_source=source,
                         name_xref=xref, compartment_id=comp)

    def loc_meta(lg, level, parent):
        cid, _, edge = lg.partition("@")
        name = m.comps.get(cid, {}).get("name") or cid
        if edge:   # a membrane compartment that encloses this space names the boundary, else "<space> boundary"
            membranes = sorted(c for c, row in m.comps.items() if row["kind"] == "membrane" and row["encloses_id"] == cid)
            if membranes:
                add_group("loc_" + lg, level, parent, m.comps[membranes[0]]["name"] or membranes[0], "localisation",
                          "compartment", membranes[0], lg)
            else:
                add_group("loc_" + lg, level, parent, name + " boundary", "localisation", "compartment", cid, lg)
        else:
            add_group("loc_" + lg, level, parent, name, "localisation", "compartment", cid, lg)

    def members(gid):
        return [n for n, p in path.items() if gid in p]

    def place_species(parent_gid, species, level):
        """R9, R10: each species goes to the child box (at `level`) of `parent_gid` holding most of its reactions.

        Only reactions inside `parent_gid` count (None = the whole model). Ties: most uses as product,
        then as reactant, then a box before "directly in the parent", then the lowest group id."""
        for s in species:
            tally = collections.defaultdict(lambda: [0, 0, 0])
            for r, role in m.by_species.get(s, []):
                p = path[r]
                if parent_gid is not None and (len(p) < level - 1 or p[level - 2] != parent_gid):
                    continue
                t = tally[p[level - 1] if len(p) >= level else ""]
                t[0] += 1
                t[1] += role == "product"
                t[2] += role == "reactant"
            if not tally:
                rule[s] = "R10"
                continue
            best = sorted(tally, key=lambda c: (-tally[c][0], -tally[c][1], -tally[c][2], c == "", c))[0]
            rule[s] = "R9"
            if best:
                path[s] = path[s][:level - 1] + [best]

    def subdivide(parent_gid, level, rxns, species, keys, ec_within="", first_type_level=False):
        """One more type level inside `parent_gid` (None = the whole model). True when a key divided it.

        Tries the ladder keys in order. For a key that divides the reactions (R4-R7), species are placed
        (R9, R10); a box left with fewer than MIN_BOX nodes is withdrawn and the key is tried again without
        that value (R8). A key that no longer gives two boxes is dropped and the next key is tried."""
        before = {n: (list(path[n]), rule[n]) for n in list(rxns) + list(species)}
        for key in keys:
            if key in ("ec2", "ec3") and not ec_within:
                continue  # an EC sub-class level only inside an EC box
            banned = set()
            while True:
                for n, (p, r) in before.items():
                    path[n], rule[n] = list(p), r
                got = m.split(rxns, key, ec_within, banned)
                if not got:
                    break
                boxes, left = got
                gid_of = {}
                for value in sorted(boxes, key=lambda v: ec_key(v) if key.startswith("ec") else (v,)):
                    gid_of[value] = (parent_gid + "/" if parent_gid else "") + m.describe(key, value)[0]
                    for r in boxes[value]:
                        path[r] = path[r][:level - 1] + [gid_of[value]]
                        rule[r] = KEY_RULE[key]
                if first_type_level:
                    for r in left:
                        rule[r] = "R8"   # at a deeper level a left-over reaction keeps the rule of the box it stays in
                place_species(parent_gid, species, level)
                size = collections.Counter(path[n][level - 1] for n in before if len(path[n]) >= level)
                small = {v for v, gid in gid_of.items() if size[gid] < MIN_BOX}
                if small:
                    banned |= small
                    continue
                for value, gid in gid_of.items():
                    _part, label, source, xref = m.describe(key, value)
                    add_group(gid, level, parent_gid or "", label, "reaction_type", source, xref)
                    todo[gid] = (KEYS[KEYS.index(key) + 1:], value if key.startswith("ec") else ec_within)
                return True
        for n, (p, r) in before.items():
            path[n], rule[n] = list(p), r
        return False

    loc_groups = sorted({m.loc_group(n["node_id"]) for n in m.nodes})
    keep_loc = keep_single or len(loc_groups) > 1
    if keep_loc:
        # R1, R2: level 1 = localisation; level 2 = type inside each localisation group (always tried)
        for lg in loc_groups:
            loc_meta(lg, 1, "")
        for n in m.nodes:
            path[n["node_id"]] = ["loc_" + m.loc_group(n["node_id"])]
            rule[n["node_id"]] = "R2" if n["loc_kind"] == "boundary" else "R1"
        for lg in loc_groups:
            gid = "loc_" + lg
            inside = members(gid)
            subdivide(gid, 2, [n for n in inside if m.node[n]["kind"] == "reaction"],
                      [n for n in inside if m.node[n]["kind"] != "reaction"],
                      KEYS if len(inside) > MAX_BOX else FIRST_KEYS, first_type_level=True)
        deeper = [3]
    else:
        # R3: one localisation group only - the level is skipped, type starts at level 1; whatever no
        # type box takes stays in the compartment's own group
        lg = loc_groups[0]
        for n in m.nodes:
            path[n["node_id"]] = []
            rule[n["node_id"]] = "R3"
        if subdivide(None, 1, m.reactions, m.species, KEYS if len(m.nodes) > MAX_BOX else FIRST_KEYS,
                     first_type_level=True):
            rest = [n["node_id"] for n in m.nodes if not path[n["node_id"]]]
            if rest:
                loc_meta(lg, 1, "")
            for n in rest:
                path[n] = ["loc_" + lg]
        else:
            loc_meta(lg, 1, "")
            for n in m.nodes:
                path[n["node_id"]] = ["loc_" + lg]
                rule[n["node_id"]] = "R1"
        deeper = [2, 3]

    # R11: one more level inside a type box that is still too large to read
    for level in deeper:
        for gid in sorted(g for g, row in meta.items() if row["level"] == level - 1 and g in todo):
            inside = members(gid)
            if len(inside) <= MAX_BOX:
                continue
            keys, ec_within = todo[gid]
            subdivide(gid, level, [n for n in inside if m.node[n]["kind"] == "reaction" and path[n][-1] == gid],
                      [n for n in inside if m.node[n]["kind"] != "reaction" and path[n][-1] == gid], keys, ec_within)

    # R12: copies of a species, only inside its own localisation group
    clones = {}
    for s in m.species:
        home = path[s]
        touched = collections.defaultdict(set)
        for r, role in m.by_species.get(s, []):
            if m.loc_group(r) != m.loc_group(s):
                continue
            touched[path[r][-1]].add(role)
        boxes = set(touched) | {home[-1]}
        out = []
        for box in sorted(touched):
            if box in home:
                continue  # its home, or a box that encloses its home
            if touched[box] == {"modifier"} or len(boxes) >= HUB_BOXES:
                out.append(box)
        clones[s] = "|".join(out)

    used = {g for p in path.values() for g in p}
    meta_rows = [meta[g] for g in sorted(meta) if g in used]
    group_rows = []
    for n in m.nodes:
        p = path[n["node_id"]] + ["", "", ""]
        group_rows.append(dict(node_id=n["node_id"], L1=p[0], L2=p[1], L3=p[2], rule=rule[n["node_id"]],
                               clone_in=clones.get(n["node_id"], "")))
    return group_rows, meta_rows


def main(argv=None):
    global MIN_RX, MIN_BOX, MAX_BOX, HUB_BOXES
    ap = argparse.ArgumentParser(description="location_type grouping: cell localisation, then reaction type")
    ap.add_argument("--features", required=True, help="folder with nodes.tsv, edges.tsv, compartments.tsv")
    ap.add_argument("--sbml", required=True, help="the SBML model (accepted for the contract; not read)")
    ap.add_argument("--out", required=True, help="output folder for groups.tsv and groups_meta.tsv")
    ap.add_argument("--keep-single-compartment", action="store_true",
                    help="variant: keep the compartment as level 1 even when it is the only one")
    ap.add_argument("--min-rx", type=int, default=MIN_RX, help="sensitivity runs only; default is the rule value")
    ap.add_argument("--min-box", type=int, default=MIN_BOX, help="sensitivity runs only; default is the rule value")
    ap.add_argument("--max-box", type=int, default=MAX_BOX, help="sensitivity runs only; default is the rule value")
    ap.add_argument("--hub-boxes", type=int, default=HUB_BOXES, help="sensitivity runs only; default is the rule value")
    a = ap.parse_args(argv)
    MIN_RX, MIN_BOX, MAX_BOX, HUB_BOXES = a.min_rx, a.min_box, a.max_box, a.hub_boxes
    groups, meta = build(a.features, a.keep_single_compartment)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "groups.tsv", GROUP_COLUMNS, groups)
    write_tsv(out / "groups_meta.tsv", META_COLUMNS, meta)
    levels = collections.Counter(row["level"] for row in meta)
    print(f"{out}: {len(groups)} nodes, groups per level {dict(sorted(levels.items()))}")


if __name__ == "__main__":
    main()
