"""The output contract of a grouping strategy: validate it, write it as SBML Level 3 groups.

Run:
    python common/contract.py validate  panel/features/MODEL  strategies/NAME/out/MODEL
    python common/contract.py write-sbml panel/models/MODEL.xml strategies/NAME/out/MODEL  -o MODEL_grouped.xml
    python common/contract.py strip-groups IN.xml -o OUT.xml

A strategy writes two files per model (tab-separated, UTF-8, "\n" line ends, rows sorted by id):

groups.tsv        one row per node of nodes.tsv - every species and reaction exactly once
    node_id
    L1, L2, L3    group id at each level. L1 is required. L2 and L3 are optional, but no gaps
                  (L3 needs L2). A group belongs to exactly one parent group.
    rule          id of the rule in the strategy's RULES.md that placed this node (e.g. "R3")
    clone_in      optional, species only: "|"-joined ids of further groups in which a copy of the
                  species is drawn. The node still has one home (L1..L3).

groups_meta.tsv   one row per group id used in groups.tsv
    group_id      unique over all levels
    level         1, 2 or 3
    parent        group id one level up; empty for level 1
    label         the group's name; empty when name_source is "none"
    basis         what the group is defined by: reaction_type | localisation | function | topology |
                  flow | mixed
    name_source   where the label comes from: sbo | ec | go | xref | notes | compartment | member | none
                  ("member": named after one member node, whose id is name_xref)
    name_xref     the identifier behind the label (SBO:0000176, ec-code:2.7.1, GO:0006096,
                  mapman:10.5, a compartment id, a node id); empty for "none"
    compartment_id  required when basis is localisation: a compartment_id of compartments.tsv, or
                  "<compartment_id>@boundary" for a group of nodes drawn on that space's boundary

Localisation is not a strategy's choice. nodes.tsv gives every node its place (loc_space, loc_kind);
a localisation group must contain only nodes whose loc_space is its compartment, and an "@boundary"
group only nodes with loc_kind "boundary". The drawing nests localisation groups by
compartments.tsv (organelles inside cytosol, extracellular around it), never side by side.
"""
from __future__ import annotations

import argparse
import collections
import csv
import pathlib
import re
import sys

GROUP_COLUMNS = ["node_id", "L1", "L2", "L3", "rule", "clone_in"]
META_COLUMNS = ["group_id", "level", "parent", "label", "basis", "name_source", "name_xref", "compartment_id"]
BASES = {"reaction_type", "localisation", "function", "topology", "flow", "mixed"}
NAME_SOURCES = {"sbo", "ec", "go", "xref", "notes", "compartment", "member", "none"}
LEVELS = ("L1", "L2", "L3")


def read_tsv(path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return [{k: (v or "") for k, v in row.items()} for row in csv.DictReader(fh, delimiter="\t")]


def load(features_dir, out_dir):
    f, o = pathlib.Path(features_dir), pathlib.Path(out_dir)
    return (read_tsv(f / "nodes.tsv"), read_tsv(f / "compartments.tsv"), read_tsv(o / "groups.tsv"),
            read_tsv(o / "groups_meta.tsv"))


def validate(features_dir, out_dir) -> list[str]:
    """Contract violations (empty list = valid)."""
    nodes, comps, groups, meta = load(features_dir, out_dir)
    err = []
    node = {n["node_id"]: n for n in nodes}
    for name, rows, cols in (("groups.tsv", groups, GROUP_COLUMNS), ("groups_meta.tsv", meta, META_COLUMNS)):
        missing = [c for c in cols if rows and c not in rows[0]]
        if not rows:
            err.append(f"{name}: empty")
        if missing:
            err.append(f"{name}: missing columns {missing}")
    if err:
        return err

    seen = collections.Counter(g["node_id"] for g in groups)
    dup = sorted(k for k, v in seen.items() if v > 1)
    if dup:
        err.append(f"{len(dup)} nodes appear more than once, e.g. {dup[:3]}")
    unknown = sorted(set(seen) - set(node))
    if unknown:
        err.append(f"{len(unknown)} ids are not in nodes.tsv, e.g. {unknown[:3]}")
    lost = sorted(set(node) - set(seen))
    if lost:
        err.append(f"{len(lost)} nodes of nodes.tsv are not grouped, e.g. {lost[:3]}")

    m = {}
    for row in meta:
        if row["group_id"] in m:
            err.append(f"group {row['group_id']}: more than one row in groups_meta.tsv")
        m[row["group_id"]] = row
    parent_seen = {}
    members = collections.defaultdict(list)
    for g in groups:
        if not g["L1"]:
            err.append(f"node {g['node_id']}: no L1")
        if g["L3"] and not g["L2"]:
            err.append(f"node {g['node_id']}: L3 without L2")
        if not g["rule"]:
            err.append(f"node {g['node_id']}: no rule")
        up = ""
        for i, lv in enumerate(LEVELS, start=1):
            gid = g[lv]
            if not gid:
                continue
            members[gid].append(g["node_id"])
            if gid not in m:
                err.append(f"group {gid}: used in groups.tsv, absent from groups_meta.tsv")
            else:
                if m[gid]["level"] != str(i):
                    err.append(f"group {gid}: used at level {i}, declared level {m[gid]['level']}")
                if m[gid]["parent"] != up:
                    err.append(f"group {gid}: parent in groups.tsv is '{up}', declared '{m[gid]['parent']}'")
            if parent_seen.setdefault(gid, up) != up:
                err.append(f"group {gid}: sits under more than one parent")
            up = gid
        for c in filter(None, g["clone_in"].split("|")):
            if node.get(g["node_id"], {}).get("kind") != "species":
                err.append(f"node {g['node_id']}: clone_in on a reaction")
            if c not in m:
                err.append(f"node {g['node_id']}: clone_in names unknown group {c}")
    for gid in sorted(set(m) - set(members)):
        err.append(f"group {gid}: declared but has no member")

    comp_ids = {c["compartment_id"] for c in comps}
    for gid, row in m.items():
        if row["basis"] not in BASES:
            err.append(f"group {gid}: basis '{row['basis']}' not in {sorted(BASES)}")
        if row["name_source"] not in NAME_SOURCES:
            err.append(f"group {gid}: name_source '{row['name_source']}' not in {sorted(NAME_SOURCES)}")
        if (row["name_source"] == "none") != (row["label"] == ""):
            err.append(f"group {gid}: label must be empty exactly when name_source is 'none'")
        if row["name_source"] not in ("none", "") and not row["name_xref"]:
            err.append(f"group {gid}: a named group needs name_xref")
        if row["basis"] == "localisation":
            cid, _, edge = row["compartment_id"].partition("@")
            if cid not in comp_ids or edge not in ("", "boundary"):
                err.append(f"group {gid}: compartment_id '{row['compartment_id']}' is not in compartments.tsv")
                continue
            wrong = [n for n in members.get(gid, []) if n in node and (
                node[n]["loc_space"] != cid or (edge and node[n]["loc_kind"] != "boundary"))]
            if wrong:
                err.append(f"group {gid}: {len(wrong)} members are not located at {row['compartment_id']}, "
                           f"e.g. {wrong[:3]}")
        elif row["compartment_id"]:
            err.append(f"group {gid}: compartment_id is only for basis 'localisation'")
    return sorted(set(err))


def _sid(text: str, used: set) -> str:
    base = "G_" + re.sub(r"[^A-Za-z0-9_]", "_", text)
    sid, i = base, 1
    while sid in used:
        i += 1
        sid = f"{base}_{i}"
    used.add(sid)
    return sid


def _miriam(xref: str) -> str:
    if re.fullmatch(r"(SBO|GO):\d{7}", xref):
        return f"https://identifiers.org/{xref}"
    m = re.fullmatch(r"ec-code:([\d.\-]+)", xref)
    if m:
        parts = (m.group(1).strip(".").split(".") + ["-"] * 4)[:4]
        return "https://identifiers.org/ec-code/" + ".".join(parts)
    return ""


def strip_groups(sbml_in, sbml_out, drop_note_keys=()) -> int:
    """Copy an SBML file without its groups (and without the named note keys). Returns groups removed."""
    import libsbml

    doc = libsbml.readSBMLFromFile(str(sbml_in))
    model = doc.getModel()
    plugin = model.getPlugin("groups")
    n = 0
    if plugin is not None:
        n = plugin.getNumGroups()
        for _ in range(n):
            plugin.removeGroup(0)
        doc.disablePackage(plugin.getURI(), "groups")
    if drop_note_keys:
        pat = re.compile(r"<p>\s*(" + "|".join(re.escape(k) for k in drop_note_keys) + r")\s*:.*?</p>\s*", re.S)
        for el in ([model.getSpecies(i) for i in range(model.getNumSpecies())]
                   + [model.getReaction(i) for i in range(model.getNumReactions())]):
            raw = el.getNotesString()
            if raw and pat.search(raw):
                el.unsetNotes()
                el.setNotes(pat.sub("", raw))
    libsbml.writeSBMLToFile(doc, str(sbml_out))
    return n


def write_sbml(sbml_in, features_dir, out_dir, sbml_out, strategy="", path_names=True) -> int:
    """The model with the grouping as nested SBML L3 groups (kind partonomy, SBO:0000633).

    Deepest groups are written first, so a reader that takes the first group containing a node
    (layout-benchmark's grouped layouts) gets the innermost group. A parent group lists its child
    groups and any nodes placed directly in it. path_names: a group's name is "L1 / L2 / L3" (readable
    when groups are drawn without their parents); otherwise its own label.
    """
    import libsbml

    errors = validate(features_dir, out_dir)
    if errors:
        raise SystemExit("not written, the grouping breaks the contract:\n  " + "\n  ".join(errors[:20]))
    _nodes, _comps, groups, meta = load(features_dir, out_dir)
    m = {r["group_id"]: r for r in meta}
    doc = libsbml.readSBMLFromFile(str(sbml_in))
    model = doc.getModel()
    if doc.getLevel() < 3:
        raise SystemExit("SBML groups need Level 3")
    plugin = model.getPlugin("groups")
    if plugin is None:
        doc.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), "groups", True)
        doc.setPackageRequired("groups", False)
        plugin = model.getPlugin("groups")
    for _ in range(plugin.getNumGroups()):
        plugin.removeGroup(0)

    direct, kids = collections.defaultdict(list), collections.defaultdict(set)
    for g in groups:
        chain = [g[lv] for lv in LEVELS if g[lv]]
        direct[chain[-1]].append(g["node_id"])
        for a, b in zip(chain, chain[1:]):
            kids[a].add(b)

    def label(gid):
        return m[gid]["label"] or f"[unnamed {gid}]"

    def path(gid):
        out = [gid]
        while m[out[-1]]["parent"]:
            out.append(m[out[-1]]["parent"])
        return " / ".join(label(x) for x in reversed(out))

    used = set()
    for i in range(model.getNumSpecies()):
        used.add(model.getSpecies(i).getId())
    for i in range(model.getNumReactions()):
        used.add(model.getReaction(i).getId())
    sid = {}
    order = sorted(m, key=lambda x: (-int(m[x]["level"]), x))
    for gid in order:
        sid[gid] = _sid(gid, used)
    for gid in order:
        row = m[gid]
        grp = plugin.createGroup()
        grp.setId(sid[gid])
        grp.setName(path(gid) if path_names else label(gid))
        grp.setKind("partonomy")
        grp.setSBOTerm("SBO:0000633")
        for ref in sorted(direct.get(gid, [])) + sorted(sid[k] for k in kids.get(gid, ())):
            grp.createMember().setIdRef(ref)
        note = (f"level: {row['level']}; basis: {row['basis']}; name_source: {row['name_source']}; "
                f"name_xref: {row['name_xref']}; compartment_id: {row['compartment_id']}; strategy: {strategy}")
        grp.setNotes(f'<body xmlns="http://www.w3.org/1999/xhtml"><p>{note}</p></body>')
        uri = _miriam(row["name_xref"])
        if uri:
            grp.setMetaId("meta_" + sid[gid])
            cv = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
            cv.setBiologicalQualifierType(libsbml.BQB_IS)
            cv.addResource(uri)
            grp.addCVTerm(cv)
    libsbml.writeSBMLToFile(doc, str(sbml_out))
    return len(order)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("features_dir")
    v.add_argument("out_dir")
    w = sub.add_parser("write-sbml")
    w.add_argument("sbml")
    w.add_argument("out_dir")
    w.add_argument("--features", required=True)
    w.add_argument("-o", "--out", required=True)
    w.add_argument("--strategy", default="")
    w.add_argument("--own-names", action="store_true", help="group name = its own label, not the path")
    s = sub.add_parser("strip-groups")
    s.add_argument("sbml")
    s.add_argument("-o", "--out", required=True)
    s.add_argument("--drop-note-key", action="append", default=[])
    a = ap.parse_args(argv)
    if a.cmd == "validate":
        errors = validate(a.features_dir, a.out_dir)
        print("\n".join(errors) if errors else "valid")
        sys.exit(1 if errors else 0)
    if a.cmd == "write-sbml":
        n = write_sbml(a.sbml, a.features, a.out_dir, a.out, a.strategy, not a.own_names)
        print(f"{n} groups -> {a.out}")
    if a.cmd == "strip-groups":
        print(f"{strip_groups(a.sbml, a.out, a.drop_note_key)} groups removed -> {a.out}")


if __name__ == "__main__":
    main()
