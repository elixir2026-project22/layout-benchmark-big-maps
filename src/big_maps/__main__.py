"""The ``big-maps`` command line.

    big-maps prepare MODEL.xml [MODEL.xml ...]   put models in panel/models/ and write panel/features/
    big-maps prepare --panel                     rewrite panel/features/ for every model already there
    big-maps options                             the layout options, and which of them can run here
    big-maps models                              the models of the panel, from panel/panel.tsv

Everything else in this repository is run as a script from the repository root;
the README gives the commands.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import shutil
import sys

from big_maps import __version__, repo_root

PREPARED_COLUMNS = [
    "model", "source", "source_sha256", "groups_removed",
    "species", "reactions", "compartments", "prepared",
]


def _common(root: pathlib.Path):
    """Import the repository's own ``common/`` modules from *root*."""
    sys.path.insert(0, str(root / "common"))
    import sbml_features  # noqa: PLC0415
    import contract  # noqa: PLC0415
    return sbml_features, contract


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def features(root: pathlib.Path, model: str, sbml: pathlib.Path) -> dict:
    """Write ``panel/features/<model>/`` from *sbml*; return its counts."""
    F, _ = _common(root)
    nodes, edges, comps = F.extract(str(sbml))
    out = root / "panel" / "features" / model
    out.mkdir(parents=True, exist_ok=True)
    F.write_tsv(out / "nodes.tsv", nodes, F.NODE_COLUMNS)
    F.write_tsv(out / "edges.tsv", edges, F.EDGE_COLUMNS)
    F.write_tsv(out / "compartments.tsv", comps, F.COMPARTMENT_COLUMNS)
    n_species = sum(1 for r in nodes if r["kind"] == "species")
    return {"species": n_species, "reactions": len(nodes) - n_species,
            "edges": len(edges), "compartments": len(comps)}


def prepare(args) -> int:
    root = repo_root(args.root)
    _, contract = _common(root)
    models_dir = root / "panel" / "models"
    jobs: list[tuple[str, pathlib.Path | None]] = []
    if args.panel:
        if not models_dir.is_dir():
            raise SystemExit(f"{models_dir} does not exist: give the SBML files to prepare")
        jobs = [(p.stem, None) for p in sorted(models_dir.glob("*.xml"))]
    for i, src in enumerate(args.sbml):
        name = args.name[i] if i < len(args.name) else pathlib.Path(src).stem
        jobs.append((name, pathlib.Path(src)))
    if not jobs:
        raise SystemExit("nothing to prepare: give SBML files, or --panel to redo panel/models/")

    rows = []
    for model, src in jobs:
        dest = models_dir / f"{model}.xml"
        removed, origin = 0, dest
        if src is not None:
            models_dir.mkdir(parents=True, exist_ok=True)
            origin = src
            if args.strip_groups:
                removed = contract.strip_groups(str(src), str(dest), args.drop_note_key)
            else:
                shutil.copy2(src, dest)
        counts = features(root, model, dest)
        rows.append({
            "model": model, "source": str(origin), "source_sha256": sha256(origin),
            "groups_removed": removed, "species": counts["species"],
            "reactions": counts["reactions"], "compartments": counts["compartments"],
            "prepared": "groups removed" if removed else ("unchanged" if src is not None else "features only"),
        })
        print(f"{model}: {counts['species']} species, {counts['reactions']} reactions, "
              f"{counts['edges']} edges, {counts['compartments']} compartments "
              f"-> panel/features/{model}")

    out = root / "panel" / "prepared.tsv"
    old = {}
    if out.exists():
        with open(out, encoding="utf-8", newline="") as fh:
            old = {r["model"]: r for r in csv.DictReader(fh, delimiter="\t")}
    old.update({r["model"]: r for r in rows})
    with open(out, "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=PREPARED_COLUMNS, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows([old[k] for k in sorted(old)])
    print(f"{len(rows)} model(s) prepared -> panel/prepared.tsv "
          f"(panel/panel.tsv records the models the published results were made from)")
    return 0


def options(args) -> int:
    root = repo_root(args.root)
    sys.path.insert(0, str(root / "layout"))
    import options as OPT  # noqa: PLC0415
    OPT.main()
    return 0


def models(args) -> int:
    root = repo_root(args.root)
    panel = root / "panel" / "panel.tsv"
    with open(panel, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = ["model", "species", "reactions", "compartments", "source_sha256"]
    print("  ".join(f"{c:<14s}" if c != "source_sha256" else c for c in cols))
    for r in rows:
        print("  ".join(f"{r[c]:<14s}" if c != "source_sha256" else r[c][:16] for c in cols))
    print(f"\n{len(rows)} models; source path and full SHA-256 of each in {panel.relative_to(root)}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="big-maps", description=__doc__.split("\n")[0])
    ap.add_argument("--version", action="version", version=f"big-maps {__version__}")
    ap.add_argument("--root", help="the repository checkout to work in (default: the current directory)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="models into panel/models/, feature tables into panel/features/")
    p.add_argument("sbml", nargs="*", help="SBML files to add to the panel")
    p.add_argument("--name", action="append", default=[], help="model name for each SBML file (default: its stem)")
    p.add_argument("--panel", action="store_true", help="also redo every model already in panel/models/")
    p.add_argument("--strip-groups", action="store_true",
                   help="remove the groups the file already carries (common/contract.py strip-groups)")
    p.add_argument("--drop-note-key", action="append", default=[],
                   help="with --strip-groups: also drop this notes key")
    p.set_defaults(func=prepare)

    p = sub.add_parser("options", help="the layout options, and which can run here")
    p.set_defaults(func=options)

    p = sub.add_parser("models", help="the panel, from panel/panel.tsv")
    p.set_defaults(func=models)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
