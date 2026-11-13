"""Collect every strategy's output into one folder the scorer can walk.

score_grouping.py --all expects <strategies>/<strategy>/<model>/groups.tsv. In this repository a strategy
keeps its outputs in strategies/<name>/out/<model>/ (and the topology variants one level deeper), so this
script copies the contract files into scoring/groupings/<strategy>/<model>/.

Run:
    python scoring/collect_groupings.py
    python scoring/score_grouping.py --all --panel panel/features --strategies scoring/groupings \
        --baseline baseline --heldout heldout --out scoring/scores.csv --summary scoring/scores_summary.csv
"""
from __future__ import annotations

import pathlib
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
SOURCES = {
    "type_location": "strategies/type_location/out",
    "location_type": "strategies/location_type/out",
    "ontology": "strategies/ontology/out",
    "topology_leiden": "strategies/topology/topology_leiden/out",
    "topology_infomap": "strategies/topology/topology_infomap/out",
    "topology_walktrap": "strategies/topology/topology_walktrap/out",
    "cascade": "strategies/cascade/out",
    "hybrid": "strategies/hybrid/out",
}

if __name__ == "__main__":
    dest_root = ROOT / "scoring" / "groupings"
    shutil.rmtree(dest_root, ignore_errors=True)
    n = 0
    for name, rel in SOURCES.items():
        for model_dir in sorted((ROOT / rel).iterdir()):
            if not (model_dir / "groups.tsv").exists():
                continue
            dest = dest_root / name / model_dir.name
            dest.mkdir(parents=True, exist_ok=True)
            for f in ("groups.tsv", "groups_meta.tsv"):
                shutil.copy(model_dir / f, dest / f)
            n += 1
    print(f"{n} groupings -> {dest_root.relative_to(ROOT)}")
