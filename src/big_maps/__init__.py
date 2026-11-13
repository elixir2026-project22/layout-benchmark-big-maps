"""layout-benchmark-big-maps: grouping and recursive layout of big pathway maps.

The working code of this repository is run as scripts from the repository root
(``python layout/placement.py ...``); this package holds only what a clean clone
needs before those scripts can run: :func:`repo_root` and the ``big-maps``
command line (``big_maps.__main__``), whose ``prepare`` step turns SBML models
into the feature tables every grouping strategy reads.
"""
from __future__ import annotations

import os
import pathlib

__version__ = "0.1.0"

#: the file that identifies a repository root
MARKER = pathlib.Path("common") / "sbml_features.py"


def repo_root(explicit: str | os.PathLike | None = None) -> pathlib.Path:
    """The root of the repository checkout to work in.

    An *explicit* root is authoritative: it must hold ``common/sbml_features.py``
    or the call fails. Otherwise the order is ``$BIG_MAPS_ROOT``, the current
    directory and its parents, then the checkout this module was imported from
    (which exists when the package is installed in editable mode or run from
    ``src/``).
    """
    if explicit:
        root = pathlib.Path(explicit)
        if not (root / MARKER).exists():
            raise SystemExit(f"{root} is not a layout-benchmark-big-maps checkout (no {MARKER})")
        return root.resolve()
    candidates: list[pathlib.Path] = []
    env = os.environ.get("BIG_MAPS_ROOT")
    if env:
        candidates.append(pathlib.Path(env))
    here = pathlib.Path.cwd().resolve()
    candidates.extend([here, *here.parents])
    candidates.append(pathlib.Path(__file__).resolve().parents[2])
    for c in candidates:
        if (c / MARKER).exists():
            return c.resolve()
    raise SystemExit(
        "no repository root found (no common/sbml_features.py in the current directory, its parents, "
        "or $BIG_MAPS_ROOT): run big-maps from a clone of layout-benchmark-big-maps"
    )
