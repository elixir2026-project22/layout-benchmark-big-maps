This is what each file and folder of this repository is.

The working code is run as scripts from the repository root (`python layout/placement.py ...`), which is
how it was written and how the documented commands read. `src/` holds only the installable part: the
`big-maps` command line, whose `prepare` step is what a clean clone needs before any script can run.

* Code
  * `common/` - what every stage shares: `sbml_features.py` (the feature extractor: SBML to
    `nodes.tsv` / `edges.tsv` / `compartments.tsv`), `compartments.py` (the compartment containment tree),
    `contract.py` (the contract validator, the SBML Level 3 groups writer and `strip-groups`),
    `SBO_OBO.obo` and `go_cellular_component.json` (the identifier tables the naming rules use).
  * `strategies/<name>/` - one grouping strategy: `strategy.py`, `RULES.md` (the numbered rules its
    output refers to), `NOTES.md`, and `out/<model>/` with the grouping it produced. `strategies/topology/`
    holds one program and three variants (`topology_leiden`, `topology_infomap`, `topology_walktrap`),
    each with its own `out/` and `stats/`.
  * `baseline/` - the project's current grouping converted to the same format, for comparison.
  * `rules/` - the three independent rule proposals: `rules_*.md` and the matching `metrics_*.py` with
    their tests.
  * `scoring/` - `score_grouping.py` (the merged rulebook of `scoring/RULES.md`), `collect_groupings.py`,
    `sensitivity.py`, their tests, and the results (`scores.csv`, `scores_summary.csv`,
    `sensitivity/`).
  * `examples/` - the grouping stage's renderer (`render_grouping.py`), its gallery builder and the
    gallery itself (`index.html`, `gallery.tsv`, `svg/`).
  * `layout/` - the joining stage: `options.py` (the eight layout options), `placement.py` (the recursive
    placement), `render_layout.py`, `build_gallery.py`, its `README.md` and `PROTOCOL.md`, and the gallery
    in `examples/`.
  * `src/big_maps/` - the `big-maps` command line (`prepare`, `options`, `models`), installed by
    `pip install -e .`; also runnable as `python -m big_maps`.
  * `tests/` - the package's own tests: `conftest.py` builds a four-species SBML model so the layout
    stage can be tested in a clean clone, `test_cli.py`, `test_layout.py`, `test_gallery_page.py`.
    The other tests live next to the code they test (`rules/`, `scoring/`); `pyproject.toml` lists all
    three as `testpaths`.
* Data that is not here
  * `panel/panel.tsv` - the provenance record of the models: source path, SHA-256, size and what was done
    to each. The models themselves (`panel/models/`) and their feature tables (`panel/features/`) are not
    in this repository; `big-maps prepare` writes both, and `panel/prepared.tsv` records what it wrote.
    `.gitignore` lists every excluded path.
* Documentation
  * `README.md` - what the repository is, the two stages, the options, and the commands that reproduce
    everything.
  * `CONTRACT.md` - what a grouping strategy reads and writes, and the rules every strategy keeps.
  * `layout/PROTOCOL.md` - the five rules of the placement, and what a run writes.
  * `doc/REPORT.md` - the result of the grouping comparison; `doc/index.md` lists the documentation.
  * `PUSH.md` - the commands that push this repository to its GitHub remote.
  * `COPYING.md` (copyright holders), `LICENSE.md` (MIT), `CITATION.cff` (citation metadata).
* Build and automation
  * `pyproject.toml` - project metadata, the dependencies, the `big-maps` entry point, the pytest and
    ruff configuration.
  * `uv.lock` - the resolved version of every dependency, as in layout-benchmark.
  * `.github/workflows/ci.yml` - the tests and the gallery-page check on every push and pull request to
    `main`; `.github/dependabot.yml` keeps the workflow's actions up to date.
  * `.pre-commit-config.yaml` - ruff, YAML/TOML checks and a file size limit of 1 MB (the gallery SVGs of
    the big models reach 0.83 MB). Not applied to the carried-over code: it is opt-in
    (`pre-commit install`) and CI does not run it.
  * `.gitattributes` - LF line endings in the repository; `.gitignore` - the excluded data paths and the
    usual Python noise.
