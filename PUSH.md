# Pushing this repository to GitHub

> **State:** `main` is pushed (`origin/main`). What remains is the website: push the `gh-pages`
> branch and switch GitHub Pages on. Jump to [The website](#the-website).

## The website

The site lives on an orphan `gh-pages` branch - `.nojekyll`, a landing `index.html`, and one folder
per gallery - which is how
[layout-benchmark](https://elixir2026-project22.github.io/layout-benchmark/) publishes its report.
The branch already exists locally, built by `scripts/build_site.py` from the committed galleries
(346 files, about 100 MB; it packs far smaller, as the SVGs compress about ninefold).

```bash
git push -u origin gh-pages
```

Then, in **Settings → Pages**: *Source* = *Deploy from a branch*, branch `gh-pages`, folder
`/ (root)`. The first build takes a minute or two, after which the galleries are at

    https://elixir2026-project22.github.io/layout-benchmark-big-maps/
    https://elixir2026-project22.github.io/layout-benchmark-big-maps/grouping/index.html
    https://elixir2026-project22.github.io/layout-benchmark-big-maps/layout/index.html

To refresh the site after rebuilding a gallery, rebuild it into the worktree and push again:

```bash
python scripts/build_site.py -o ../_gh-pages
git -C ../_gh-pages add -A
git -C ../_gh-pages commit -m "Site: galleries at $(git rev-parse --short main)"
git -C ../_gh-pages push
```

`git worktree remove ../_gh-pages` drops the worktree; the branch stays. If the worktree is gone,
recreate it with `git worktree add ../_gh-pages gh-pages`.

## The topic

```bash
gh repo edit elixir2026-project22/layout-benchmark-big-maps --add-topic biohackeu26
```

or the gear next to *About* on the repository's main page.

---

## Reference: the original push



This checkout has one commit and **no remote**. Nothing has been pushed: the session that assembled it had
no GitHub credential, so it could neither create the repository nor push.

The remote repository is created by hand: **public**, in the
[`elixir2026-project22`](https://github.com/orgs/elixir2026-project22/repositories) organisation, named
`layout-benchmark-big-maps`, **without** a README, licence or `.gitignore` (the remote must be empty, or
the first push is rejected as a non-fast-forward).

## Before you push: the commit author

The single commit was made without a configured git identity, so git derived one from the machine, the
same way the commits of the repository this was assembled from were made. Check it, and set it to the
identity you want to appear publicly:

```bash
git log -1 --format='%an <%ae>'          # what the commit carries now
git config user.name "Your Name"         # for this repository only
git config user.email "you@example.org"  # e.g. your GitHub noreply address
git commit --amend --reset-author --no-edit
```

## Push

Run these from the repository root:

```bash
cd "layout-benchmark-big-maps"
git remote add origin https://github.com/elixir2026-project22/layout-benchmark-big-maps.git
git branch -M main
git push -u origin main
```

Over SSH instead of HTTPS, the remote is
`git@github.com:elixir2026-project22/layout-benchmark-big-maps.git`.

`main` is the branch name the CI workflow triggers on (`.github/workflows/ci.yml`), and is the convention
of [layout-benchmark](https://github.com/elixir2026-project22/layout-benchmark).

The commit holds 547 files and 103.6 MB of working tree, which pack into 9.8 MiB of git objects (the
gallery SVGs compress about 9-fold) - that is what the push uploads and what a clone downloads. No single
file is near GitHub's 100 MB per-file limit; the largest is 0.83 MB.

## The topic

The repository carries the single topic `biohackeu26`. Either set it in the GitHub UI - the repository's
main page, the gear next to *About*, *Topics* - or, with the [GitHub CLI](https://cli.github.com/) signed
in (`gh auth login`):

```bash
gh repo edit elixir2026-project22/layout-benchmark-big-maps --add-topic biohackeu26
```

## Check afterwards

```bash
git remote -v                     # origin, both fetch and push
git log --oneline -1              # the commit that was pushed
git status                        # clean, and no ignored data path staged
gh repo view elixir2026-project22/layout-benchmark-big-maps   # visibility and topics
```

The CI workflow runs on the first push to `main`: it installs Graphviz, runs `pytest`, and checks that
`layout/examples/index.html` is still what `layout/build_gallery.py --no-render` writes from the committed
`layout/examples/gallery.tsv`.
