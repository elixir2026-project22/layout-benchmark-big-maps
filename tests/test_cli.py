"""The ``big-maps`` command line."""
from __future__ import annotations

import pytest

from big_maps import __version__, repo_root
from big_maps.__main__ import main


def test_repo_root_finds_the_checkout(root, monkeypatch):
    monkeypatch.chdir(root / "layout")
    assert repo_root() == root
    assert repo_root(root) == root


def test_repo_root_refuses_a_root_that_is_not_one(tmp_path):
    with pytest.raises(SystemExit):
        repo_root(tmp_path)


def test_repo_root_falls_back_to_the_imported_checkout(root, tmp_path, monkeypatch):
    """Away from any checkout, the one this package was imported from is used."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("BIG_MAPS_ROOT", raising=False)
    assert repo_root() == root


def test_version():
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0
    assert __version__ == "0.1.0"


def test_options_lists_every_option(root, capsys):
    import options as OPT

    assert main(["--root", str(root), "options"]) == 0
    out = capsys.readouterr().out
    assert len(OPT.OPTIONS) == 8
    for name in OPT.OPTIONS:
        assert name in out


def test_models_reads_the_panel_record(root, capsys):
    assert main(["--root", str(root), "models"]) == 0
    out = capsys.readouterr().out
    for model in ("brenda_merged", "e_coli_core", "R-HSA-77289", "skm_pss"):
        assert model in out


def test_prepare_needs_something_to_do(root):
    with pytest.raises(SystemExit):
        main(["--root", str(root), "prepare"])


def test_prepare_writes_feature_tables(root, tiny_model, tmp_path):
    """prepare on a scratch checkout: the model is copied and its tables written."""
    import shutil

    scratch = tmp_path / "clone"
    shutil.copytree(root / "common", scratch / "common")
    assert main(["--root", str(scratch), "prepare", str(tiny_model), "--name", "tiny"]) == 0
    feats = scratch / "panel" / "features" / "tiny"
    assert (feats / "nodes.tsv").exists()
    assert (feats / "edges.tsv").exists()
    assert (feats / "compartments.tsv").exists()
    assert (scratch / "panel" / "models" / "tiny.xml").exists()
    prepared = (scratch / "panel" / "prepared.tsv").read_text(encoding="utf-8")
    assert "tiny" in prepared and "source_sha256" in prepared
