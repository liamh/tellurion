"""Tests for the Sphinx project version configuration."""

import importlib.metadata
from pathlib import Path
import runpy

import pytest


def test_sphinx_version_uses_distribution_metadata(monkeypatch):
    pytest.importorskip("sphinx_astropy.conf.v2")
    original = importlib.metadata.version
    monkeypatch.setattr(
        importlib.metadata, "version",
        lambda name: "1.2.3" if name == "tellurion" else original(name),
    )
    config = runpy.run_path(str(Path(__file__).resolve().parents[1] /
                                "docs" / "source" / "conf.py"))
    assert config["release"] == "1.2.3"
    assert config["version"] == "1.2.3"
