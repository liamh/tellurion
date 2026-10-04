"""Tests for Orekit data location helper."""

from tellurion import orekit_data as od


def _make_data(path):
    path.mkdir(parents=True)
    (path / "UTC-TAI.history").write_text("x")
    return path


def test_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv(od.ENV_VAR, str(tmp_path))
    assert od.orekit_data_path() == tmp_path


def test_bundled_preferred(monkeypatch, tmp_path):
    monkeypatch.delenv(od.ENV_VAR, raising=False)
    bundled = _make_data(tmp_path / "b")
    monkeypatch.setattr(od, "bundled_data_path", lambda: bundled)
    assert od.orekit_data_path() == bundled


def test_cache_fallback(monkeypatch, tmp_path):
    monkeypatch.delenv(od.ENV_VAR, raising=False)
    monkeypatch.setenv("TELLURION_CACHE", str(tmp_path))
    monkeypatch.setattr(od, "bundled_data_path", lambda: None)
    cached = _make_data(tmp_path / "orekit-data")
    assert od.orekit_data_path() == cached


def test_download_fallback(monkeypatch, tmp_path):
    monkeypatch.delenv(od.ENV_VAR, raising=False)
    monkeypatch.setenv("TELLURION_CACHE", str(tmp_path))
    monkeypatch.setattr(od, "bundled_data_path", lambda: None)
    called = []
    monkeypatch.setattr(
        od, "download_orekit_data", lambda d: called.append(d) or d
    )
    assert od.orekit_data_path() == tmp_path / "orekit-data"
    assert called
