import importlib
import subprocess
import sys
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.table import QTable
from astropy.time import Time

import tellurion as tell
horizons_reader = importlib.import_module("tellurion.core.horizons")


def _tdb_iso(utc):
    tdb = Time(utc, scale="utc").tdb
    tdb.precision = 3
    return tdb.isot


def _vectors():
    return QTable(
        {
            "targetname": ["ORION (spacecraft)"] * 2,
            "datetime_jd": [2461131.5, 2461131.5069444444] * u.day,
            "x": [1.0, 2.0] * u.au,
            "y": [2.0, 3.0] * u.au,
            "z": [3.0, 4.0] * u.au,
            "vx": [1.0, 2.0] * u.au / u.day,
            "vy": [2.0, 3.0] * u.au / u.day,
            "vz": [3.0, 4.0] * u.au / u.day,
        }
    )


def test_horizons_returns_pvt_with_center_frame_and_tdb_epochs(monkeypatch):
    seen = {}

    class FakeHorizons:
        def __init__(self, **kwargs):
            seen["query"] = kwargs

        def vectors(self, **kwargs):
            seen["vectors"] = kwargs
            return _vectors()

    monkeypatch.setattr(
        horizons_reader.importlib,
        "import_module",
        lambda name: SimpleNamespace(Horizons=FakeHorizons),
    )

    pvt = tell.horizons(
        target_id="-117",
        center="@399",
        start_time="2026-04-01",
        stop_time="2026-04-05",
        step="10m",
    )

    assert pvt.frame == tell.FrameDescriptor(
        "Earth", "ICRF/Earth-equator", kind=tell.FrameKind.INERTIAL
    )
    assert pvt.name == "ORION"
    assert pvt.time.scale == "tdb"
    assert len(pvt) == 2
    assert pvt.position_vector.unit == u.km
    assert pvt.velocity_vector.unit == u.km / u.s
    assert seen["query"] == {
        "id": "-117",
        "location": "@399",
        "epochs": {
            "start": _tdb_iso("2026-04-01"),
            "stop": _tdb_iso("2026-04-05"),
            "step": "10m",
        },
    }
    assert seen["vectors"] == {"refplane": "earth", "aberrations": "geometric"}
    np.testing.assert_allclose(
        pvt.position_vector[0].to_value(u.km),
        (1 * u.au).to_value(u.km) * np.array([1, 2, 3]),
    )


def test_default_target_is_artemis_ii_orion(monkeypatch):
    seen = {}

    class FakeHorizons:
        def __init__(self, **kwargs):
            seen.update(kwargs)

        def vectors(self, **kwargs):
            return _vectors()

    monkeypatch.setattr(
        horizons_reader.importlib,
        "import_module",
        lambda name: SimpleNamespace(Horizons=FakeHorizons),
    )

    tell.horizons()

    assert seen["id"] == "-1024"


@pytest.mark.parametrize(
    ("center", "origin"),
    [
        ("@301", "Moon"),
        ("@0", "Solar System Barycenter"),
        ("@3", "Earth-Moon Barycenter"),
        ("500@399", "500@399"),
    ],
)
def test_center_aliases_are_preserved(center, origin, monkeypatch):
    monkeypatch.setattr(
        horizons_reader.importlib,
        "import_module",
        lambda name: SimpleNamespace(
            Horizons=lambda **kwargs: SimpleNamespace(
                vectors=lambda **options: _vectors()
            )
        ),
    )
    pvt = tell.horizons(center=center, refplane="ecliptic")
    assert pvt.frame.origin == origin
    assert pvt.frame.axes == "ICRF/Ecliptic"
    assert pvt.frame.kind is tell.FrameKind.INERTIAL


def test_body_plane_is_custom_and_invalid_plane_is_rejected(monkeypatch):
    monkeypatch.setattr(
        horizons_reader.importlib,
        "import_module",
        lambda name: SimpleNamespace(
            Horizons=lambda **kwargs: SimpleNamespace(
                vectors=lambda **options: _vectors()
            )
        ),
    )
    pvt = tell.horizons(refplane="body")
    assert pvt.frame.kind is tell.FrameKind.CUSTOM
    with pytest.raises(ValueError, match="refplane"):
        tell.horizons(refplane="equatorial")


def test_missing_astroquery_has_optional_install_message(monkeypatch):
    original = importlib.import_module

    def missing_astroquery(name):
        if name == "astroquery.jplhorizons":
            raise ModuleNotFoundError("missing", name="astroquery")
        return original(name)

    monkeypatch.setattr(
        horizons_reader.importlib, "import_module", missing_astroquery
    )
    with pytest.raises(ImportError, match=r"tellurion\[horizons\]"):
        tell.horizons()


def test_tellurion_import_does_not_import_astroquery():
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, tellurion; assert 'astroquery' not in sys.modules",
        ],
        check=True,
    )


def test_epochs_accept_any_time_form_and_are_sent_as_iso8601_tdb(monkeypatch):
    seen = {}

    class FakeHorizons:
        def __init__(self, **kwargs):
            seen.update(kwargs)

        def vectors(self, **kwargs):
            return _vectors()

    monkeypatch.setattr(
        horizons_reader.importlib,
        "import_module",
        lambda name: SimpleNamespace(Horizons=FakeHorizons),
    )
    start = tell.abstime("2026-04-02T01:30:00")
    tell.horizons(start_time=start, stop_time="2026-04-03T01:30:00")

    assert seen["epochs"]["start"] == _tdb_iso("2026-04-02T01:30:00")
    assert seen["epochs"]["stop"] == _tdb_iso("2026-04-03T01:30:00")
    assert seen["epochs"]["start"].startswith("2026-04-02T01:31:09.")
