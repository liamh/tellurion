# conftest.py
"""Shared test fixtures for orbit and element tests."""

from __future__ import annotations
import numpy as np
import astropy.units as u
import tellurion as tell
import pytest
import os
from pathlib import Path



def keppvt(dictels, time):
    """From the dictionary of element values, compute the Kepler
    elements (sma, ecc, inc, argper, raan, and ma or ta), the PVT
    transformation of those elements, and the spherical coordinates."""
    d = {'els': dictels}
    d['kep'] = tell.kepler(tell.allplane(d['els']), time)
    d['pvt'] = tell.pvt(d['kep'])
    d['pvt'].spherical
    return d


# -----------------------------------------------------------------------
# Orbit Fixtures
# -----------------------------------------------------------------------

@pytest.fixture
def leo1():
    """A LEO circular orbit like the ISS"""
    return keppvt({"altper": 350*u.km, "altapo": 350*u.km,
                   "inc":55.0*u.deg, "argper": 120.0*u.deg,
                   "raan": 20.0*u.deg, "ma": 30.0*u.deg},
                  tell.abstime('2026-01-01 05:55:00'))

@pytest.fixture
def leo2():
    """A LEO near-circular orbit like the HST"""
    return keppvt({"altper": 525.0*u.km, "altapo": 555.0*u.km,
                   "inc":28.5*u.deg, "argper": 40.0*u.deg,
                   "raan": 40.0*u.deg, "ma": -30.0*u.deg},
                  tell.abstime('2026-01-01 09:30:00'))

@pytest.fixture
def geo1():
    """A GEO orbit"""
    return keppvt({"memo":1.0*u.rev/u.sday,
                   "ecc":0.0,
                   "inc":0.0*u.deg, "argper": 120.0*u.deg,
                   "raan": 0.0*u.deg, "ma": 0.0*u.deg},
                  tell.abstime('2026-01-01 20:30:00'))

@pytest.fixture
def ell1():
    """A GEO transfer orbit"""
    return keppvt({"altper": 350*u.km, "altapo": tell.sma(1.0, True),
                   "inc":0.0*u.deg, "argper": 120.0*u.deg,
                   "raan": 0.0*u.deg, "ma": 90.0*u.deg},
                  tell.abstime('2026-01-01 05:55:00'))

@pytest.fixture
def ell2():
    """An elliptical orbit"""
    return keppvt({"altper":160*u.km, "altapo":20250*u.km,
                   "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg},
                  tell.abstime('2026-01-01 08:05:00'))

@pytest.fixture
def vang1():
    """Like Vanguard 1"""
    return keppvt({"altper": 600.0*u.km, "altapo": 12000.0*u.km,
                   "inc":36.0*u.deg, "argper": 140.0*u.deg,
                   "raan": 0.0*u.deg, "ma": 100.0*u.deg},
                  tell.abstime('2026-01-01 14:45:00'))

@pytest.fixture
def gps1():
    """A semisynchronous orbit like GPS"""
    return keppvt({"sma": tell.sma(2.0), "ecc": 0.0*u.dimensionless_unscaled,
                   "inc":55.0*u.deg, "argper": 0.0*u.deg,
                   "raan": 120.0*u.deg, "ma": 77.0*u.deg},
                  tell.abstime('2026-01-01 12:20:00'))

@pytest.fixture
def all_orbits(leo1, leo2, geo1, ell1, ell2, vang1, gps1):
    """All orbit fixtures as a list"""
    return [leo1, leo2, geo1, ell1, ell2, vang1, gps1]

# -----------------------------------------------------------------------
# Orekit data fixtures
# -----------------------------------------------------------------------

def _candidate_paths(repo_root: Path) -> list[Path]:
    return [
        # explicit override
        Path(os.environ["OREKIT_DATA_PATH"]).expanduser().resolve()
        if "OREKIT_DATA_PATH" in os.environ
        else Path("__missing__"),
        # package-bundled
        repo_root / "tellurion" / "ork" / "_orekit_data",
        repo_root / "tellurion" / "ork" / "_orekit_data.zip",
        # cache
        Path.home() / ".cache" / "tellurion",
        Path.home() / ".cache" / "tellurion" / "orekit-data.zip",
    ]


def _find_orekit_data(repo_root: Path) -> Path | None:
    for p in _candidate_paths(repo_root):
        if p.exists():
            return p
    return None


@pytest.fixture(scope="session")
def orekit_data_path() -> Path:
    repo_root = Path(__file__).resolve().parents[1]
    p = _find_orekit_data(repo_root)
    if p is None:
        pytest.skip(
            "Orekit data not found. Checked tellurion/ork/_orekit_data, "
            "tellurion/ork/_orekit_data.zip, ~/.cache/tellurion, "
            "~/.cache/tellurion/orekit-data.zip, and OREKIT_DATA_PATH."
        )
    return p


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_orekit_data: mark test as requiring local Orekit data files",
    )


def pytest_collection_modifyitems(config, items):
    repo_root = Path(__file__).resolve().parents[1]
    data_path = _find_orekit_data(repo_root)
    if data_path is not None:
        return

    skip_marker = pytest.mark.skip(
        reason="Orekit data not found (set OREKIT_DATA_PATH or install cached/package data)."
    )
    for item in items:
        if "requires_orekit_data" in item.keywords:
            item.add_marker(skip_marker)
