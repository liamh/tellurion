"""Tests for the space-track query interface using a mocked client."""

import datetime
from unittest import mock

import astropy.units as u
import numpy as np
import pytest

import tellurion as tell


def _rec(catid, name, epoch, inc="51.6", bstar=None):
    return {
        "OBJECT_NAME": name, "OBJECT_TYPE": "PAYLOAD", "NORAD_CAT_ID": str(catid),
        "OBJECT_ID": "1998-067A", "EPOCH": epoch,
        "SEMIMAJOR_AXIS": "6780.0", "ECCENTRICITY": "0.0005",
        "INCLINATION": inc, "RA_OF_ASC_NODE": "10.0",
        "ARG_OF_PERICENTER": "20.0", "MEAN_ANOMALY": "30.0",
        "MEAN_MOTION": "15.5", "MEAN_MOTION_DOT": "0.0001",
        "MEAN_MOTION_DDOT": "0", "PERIOD": "92.9", "PERIAPSIS": "400",
        "APOAPSIS": "410", "BSTAR": bstar, "MEAN_ELEMENT_THEORY": "SGP4",
        "TLE_LINE1": "1 " + epoch, "TLE_LINE2": "2 " + epoch,
    }


def _client(records):
    c = mock.Mock()
    c.gp.return_value = records
    c.gp_history.return_value = records
    return c


def test_latest_scalar_keyed_by_catid():
    c = _client([_rec(25544, "ISS", "2026-01-01T00:00:00.000000"),
                 _rec(41335, "SENTINEL 3A", "2026-01-02T00:00:00.000000")])
    out = tell.spacetrack_latest([25544, 41335], c)
    assert set(out) == {25544, 41335}
    assert out[25544].time.isscalar
    assert out[25544].tle[0].startswith("1 ")
    assert out[41335].scdata["name"] == "SENTINEL 3A"
    assert c.gp.call_args.kwargs["norad_cat_id"] == [25544, 41335]


def test_history_series_sorted_vectors():
    c = _client([_rec(25544, "ISS", "2026-01-03T00:00:00.000000", inc="51.7"),
                 _rec(25544, "ISS", "2026-01-01T00:00:00.000000", inc="51.5")])
    out = tell.spacetrack_history(
        25544, tell.abstime("2026-01-01T00:00:00"),
        "2026-01-04T00:00:00", c)
    ser = out[25544]
    assert ser.time.shape == (2,)
    assert ser.time[0] < ser.time[1]
    assert ser.els["inc"].unit == u.deg
    assert list(ser.els["inc"].value) == [51.5, 51.7]
    assert ser.tle[0] == ["1 2026-01-01T00:00:00.000000",
                          "1 2026-01-03T00:00:00.000000"]
    rng = c.gp_history.call_args.kwargs["epoch"]
    assert "2026-01-01" in str(rng) and "2026-01-04" in str(rng)


def test_missing_fields_become_nan():
    c = _client([_rec(1, "X", "2026-01-01T00:00:00.000000", bstar=None)])
    assert np.isnan(tell.spacetrack_latest(1, c)[1].els["B"].value)


def test_query_other_class_returns_raw():
    c = mock.Mock()
    c.satcat.return_value = [{"NORAD_CAT_ID": "1"}]
    assert tell.spacetrack_query("satcat", c, norad_cat_id=1) == [{"NORAD_CAT_ID": "1"}]


def test_series_hdf5_roundtrip(tmp_path):
    fsc = pytest.importorskip("fsc.hdf5_io")
    pytest.importorskip("astropy_hdf5io")
    c = _client([_rec(25544, "ISS", "2026-01-01T00:00:00.000000"),
                 _rec(25544, "ISS", "2026-01-03T00:00:00.000000")])
    ser = tell.spacetrack_history(25544, "2026-01-01", "2026-01-04", c)[25544]
    f = str(tmp_path / "s.h5")
    fsc.save(ser, f)
    back = fsc.load(f)
    assert list(back.tle) == ser.tle
    assert back.time.shape == (2,)
    assert (back.els["inc"] == ser.els["inc"]).all()
