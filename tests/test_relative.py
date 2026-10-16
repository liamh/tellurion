"""Tests for tellurion.ork.relative (ntw, lvlh and inverses)."""

from __future__ import annotations

import numpy as np
import pytest
import astropy.units as u

import tellurion as tell

pytestmark = pytest.mark.requires_orekit_data


def _pos(p):
    return np.asarray(p.position_vector.to_array(), dtype=float)


def _vel(p):
    return np.asarray(p.velocity_vector.to_array(), dtype=float)


def _offset(pvt, dpos_km):
    pv = np.concatenate([_pos(pvt) + dpos_km, _vel(pvt)])
    return tell.pvtcart(pv, pvt.time)


ORBITS = ["leo1", "leo2", "geo1", "ell1", "ell2", "vang1", "gps1"]


@pytest.fixture(params=ORBITS)
def ref(request):
    return request.getfixturevalue(request.param)["pvt"]


@pytest.mark.parametrize("func", [tell.ntw, tell.lvlh])
def test_identical_states_give_zero(ref, func):
    rel = func(ref, ref, tell.siunits)
    assert np.allclose(np.asarray(rel.position.position_vector.to_array(), dtype=float), 0.0, atol=1e-6)
    assert np.allclose(np.asarray(rel.velocity_vector.to_array(), dtype=float), 0.0, atol=1e-9)


@pytest.mark.parametrize("func", [tell.ntw, tell.lvlh])
def test_range_is_frame_invariant(ref, func):
    sec = _offset(ref, np.array([1.0, -2.0, 0.5]))
    rel = func(sec, ref, tell.siunits)
    assert np.linalg.norm(np.asarray(rel.position.position_vector.to_array(), dtype=float)) == pytest.approx(
        np.linalg.norm([1.0, -2.0, 0.5]), rel=1e-9
    )


def test_names(leo1):
    sec = _offset(leo1["pvt"], np.array([1.0, 0.0, 0.0]))
    assert tell.ntw(sec, leo1["pvt"]).name == "NTW"
    assert tell.lvlh(sec, leo1["pvt"]).name == "LVLH"


def test_lvlh_radial_offset_is_first_axis(leo1):
    ref = leo1["pvt"]
    rhat = _pos(ref) / np.linalg.norm(_pos(ref))
    rel = tell.lvlh(_offset(ref, 1.0 * rhat), ref, tell.siunits)
    x, y, z = np.asarray(rel.position.position_vector.to_array(), dtype=float).ravel()
    assert x == pytest.approx(1.0, rel=1e-9)
    assert abs(y) < 1e-9 and abs(z) < 1e-9


def test_ntw_crosstrack_offset_is_third_axis(leo1):
    ref = leo1["pvt"]
    h = np.cross(_pos(ref), _vel(ref))
    what = h / np.linalg.norm(h)
    rel = tell.ntw(_offset(ref, 1.0 * what), ref, tell.siunits)
    x, y, z = np.asarray(rel.position.position_vector.to_array(), dtype=float).ravel()
    assert z == pytest.approx(1.0, rel=1e-9)
    assert abs(x) < 1e-9 and abs(y) < 1e-9


def test_ntw_intrack_offset_is_second_axis(leo1):
    ref = leo1["pvt"]
    that = _vel(ref) / np.linalg.norm(_vel(ref))
    rel = tell.ntw(_offset(ref, 1.0 * that), ref, tell.siunits)
    x, y, z = np.asarray(rel.position.position_vector.to_array(), dtype=float).ravel()
    assert y == pytest.approx(1.0, rel=1e-9)
    assert abs(x) < 1e-9 and abs(z) < 1e-9


@pytest.mark.parametrize(
    "fwd, inv",
    [(tell.ntw, tell.inv_ntw), (tell.lvlh, tell.inv_lvlh)],
)
def test_roundtrip(ref, fwd, inv):
    sec = _offset(ref, np.array([3.0, -1.0, 2.0]))
    back = inv(fwd(sec, ref, tell.siunits), ref)
    assert _pos(back) == pytest.approx(_pos(sec), abs=1e-6)
    assert _vel(back) == pytest.approx(_vel(sec), abs=1e-9)


def test_curtis_example_7_1():
    t = tell.abstime("2025-01-01T00:00:00")
    a = tell.pvtcart([-266.77, 3865.8, 5426.2, -6.4836, -3.6198, 2.4156], t)
    b = tell.pvtcart([-5890.71, -2979.8, 1792.2, 0.93583, -5.2403, -5.5009], t)
    rel = tell.lvlh(b, a, tell.siunits)
    assert np.linalg.norm(_pos(rel)) == pytest.approx(
        np.linalg.norm(_pos(b) - _pos(a)), rel=1e-9
    )


def test_time_series_roundtrip(leo1):
    ref0 = leo1["pvt"]
    gen_ref = tell.prepare(ref0, 1 * u.hour)
    gen_sec = tell.prepare(_offset(ref0, np.array([0.5, 0.0, 0.0])), 1 * u.hour)
    steps = np.linspace(5 * u.minute, 30 * u.minute, 6)
    ref = tell.propagate(gen_ref, steps, True, output="pvt")
    sec = tell.propagate(gen_sec, steps, True, output="pvt")
    back = tell.inv_lvlh(tell.lvlh(sec, ref, tell.siunits), ref)
    assert np.allclose(_pos(back), _pos(sec), atol=1e-6)
    assert np.allclose(_vel(back), _vel(sec), atol=1e-9)