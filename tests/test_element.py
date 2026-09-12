# Create element sets and convert to both Cartesian and spherical
# (right ascension, declination, distance) position-velocity specification

# Run with pytest -q test_element.py

import astropy.units as u
import importlib
from types import SimpleNamespace

import numpy as np
import pytest
import tellurion as tell
from astropy.time import Time
from tellurion.core.element import (ElementSetT, iskepels,
                                    isequels, iscircels,
                                    equeltma_names, equeltta_names,
                                    circeltma_names, circeltta_names)

# conftest.py is auto-discovered by pytest, so fixtures are available

def _make_eq_from_pvt(pvtobj):
    """Equinoctial ElementSetT from a PVT via the monkey-patch."""
    return pvtobj.equinoctial()

def _make_circ_from_pvt(pvtobj):
    """Circular ElementSetT from a PVT via the monkey-patch."""
    return pvtobj.circular()

# -----------------------------------------------------------------------
# ElementSetT structure tests
# -----------------------------------------------------------------------

class TestElementSetT:
    def test_is_class_not_namedtuple(self, leo1):
        est = leo1['kep']
        assert isinstance(est, ElementSetT)

    def test_has_elements_attribute(self, leo1):
        est = leo1['kep']
        assert isinstance(est.elements, u.Quantity)

    def test_has_time_attribute(self, leo1):
        est = leo1['kep']
        assert isinstance(est.time, Time)

    def test_iskepels_true_elements_only(self, leo1):
        est = leo1['kep']
        assert isinstance(est, ElementSetT)

    def test_unpack_iteration(self, leo1):
        est = leo1['kep']
        elements, time = est
        assert isinstance(elements, u.Quantity)
        assert isinstance(time, Time)

    def test_getitem(self, leo1):
        est = leo1['kep']
        assert est[0] is est.elements
        assert est[1] is est.time

    def test_pvt_method_delegates_to_ork_converter(self, monkeypatch):
        est = tell.kepler(
            {
                "sma": 7000 * u.km,
                "ecc": 0.01 * u.dimensionless_unscaled,
                "inc": 30 * u.deg,
                "argper": 40 * u.deg,
                "raan": 50 * u.deg,
                "ma": 60 * u.deg,
            },
            tell.abstime("2026-01-01 00:00:00"),
        )
        sentinel = object()
        calls = []

        def fake_import(name):
            calls.append(("import", name))
            return SimpleNamespace(pvt=lambda arg: calls.append(("pvt", arg)) or sentinel)

        monkeypatch.setattr(
            importlib,
            "import_module",
            fake_import,
        )

        assert est.pvt() is sentinel
        assert calls == [("import", "tellurion.ork.element"), ("pvt", est)]

    def test_iskepels_true_with_epoch(self, leo1):
        est = leo1['kep']
        assert iskepels(est, est=True)

    def test_iskepels_false_for_plain_quantity(self, leo1):
        assert not iskepels(1.0*u.km, est=False)

    def test_repr(self, leo1):
        est = leo1['kep']
        r = repr(est)
        assert 'ElementSetT' in r

    def test_equality(self, leo1):
        est1 = leo1['kep']
        est2 = leo1['kep']
        assert est1 == est2

    def test_inequality(self, leo1):
        est1 = leo1['kep']
        est2 = tell.kepler(
            tell.allplane({"altper": 400*u.km, "altapo": 400*u.km,
                           "inc":55.0*u.deg, "argper": 120.0*u.deg,
                           "raan": 20.0*u.deg, "ma": 30.0*u.deg}),
            tell.abstime('2026-01-01 05:55:00'))
        assert est1 != est2

    def test_constructor_rejects_bad_t(self, leo1):
        est = leo1['kep']
        with pytest.raises(TypeError):
            ElementSetT(est.elements, "not-a-time")


# -----------------------------------------------------------------------
# allplane / sma conversion tests
# -----------------------------------------------------------------------

def test_allplane(geo1, ell1):
    """Add fixtures as parameters"""
    smaecc = tell.allplane({"altper":160*u.km, "altapo":20250*u.km,
                            "inc":28.5*u.deg, "argper": 0.0*u.deg,
                            "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    alts = tell.allplane({"sma":8000*u.km, "ecc":0.1*u.dimensionless_unscaled,
                          "inc":45*u.deg, "argper": 120.0*u.deg,
                          "raan": 80.0*u.deg, "ma": 0.0*u.deg})

    np.testing.assert_allclose(smaecc['sma'], 16583.13646*u.km)
    np.testing.assert_allclose(smaecc['ecc'], 0.60573583*u.dimensionless_unscaled)
    np.testing.assert_allclose(alts['altper'], 821.86354*u.km)
    np.testing.assert_allclose(alts['altapo'], 2421.86354*u.km)
    np.testing.assert_allclose(geo1['kep'].elements['sma'], 42164.1696233*u.km)
    # Test conversion to semimajor axis
    np.testing.assert_allclose(tell.sma(1e4*u.s), 10032.11910363*u.km)
    np.testing.assert_allclose(tell.sma(-20*(u.km/u.s)**2), 9965.0110375*u.km)
    np.testing.assert_allclose(ell1['kep'].elements['sma'], 24446.15304165*u.km)


# -----------------------------------------------------------------------
# Position-velocity tests
# -----------------------------------------------------------------------

def test_posvel(leo1, leo2, geo1, ell1, ell2, vang1, gps1):
    """Add all fixtures as parameters"""
    np.testing.assert_allclose(leo1['pvt'].cartesian['position'].si.value,
                               np.array([-6135286.90982537, -179677.30883764, 2755683.36773215]))
    np.testing.assert_allclose(leo1['pvt'].cartesian['velocity'].si.value,
                               np.array([-2308.74628271, -4909.03309813, -5460.30174535]))
    np.testing.assert_allclose(np.array(leo1['pvt'].spherical.si.value.tolist()),
                               np.array([3.17087017e+00, 4.21989266e-01, 6.72813646e+06,
                                         7.88434305e-04, -8.89601706e-04, -2.66453526e-12]),
                               atol=1.0e-7, rtol=0.0)
    np.testing.assert_allclose(leo2['pvt'].cartesian['position'].si.value,
                               np.array([4542282.92387955, 5170057.21334251, 565092.22923442]))
    np.testing.assert_allclose(leo2['pvt'].cartesian['velocity'].si.value,
                               np.array([-5236.82832614, 4199.24378643, 3574.26419318]))
    np.testing.assert_allclose(np.array(leo2['pvt'].spherical.si.value.tolist()),
                               np.array([8.49945142e-01, 8.19279154e-02, 6.90515423e+06,
                                         9.74389282e-04, 5.19462927e-04, -8.25996492e+00]))
    np.testing.assert_allclose(geo1['pvt'].cartesian['position'].si.value,
                               np.array([-21082084.8116475, 36515242.02324963, 0.]))
    np.testing.assert_allclose(geo1['pvt'].cartesian['velocity'].si.value,
                               np.array([-2662.73375321, -1537.33004919, -0.]))
    np.testing.assert_allclose(np.array(geo1['pvt'].spherical.si.value.tolist()),
                               np.array([2.0943951023931953, 0.0, 42164169.623295024,
                                         7.292115855377073e-05, 0.0, 0.0]),
                               atol=1.0e-7, rtol=0.0)
    np.testing.assert_allclose(ell1['pvt'].cartesian['position'].si.value,
                               np.array([3698345.45792049, -34232447.31322606, -0.]))
    np.testing.assert_allclose(ell1['pvt'].cartesian['velocity'].si.value,
                               np.array([2148.20304004, -1494.36501206, 0.]))
    np.testing.assert_allclose(np.array(ell1['pvt'].spherical.si.value.tolist()),
                               np.array([4.82000783e+00, 0.00000000e+00, 3.44316454e+07,
                                         5.73676739e-05, 0.00000000e+00, 1.71646077e+03]))
    np.testing.assert_allclose(ell2['pvt'].cartesian['position'].si.value,
                               np.array([6538136.46, 0., 0.]))
    np.testing.assert_allclose(ell2['pvt'].cartesian['velocity'].si.value,
                               np.array([-0., 8695.15748257, 4721.08531441]))
    np.testing.assert_allclose(np.array(ell2['pvt'].spherical.si.value.tolist()),
                               np.array([0.00000000e+00, 0.00000000e+00, 6.53813646e+06,
                                         1.32991374e-03, 7.22084243e-04, 0.00000000e+00]))
    np.testing.assert_allclose(vang1['pvt'].cartesian['position'].si.value,
                               np.array([3313526.7287987, -12405644.53079198, -9013228.33893749]))
    np.testing.assert_allclose(vang1['pvt'].cartesian['velocity'].si.value,
                               np.array([4321.62378008, -676.43293734, -491.45729632]))
    np.testing.assert_allclose(np.array(vang1['pvt'].spherical.si.value.tolist()),
                               np.array([4.973394318345154, -0.6120236334818134, 15688140.76609011,
                                         0.00031156788681119144, 3.91331269159631e-05, 1730.0341520321017]))
    np.testing.assert_allclose(gps1['pvt'].cartesian['position'].si.value,
                               np.array([-15843456.17405487, -2247776.5833332, 21200462.73723146]))
    np.testing.assert_allclose(gps1['pvt'].cartesian['velocity'].si.value,
                               np.array([1454.40855181, -3518.76365779, 713.82704166]))
    np.testing.assert_allclose(np.array(gps1['pvt'].spherical.si.value.tolist()),
                               np.array([3.28252622e+00, 9.24230194e-01, 2.65617624e+07,
                                         2.30480399e-04, 4.46083005e-05, 1.11022302e-13]),
                               atol=1.0e-7, rtol=1.0e-8)


# -----------------------------------------------------------------------
# Equinoctial: structural tests
# -----------------------------------------------------------------------

class TestEquinoctialStructure:

    def test_returns_elementset_t(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert isinstance(est, ElementSetT)

    def test_has_t_attribute(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert isinstance(est.time, Time)

    def test_els_is_quantity(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert isinstance(est.elements, u.Quantity)

    def test_field_names_mean(self, leo1):
        """Round-trip via mean longitude should produce mean-longitude field names."""
        est = _make_eq_from_pvt(leo1['pvt'])   # default is mean
        assert set(equeltma_names) <= set(est.elements.dtype.names)

    def test_unpack(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        els, t = est
        assert isinstance(els, u.Quantity)
        assert isinstance(t, Time)

    def test_getitem(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert est[0] is est.elements
        assert est[1] is est.time

    def test_repr(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert 'ElementSetT' in repr(est)

    def test_equality(self, leo1):
        est1 = _make_eq_from_pvt(leo1['pvt'])
        est2 = _make_eq_from_pvt(leo1['pvt'])
        assert est1 == est2

    def test_inequality(self, leo1, leo2):
        eq_leo1 = _make_eq_from_pvt(leo1['pvt'])
        eq_leo2 = _make_eq_from_pvt(leo2['pvt'])
        assert eq_leo1 != eq_leo2


# -----------------------------------------------------------------------
# Equinoctial: predicate tests
# -----------------------------------------------------------------------

class TestIsequels:

    def test_true_with_epoch(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert isequels(est, est=True)

    def test_true_els_only(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert isequels(est.elements, est=False)

    def test_false_for_keplerian_est(self, leo1):
        assert not isequels(leo1['kep'])

    def test_false_for_plain_quantity(self):
        assert not isequels(1.0 * u.km, est=False)

    def test_false_for_circular_est(self, leo1):
        circ = _make_circ_from_pvt(leo1['pvt'])
        assert not isequels(circ)

    def test_iskepels_false_for_equinoctial(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert not iskepels(est)


# -----------------------------------------------------------------------
# Equinoctial: numerical consistency
# -----------------------------------------------------------------------

class TestEquinoctialValues:

    def test_sma_matches_keplerian(self, leo1):
        """Equinoctial sma must equal Keplerian sma for the same orbit."""
        eq_sma  = tell.elementval(leo1['pvt'].equinoctial(), 'sma')
        kep_sma = tell.elementval(leo1['pvt'].kepler(),      'sma')
        np.testing.assert_allclose(eq_sma.si.value, kep_sma.si.value, rtol=1e-10)

    def test_circular_orbit_ex_ey_near_zero(self, leo1):
        """For a near-circular orbit (leo1, altper=altapo=350km), ex and ey ≈ 0."""
        est = _make_eq_from_pvt(leo1['pvt'])
        np.testing.assert_allclose(est.elements['ex'].value, 0.0, atol=1e-6)
        np.testing.assert_allclose(est.elements['ey'].value, 0.0, atol=1e-6)

    def test_geo_hx_hy_near_zero(self, geo1):
        """For a zero-inclination GEO orbit, hx=tan(i/2)cos(Ω)≈0, hy≈0."""
        est = _make_eq_from_pvt(geo1['pvt'])
        np.testing.assert_allclose(est.elements['hx'].value, 0.0, atol=1e-8)
        np.testing.assert_allclose(est.elements['hy'].value, 0.0, atol=1e-8)

    def test_elementval_sma_from_equinoctial(self, gps1):
        """elementval dispatches correctly to _eqdict for an equinoctial ElementSetT."""
        est = _make_eq_from_pvt(gps1['pvt'])
        sma = tell.elementval(est, 'sma')
        np.testing.assert_allclose(sma.to(u.km).value,
                                   tell.sma(2.0).to(u.km).value,
                                   rtol=1e-6)

    def test_elementval_list_from_equinoctial(self, leo2):
        """elementval with a list of element names works for equinoctial."""
        est  = _make_eq_from_pvt(leo2['pvt'])
        vals = tell.elementval(est, ['sma', 'ex', 'ey', 'hx', 'hy', 'ml'])
        assert len(vals) == 6
        assert all(isinstance(v, u.Quantity) for v in vals)

    def test_sma_multiple_orbits(self, all_orbits):
        """Equinoctial sma matches Keplerian sma across a range of orbit types."""
        for d in all_orbits:
            eq_sma  = tell.elementval(_make_eq_from_pvt(d['pvt']), 'sma')
            kep_sma = tell.elementval(d['kep'], 'sma')
            np.testing.assert_allclose(eq_sma.si.value, kep_sma.si.value,
                                       rtol=1e-9,
                                       err_msg=f"sma mismatch for {d}")


# -----------------------------------------------------------------------
# Equinoctial: round-trip tests
# -----------------------------------------------------------------------

class TestEquinoctialRoundTrip:

    def _pvt_array(self, pvtobj):
        return np.concatenate([pvtobj.cartesian['position'].si.value,
                               pvtobj.cartesian['velocity'].si.value])

    def test_pvt_to_equinoctial_to_pvt(self, leo1):
        """PVT → equinoctial ElementSetT → PVT should recover original."""
        eq_est = _make_eq_from_pvt(leo1['pvt'])
        recovered = eq_est.pvt()
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(leo1['pvt']),
                                   rtol=1e-9)

    def test_round_trip_elliptical(self, ell2):
        eq_est = _make_eq_from_pvt(ell2['pvt'])
        recovered = tell.pvt(eq_est)
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(ell2['pvt']),
                                   rtol=1e-9)

    def test_round_trip_geo(self, geo1):
        eq_est = _make_eq_from_pvt(geo1['pvt'])
        recovered = tell.pvt(eq_est)
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(geo1['pvt']),
                                   rtol=1e-9)


# -----------------------------------------------------------------------
# Circular: structural tests
# -----------------------------------------------------------------------

class TestCircularStructure:

    def test_returns_elementset_t(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert isinstance(est, ElementSetT)

    def test_has_t_attribute(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert isinstance(est.time, Time)

    def test_els_is_quantity(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert isinstance(est.elements, u.Quantity)

    def test_field_names_mean(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert set(circeltma_names) <= set(est.elements.dtype.names)

    def test_unpack(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        els, t = est
        assert isinstance(els, u.Quantity)
        assert isinstance(t, Time)

    def test_repr(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert 'ElementSetT' in repr(est)

    def test_equality(self, leo1):
        est1 = _make_circ_from_pvt(leo1['pvt'])
        est2 = _make_circ_from_pvt(leo1['pvt'])
        assert est1 == est2

    def test_inequality(self, leo1, leo2):
        c1 = _make_circ_from_pvt(leo1['pvt'])
        c2 = _make_circ_from_pvt(leo2['pvt'])
        assert c1 != c2


# -----------------------------------------------------------------------
# Circular: predicate tests
# -----------------------------------------------------------------------

class TestIscircels:

    def test_true_with_epoch(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert iscircels(est, est=True)

    def test_true_els_only(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert iscircels(est.elements, est=False)

    def test_false_for_keplerian_est(self, leo1):
        assert not iscircels(leo1['kep'])

    def test_false_for_plain_quantity(self):
        assert not iscircels(1.0 * u.km, est=False)

    def test_false_for_equinoctial_est(self, leo1):
        eq = _make_eq_from_pvt(leo1['pvt'])
        assert not iscircels(eq)

    def test_iskepels_false_for_circular(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert not iskepels(est)


# -----------------------------------------------------------------------
# Circular: numerical consistency
# -----------------------------------------------------------------------

class TestCircularValues:

    def test_sma_matches_keplerian(self, leo1):
        circ_sma = tell.elementval(_make_circ_from_pvt(leo1['pvt']), 'sma')
        kep_sma  = tell.elementval(leo1['kep'], 'sma')
        np.testing.assert_allclose(circ_sma.si.value, kep_sma.si.value, rtol=1e-10)

    def test_inc_matches_keplerian(self, leo2):
        circ_inc = tell.elementval(_make_circ_from_pvt(leo2['pvt']), 'inc')
        kep_inc  = tell.elementval(leo2['kep'], 'inc')
        np.testing.assert_allclose(circ_inc.si.value, kep_inc.si.value, rtol=1e-10)

    def test_raan_matches_keplerian(self, leo2):
        circ_raan = tell.elementval(_make_circ_from_pvt(leo2['pvt']), 'raan')
        kep_raan  = tell.elementval(leo2['kep'], 'raan')
        np.testing.assert_allclose(circ_raan.si.value, kep_raan.si.value, rtol=1e-10)

    def test_circular_orbit_cex_cey_near_zero(self, leo1):
        """For a near-circular orbit, cex and cey ≈ 0."""
        est = _make_circ_from_pvt(leo1['pvt'])
        np.testing.assert_allclose(est.elements['cex'].value, 0.0, atol=1e-6)
        np.testing.assert_allclose(est.elements['cey'].value, 0.0, atol=1e-6)

    def test_elementval_list_from_circular(self, leo2):
        est  = _make_circ_from_pvt(leo2['pvt'])
        vals = tell.elementval(est, ['sma', 'cex', 'cey', 'inc', 'raan', 'mla'])
        assert len(vals) == 6
        assert all(isinstance(v, u.Quantity) for v in vals)

    def test_sma_multiple_orbits(self, all_orbits):
        for d in all_orbits:
            circ_sma = tell.elementval(_make_circ_from_pvt(d['pvt']), 'sma')
            kep_sma  = tell.elementval(d['kep'], 'sma')
            np.testing.assert_allclose(circ_sma.si.value, kep_sma.si.value,
                                       rtol=1e-9,
                                       err_msg=f"sma mismatch for {d}")


# -----------------------------------------------------------------------
# Circular: round-trip tests
# -----------------------------------------------------------------------

class TestCircularRoundTrip:

    def _pvt_array(self, pvtobj):
        return np.concatenate([pvtobj.cartesian['position'].si.value,
                               pvtobj.cartesian['velocity'].si.value])

    def test_pvt_to_circular_to_pvt(self, leo1):
        circ_est  = _make_circ_from_pvt(leo1['pvt'])
        recovered = circ_est.pvt()
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(leo1['pvt']),
                                   rtol=1e-9)

    def test_round_trip_elliptical(self, ell2):
        circ_est  = _make_circ_from_pvt(ell2['pvt'])
        recovered = tell.pvt(circ_est)
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(ell2['pvt']),
                                   rtol=1e-9)

    def test_round_trip_gps(self, gps1):
        circ_est  = _make_circ_from_pvt(gps1['pvt'])
        recovered = tell.pvt(circ_est)
        np.testing.assert_allclose(self._pvt_array(recovered),
                                   self._pvt_array(gps1['pvt']),
                                   rtol=1e-9)


# -----------------------------------------------------------------------
# Cross-type predicate exclusivity
# -----------------------------------------------------------------------

class TestPredicateExclusivity:
    """Each is*els predicate must be True for exactly one element type."""

    def test_keplerian_only(self, leo1):
        est = leo1['kep']
        assert     iskepels(est)
        assert not isequels(est)
        assert not iscircels(est)

    def test_equinoctial_only(self, leo1):
        est = _make_eq_from_pvt(leo1['pvt'])
        assert not iskepels(est)
        assert     isequels(est)
        assert not iscircels(est)

    def test_circular_only(self, leo1):
        est = _make_circ_from_pvt(leo1['pvt'])
        assert not iskepels(est)
        assert not isequels(est)
        assert     iscircels(est)
