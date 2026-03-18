# HDF5 round-trip tests for ElementSetT.
#
# Requires fsc.hdf5-io and astropy-hdf5io to be installed.
# All tests are automatically skipped if either library is absent.
#
# Run with: pytest -q test_element_hdf5.py

import tempfile
import os

import pytest
import numpy as np
import astropy.units as u
from astropy.time import Time
import tellurion as tell
from tellurion.core.element import ElementSetT

# Importing element_hdf5 registers the ElementSetT serializer as a side-effect.
# The import is guarded so the whole module still collects cleanly even when
# the HDF5 libraries are absent; individual tests skip themselves as needed.
try:
    from fsc.hdf5_io import save, load
    from tellurion.core import element_hdf5  # noqa: F401 – registers serializer
    _HDF5IO_AVAILABLE = True
except ImportError:
    _HDF5IO_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not _HDF5IO_AVAILABLE,
    reason="fsc.hdf5-io / astropy-hdf5io not installed"
)


# -----------------------------------------------------------------------
# Shared orbit fixtures
# -----------------------------------------------------------------------

def _leo1_est():
    """ISS-like LEO circular orbit element set."""
    return tell.kepler(
        tell.allplane({"altper": 350*u.km, "altapo": 350*u.km,
                       "inc": 55.0*u.deg, "argper": 120.0*u.deg,
                       "raan": 20.0*u.deg, "ma": 30.0*u.deg}),
        tell.abstime('2026-01-01 05:55:00'))


def _geo1_est():
    """GEO orbit element set (uses memo instead of sma)."""
    return tell.kepler(
        tell.allplane({"memo": 1.0*u.rev/u.sday,
                       "ecc": 0.0,
                       "inc": 0.0*u.deg, "argper": 120.0*u.deg,
                       "raan": 0.0*u.deg, "ma": 0.0*u.deg}),
        tell.abstime('2026-01-01 20:30:00'))


def _vang1_est():
    """Vanguard-1-like highly elliptical orbit element set."""
    return tell.kepler(
        tell.allplane({"altper": 600.0*u.km, "altapo": 12000.0*u.km,
                       "inc": 36.0*u.deg, "argper": 140.0*u.deg,
                       "raan": 0.0*u.deg, "ma": 100.0*u.deg}),
        tell.abstime('2026-01-01 14:45:00'))


def _gps1_est():
    """GPS semisynchronous orbit element set (sma given directly)."""
    return tell.kepler(
        tell.allplane({"sma": tell.sma(2.0), "ecc": 0.0*u.dimensionless_unscaled,
                       "inc": 55.0*u.deg, "argper": 0.0*u.deg,
                       "raan": 120.0*u.deg, "ma": 77.0*u.deg}),
        tell.abstime('2026-01-01 12:20:00'))


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def _roundtrip(est):
    """Save *est* to a temp file, load it back, and return the loaded object."""
    with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
        fname = f.name
    try:
        save(est, fname)
        return load(fname)
    finally:
        os.unlink(fname)


def _assert_els_allclose(loaded, original, rtol=1e-10):
    """Assert all element fields match after a round-trip."""
    for name in original.els.dtype.names:
        np.testing.assert_allclose(
            loaded.els[name].si.value,
            original.els[name].si.value,
            rtol=rtol,
            err_msg=f"field '{name}' does not match after round-trip")


# -----------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------

class TestElementSetTHdf5:
    """HDF5 round-trip tests for ElementSetT."""

    def test_has_to_hdf5_method(self):
        """Importing element_hdf5 monkey-patches to_hdf5 onto ElementSetT."""
        est = _leo1_est()
        assert callable(getattr(est, 'to_hdf5', None))

    def test_roundtrip_returns_element_set_t(self):
        loaded = _roundtrip(_leo1_est())
        assert isinstance(loaded, ElementSetT)

    def test_roundtrip_epoch_preserved(self):
        est = _leo1_est()
        loaded = _roundtrip(est)
        assert loaded.t == est.t

    def test_roundtrip_leo1_all_fields(self):
        est = _leo1_est()
        _assert_els_allclose(_roundtrip(est), est)

    def test_roundtrip_geo1_sma(self):
        """GEO orbit uses memo; check the reconstructed sma is correct."""
        est = _geo1_est()
        loaded = _roundtrip(est)
        np.testing.assert_allclose(
            loaded.els['sma'].to(u.km).value,
            est.els['sma'].to(u.km).value,
            rtol=1e-10)

    def test_roundtrip_geo1_all_fields(self):
        est = _geo1_est()
        _assert_els_allclose(_roundtrip(est), est)

    def test_roundtrip_vang1_all_fields(self):
        """Vanguard-like eccentric orbit."""
        est = _vang1_est()
        _assert_els_allclose(_roundtrip(est), est)

    def test_roundtrip_gps1_all_fields(self):
        """GPS orbit with sma specified directly."""
        est = _gps1_est()
        _assert_els_allclose(_roundtrip(est), est)

    def test_roundtrip_field_names_preserved(self):
        """The structured dtype field names must survive the round-trip."""
        est = _leo1_est()
        loaded = _roundtrip(est)
        assert loaded.els.dtype.names == est.els.dtype.names

    def test_roundtrip_units_preserved(self):
        """Each field's unit must survive the round-trip."""
        est = _leo1_est()
        loaded = _roundtrip(est)
        for name in est.els.dtype.names:
            assert loaded.els[name].unit.is_equivalent(est.els[name].unit), \
                f"unit for field '{name}' changed: " \
                f"{loaded.els[name].unit} vs {est.els[name].unit}"

    def test_roundtrip_pvt_unchanged(self):
        """The PVT derived from the loaded element set should match the original."""
        est = _leo1_est()
        loaded = _roundtrip(est)
        pvt_orig   = tell.pvt(est)
        pvt_loaded = tell.pvt(loaded)
        np.testing.assert_allclose(
            pvt_loaded.cartesian['position'].si.value,
            pvt_orig.cartesian['position'].si.value,
            rtol=1e-10)
        np.testing.assert_allclose(
            pvt_loaded.cartesian['velocity'].si.value,
            pvt_orig.cartesian['velocity'].si.value,
            rtol=1e-10)