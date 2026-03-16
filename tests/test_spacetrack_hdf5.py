"""
Tests for spacetrack_hdf5.py serialization of MeanElementSetT.

All tests construct MeanElementSetT objects entirely in-memory using
hard-coded TLE data, so they run without any access to space-track.org.

An optional fixture-based test loads a pre-recorded HDF5 file produced by
make_spacetrack_fixtures.py (if it exists), providing a higher-fidelity
integration check.

MeanElementSetT is a plain class (not a namedtuple) with fields:
    els, t, tle, model, scdata
and __eq__ based on tle, model, and scdata only.
"""
import os
import pathlib
import tempfile

import pytest
import numpy as np
import astropy.units as u
import astropy.time

from fsc.hdf5_io import save, load
import astropy_hdf5io        # registers Quantity, Time, etc.
import tellurion as tell

# ---------------------------------------------------------------------------
# Pre-recorded fixture (optional – only used if the file exists)
# ---------------------------------------------------------------------------
_DATA_DIR = pathlib.Path(__file__).parent / 'data'
_FIXTURE_FILE = _DATA_DIR / 'spacetrack_fixtures.h5'

# ---------------------------------------------------------------------------
# Hard-coded reference data – Sentinel-3A TLE from 2024-11-01 (public domain)
# ---------------------------------------------------------------------------
_TLE_LINE1 = '1 41335U 16011A   24306.50000000  .00000072  00000-0  45678-4 0  9993'
_TLE_LINE2 = '2 41335  98.5693 312.3456 0001234  90.1234 270.0123 14.26734567456789'

_EPOCH = astropy.time.Time('2024-11-01T12:00:00', scale='utc')

_ELS = {
    'sma':    7167.744 * u.km,
    'ecc':    0.0001234 * u.dimensionless_unscaled,
    'inc':    98.5693 * u.deg,
    'raan':   312.3456 * u.deg,
    'argper': 90.1234 * u.deg,
    'ma':     270.0123 * u.deg,
    'memo':   14.26734567 * u.rev / u.day,
    'memod':  7.2e-8 * u.rev / (u.day * u.day),
    'memodd': 0.0 * u.rev / (u.day * u.day * u.day),
    'period': 100.93 * u.min,
    'peralt': 789.5 * u.km,
    'apoalt': 792.1 * u.km,
    'B':      0.000456 * u.m**2 / u.kg,
}

_SCDATA = {
    'name':    'SENTINEL 3A',
    'type':    'PAYLOAD',
    'catid':   41335,
    'intldes': '2016-011A',
}

_TLE = (_TLE_LINE1, _TLE_LINE2)


def _make_sentinel3a() -> tell.MeanElementSetT:
    """Construct a MeanElementSetT for Sentinel-3A entirely in memory."""
    return tell.MeanElementSetT(
        els=_ELS,
        t=_EPOCH,
        tle=_TLE,
        model='SGP4',
        scdata=_SCDATA,
    )


def _make_iss() -> tell.MeanElementSetT:
    """Construct a MeanElementSetT for ISS entirely in memory."""
    tle1 = '1 25544U 98067A   24306.50000000  .00001234  00000-0  23456-4 0  9991'
    tle2 = '2 25544  51.6436 150.1234 0001234  90.1234 270.0123 15.49567890123456'
    return tell.MeanElementSetT(
        els={**_ELS, 'inc': 51.6436 * u.deg},
        t=astropy.time.Time('2024-11-02T00:00:00', scale='utc'),
        tle=(tle1, tle2),
        model='SGP4',
        scdata={**_SCDATA, 'name': 'ISS', 'catid': 25544, 'intldes': '1998-067A'},
    )


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
def _roundtrip(obj: tell.MeanElementSetT) -> tell.MeanElementSetT:
    """Save obj to a temp HDF5 file and reload it."""
    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
        fname = f.name
    try:
        save(obj, fname)
        return load(fname)
    finally:
        os.unlink(fname)


# ---------------------------------------------------------------------------
# Tests: class construction
# ---------------------------------------------------------------------------
class TestMeanElementSetTConstruction:
    """Verify the plain-class constructor behaves as expected."""

    def test_fields_accessible(self):
        """All five fields are accessible as attributes."""
        obj = _make_sentinel3a()
        assert obj.els is _ELS
        assert obj.t is _EPOCH
        assert obj.tle == _TLE
        assert obj.model == 'SGP4'
        assert obj.scdata is _SCDATA

    def test_repr(self):
        """repr includes name, time, and model."""
        obj = _make_sentinel3a()
        r = repr(obj)
        assert 'SENTINEL 3A' in r
        assert 'SGP4' in r

    def test_equality_based_on_tle_model_scdata(self):
        """__eq__ uses tle, model, scdata — not els or t."""
        obj1 = _make_sentinel3a()
        # Same tle/model/scdata, different els and t
        obj2 = tell.MeanElementSetT(
            els={'sma': 9999.0 * u.km},          # different
            t=astropy.time.Time('2020-01-01'),     # different
            tle=_TLE,
            model='SGP4',
            scdata=_SCDATA,
        )
        assert obj1 == obj2

    def test_inequality_different_tle(self):
        """Objects with different TLEs are not equal."""
        obj1 = _make_sentinel3a()
        obj2 = _make_iss()
        assert obj1 != obj2

    def test_inequality_different_model(self):
        """Objects with different models are not equal."""
        obj1 = _make_sentinel3a()
        obj2 = tell.MeanElementSetT(
            els=_ELS, t=_EPOCH, tle=_TLE,
            model='J2',          # different
            scdata=_SCDATA,
        )
        assert obj1 != obj2

    def test_inequality_different_scdata(self):
        """Objects with different scdata are not equal."""
        obj1 = _make_sentinel3a()
        obj2 = tell.MeanElementSetT(
            els=_ELS, t=_EPOCH, tle=_TLE, model='SGP4',
            scdata={**_SCDATA, 'catid': 99999},   # different
        )
        assert obj1 != obj2

    def test_not_equal_to_other_types(self):
        """Comparison with non-MeanElementSetT returns NotImplemented."""
        obj = _make_sentinel3a()
        assert obj.__eq__("not a MeanElementSetT") is NotImplemented


# ---------------------------------------------------------------------------
# Tests: HDF5 round-trip
# ---------------------------------------------------------------------------
class TestMeanElementSetTRoundtrip:
    """Round-trip serialization tests for MeanElementSetT."""

    def test_type_preserved(self):
        """Loaded object is a MeanElementSetT."""
        loaded = _roundtrip(_make_sentinel3a())
        assert isinstance(loaded, tell.MeanElementSetT)

    def test_equality_after_roundtrip(self):
        """Round-tripped object compares equal to original (tle/model/scdata)."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert loaded == original

    def test_model_preserved(self):
        """Propagation model string is preserved."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert loaded.model == original.model

    def test_tle_line1_preserved(self):
        """TLE line 1 is preserved exactly."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert loaded.tle[0] == _TLE_LINE1

    def test_tle_line2_preserved(self):
        """TLE line 2 is preserved exactly."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert loaded.tle[1] == _TLE_LINE2

    def test_tle_is_two_tuple(self):
        """TLE round-trips as a two-element sequence."""
        loaded = _roundtrip(_make_sentinel3a())
        assert len(loaded.tle) == 2

    def test_epoch_preserved(self):
        """Epoch Time is preserved to sub-second precision."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert isinstance(loaded.t, astropy.time.Time)
        dt_s = abs((loaded.t - original.t).to(u.s).value)
        assert dt_s < 1e-6

    def test_scdata_name_preserved(self):
        """scdata 'name' string field is preserved."""
        loaded = _roundtrip(_make_sentinel3a())
        assert loaded.scdata['name'] == _SCDATA['name']

    def test_scdata_type_preserved(self):
        """scdata 'type' string field is preserved."""
        loaded = _roundtrip(_make_sentinel3a())
        assert loaded.scdata['type'] == _SCDATA['type']

    def test_scdata_intldes_preserved(self):
        """scdata 'intldes' string field is preserved."""
        loaded = _roundtrip(_make_sentinel3a())
        assert loaded.scdata['intldes'] == _SCDATA['intldes']

    def test_scdata_catid_value_preserved(self):
        """scdata 'catid' integer value is preserved."""
        loaded = _roundtrip(_make_sentinel3a())
        assert loaded.scdata['catid'] == _SCDATA['catid']

    def test_scdata_catid_is_int(self):
        """catid is explicitly coerced to int on load (not numpy int64 etc.)."""
        loaded = _roundtrip(_make_sentinel3a())
        assert isinstance(loaded.scdata['catid'], int)

    def test_all_els_keys_present(self):
        """All expected orbital element keys survive the round-trip."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        assert set(loaded.els.keys()) == set(original.els.keys())

    def test_orbital_elements_values(self):
        """Each orbital element's numerical value is preserved."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        for key, orig_q in original.els.items():
            loaded_q = loaded.els[key]
            assert np.isclose(
                loaded_q.to(orig_q.unit).value,
                orig_q.value,
                rtol=1e-12,
            ), f"els['{key}'] mismatch: {loaded_q} != {orig_q}"

    def test_orbital_element_units_compatible(self):
        """Each orbital element can be converted to the original unit."""
        original = _make_sentinel3a()
        loaded = _roundtrip(original)
        for key, orig_q in original.els.items():
            loaded_q = loaded.els[key]
            # Raises UnitConversionError if physically incompatible
            loaded_q.to(orig_q.unit)

    def test_els_values_are_quantities(self):
        """All orbital elements are astropy Quantities after round-trip."""
        loaded = _roundtrip(_make_sentinel3a())
        for key, val in loaded.els.items():
            assert isinstance(val, u.Quantity), \
                f"els['{key}'] is {type(val)}, expected Quantity"


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------
class TestMeanElementSetTEdgeCases:
    """Edge cases and boundary conditions."""

    def test_zero_valued_memodd_preserved(self):
        """Zero-valued memodd (mean motion double-dot) is preserved."""
        loaded = _roundtrip(_make_sentinel3a())
        assert loaded.els['memodd'].value == pytest.approx(0.0)

    def test_dimensionless_eccentricity(self):
        """Dimensionless eccentricity survives unit handling."""
        loaded = _roundtrip(_make_sentinel3a())
        ecc = loaded.els['ecc']
        assert ecc.unit.is_equivalent(u.dimensionless_unscaled)
        assert np.isclose(ecc.value, _ELS['ecc'].value, rtol=1e-12)

    def test_tle_whitespace_preserved(self):
        """TLE lines with internal spaces are stored and retrieved exactly."""
        loaded = _roundtrip(_make_sentinel3a())
        # Standard TLE lines contain fixed-width space-separated fields
        assert loaded.tle[0] == _TLE_LINE1
        assert loaded.tle[1] == _TLE_LINE2

    def test_two_independent_objects(self):
        """Two different objects round-trip independently."""
        loaded1 = _roundtrip(_make_sentinel3a())
        loaded2 = _roundtrip(_make_iss())

        assert loaded1.scdata['catid'] == 41335
        assert loaded2.scdata['catid'] == 25544
        assert loaded1.scdata['name'] == 'SENTINEL 3A'
        assert loaded2.scdata['name'] == 'ISS'
        assert not np.isclose(
            loaded1.els['inc'].to(u.deg).value,
            loaded2.els['inc'].to(u.deg).value,
        )

    def test_double_roundtrip(self):
        """A second round-trip is identical to the first."""
        once = _roundtrip(_make_sentinel3a())
        twice = _roundtrip(once)
        assert twice == once
        assert twice.scdata['catid'] == 41335
        dt_s = abs((twice.t - once.t).to(u.s).value)
        assert dt_s < 1e-6


# ---------------------------------------------------------------------------
# Optional fixture-based integration test
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not _FIXTURE_FILE.exists(),
    reason=(
        f"Pre-recorded fixture not found at {_FIXTURE_FILE}. "
        "Run tests/make_spacetrack_fixtures.py once to create it."
    ),
)
class TestSpacetrackHdf5Fixture:
    """
    Integration tests using a pre-recorded HDF5 file.

    These tests only run when tests/data/spacetrack_fixtures.h5 exists
    (created by make_spacetrack_fixtures.py with real space-track.org data).
    """

    def test_fixture_loads_as_dict(self):
        """Fixture file loads as a dict of MeanElementSetT objects."""
        data = load(str(_FIXTURE_FILE))
        assert isinstance(data, dict)
        assert len(data) > 0

    def test_fixture_values_are_mean_element_sets(self):
        """Each value in the loaded dict is a MeanElementSetT."""
        data = load(str(_FIXTURE_FILE))
        for name, obj in data.items():
            assert isinstance(obj, tell.MeanElementSetT), \
                f"Expected MeanElementSetT for '{name}', got {type(obj)}"

    def test_fixture_sentinel3a_present(self):
        """Sentinel-3A is present in the fixture."""
        data = load(str(_FIXTURE_FILE))
        assert 'SENTINEL 3A' in data

    def test_fixture_sentinel3a_catid(self):
        """Sentinel-3A has the correct NORAD catalog ID and it is an int."""
        data = load(str(_FIXTURE_FILE))
        sentinel = data['SENTINEL 3A']
        assert sentinel.scdata['catid'] == 41335
        assert isinstance(sentinel.scdata['catid'], int)

    def test_fixture_sentinel3a_model(self):
        """Sentinel-3A propagation model is SGP4."""
        data = load(str(_FIXTURE_FILE))
        assert data['SENTINEL 3A'].model == 'SGP4'

    def test_fixture_epoch_is_time(self):
        """Fixture epoch is an astropy.time.Time object for every satellite."""
        data = load(str(_FIXTURE_FILE))
        for name, obj in data.items():
            assert isinstance(obj.t, astropy.time.Time), \
                f"Epoch for '{name}' is not a Time object"

    def test_fixture_els_have_units(self):
        """All orbital elements in fixture have astropy Quantity units."""
        data = load(str(_FIXTURE_FILE))
        for sat_name, obj in data.items():
            for key, val in obj.els.items():
                assert isinstance(val, u.Quantity), \
                    f"els['{key}'] for '{sat_name}' is not a Quantity"

    def test_fixture_roundtrip(self):
        """Objects loaded from fixture survive a second round-trip unchanged."""
        data = load(str(_FIXTURE_FILE))
        for name, obj in data.items():
            reloaded = _roundtrip(obj)
            assert reloaded == obj        # uses tle/model/scdata equality
            assert reloaded.scdata['catid'] == obj.scdata['catid']


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
