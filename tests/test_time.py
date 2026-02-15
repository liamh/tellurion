"""
Comprehensive test suite for time module.

Tests time handling with AstroPy-style conventions.
"""
import pytest
import numpy as np
import astropy.units as u
import astropy.time
import astropy.timeseries
import datetime
from numpy.testing import assert_allclose, assert_array_equal

# Import the module to test
from tellurion.astro import time as atime


class TestAbstime:
    """Tests for abstime() function."""

    def test_abstime_from_string_iso(self):
        """Test creating absolute time from ISO string."""
        t = atime.abstime('2025-01-01T00:00:00')

        assert isinstance(t, astropy.time.Time)
        assert t.iso == '2025-01-01 00:00:00.000'

    def test_abstime_from_time_object(self):
        """Test that Time objects pass through unchanged."""
        orig = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime(orig)

        assert t is orig

    def test_abstime_from_datetime(self):
        """Test creating from Python datetime."""
        dt = datetime.datetime(2025, 1, 1, 12, 30, 45)
        t = atime.abstime(dt)

        assert isinstance(t, astropy.time.Time)
        assert t.datetime.year == 2025
        assert t.datetime.month == 1
        assert t.datetime.day == 1

    def test_abstime_from_numpy_datetime64(self):
        """Test creating from numpy datetime64."""
        dt = np.datetime64('2025-01-01T12:30:45')
        t = atime.abstime(dt)

        assert isinstance(t, astropy.time.Time)

    def test_abstime_relative_quantity(self):
        """Test relative time from quantity."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime(5*u.hour, reftime)

        expected = astropy.time.Time('2025-01-01T05:00:00')
        assert_allclose(t.mjd, expected.mjd, rtol=1e-10)

    def test_abstime_relative_seconds_int(self):
        """Test relative time from integer seconds."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime(3600, reftime)

        expected = astropy.time.Time('2025-01-01T01:00:00')
        assert_allclose(t.mjd, expected.mjd, rtol=1e-10)

    def test_abstime_relative_seconds_float(self):
        """Test relative time from float seconds."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime(3600.5, reftime)

        expected = reftime + 3600.5*u.s
        assert_allclose(t.mjd, expected.mjd, rtol=1e-10)

    def test_abstime_relative_string(self):
        """Test relative time from duration string."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime('12d 17hr 23min 33.1s', reftime)

        # 12 days + 17 hours + 23 minutes + 33.1 seconds
        total_seconds = 12*86400 + 17*3600 + 23*60 + 33.1
        expected = reftime + total_seconds*u.s
        assert_allclose(t.mjd, expected.mjd, rtol=1e-10)

    def test_abstime_list_of_strings(self):
        """Test creating from list of ISO strings."""
        t = atime.abstime(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])

        assert isinstance(t, astropy.time.Time)
        assert len(t) == 2
        assert t[0].iso == '2025-01-01 00:00:00.000'
        assert t[1].iso == '2025-01-02 00:00:00.000'

    def test_abstime_list_of_quantities(self):
        """Test creating from list of quantities (relative times)."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime([5*u.hour, 10*u.hour, 15*u.hour], reftime)

        assert isinstance(t, astropy.time.Time)
        assert len(t) == 3
        assert_allclose(t[0].mjd, (reftime + 5*u.hour).mjd, rtol=1e-10)
        assert_allclose(t[1].mjd, (reftime + 10*u.hour).mjd, rtol=1e-10)
        assert_allclose(t[2].mjd, (reftime + 15*u.hour).mjd, rtol=1e-10)

    def test_abstime_list_of_time_objects(self):
        """Test creating from list of Time objects."""
        t1 = astropy.time.Time('2025-01-01T00:00:00')
        t2 = astropy.time.Time('2025-01-02T00:00:00')
        t = atime.abstime([t1, t2])

        assert isinstance(t, astropy.time.Time)
        assert len(t) == 2
        assert_allclose(t[0].mjd, t1.mjd, rtol=1e-10)
        assert_allclose(t[1].mjd, t2.mjd, rtol=1e-10)

    def test_abstime_mixed_list(self):
        """Test creating from mixed list of types."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        t = atime.abstime([5*u.hour, '12d 17hr 23min 33.1s'], reftime)

        assert isinstance(t, astropy.time.Time)
        assert len(t) == 2

    def test_abstime_numpy_array_of_quantities(self):
        """Test creating from numpy array of quantities."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        intervals = np.linspace(5.0, 60.0, 12) * u.minute
        t = atime.abstime(intervals, reftime)

        assert isinstance(t, astropy.time.Time)
        assert len(t) == 12

    def test_abstime_default_reftime_now(self):
        """Test that default reftime is 'now'."""
        before = astropy.time.Time.now()
        t = atime.abstime(0)  # 0 seconds from now
        after = astropy.time.Time.now()

        # Should be between before and after
        assert before.mjd <= t.mjd <= after.mjd

    def test_abstime_invalid_string(self):
        """Test that invalid string raises ValueError."""
        with pytest.raises(ValueError, match="Cannot interpret string"):
            atime.abstime('invalid time string')


class TestTimeConcat:
    """Tests for time_concat() function."""

    def test_concat_two_scalar_times(self):
        """Test concatenating two scalar Time objects."""
        t1 = astropy.time.Time('2025-01-01T00:00:00')
        t2 = astropy.time.Time('2025-01-02T00:00:00')

        result = atime.time_concat(t1, t2)

        assert len(result) == 2
        assert_allclose(result[0].mjd, t1.mjd, rtol=1e-10)
        assert_allclose(result[1].mjd, t2.mjd, rtol=1e-10)

    def test_concat_array_and_scalar(self):
        """Test concatenating array and scalar."""
        t1 = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        t2 = astropy.time.Time('2025-01-03T00:00:00')

        result = atime.time_concat(t1, t2)

        assert len(result) == 3

    def test_concat_scalar_and_array(self):
        """Test concatenating scalar and array."""
        t1 = astropy.time.Time('2025-01-01T00:00:00')
        t2 = astropy.time.Time(['2025-01-02T00:00:00', '2025-01-03T00:00:00'])

        result = atime.time_concat(t1, t2)

        assert len(result) == 3

    def test_concat_two_arrays(self):
        """Test concatenating two arrays."""
        t1 = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        t2 = astropy.time.Time(['2025-01-03T00:00:00', '2025-01-04T00:00:00'])

        result = atime.time_concat(t1, t2)

        assert len(result) == 4


class TestTimesec:
    """Tests for timesec() function."""

    def test_timesec_from_quantity(self):
        """Test converting Quantity to seconds."""
        t = 5 * u.hour
        result = atime.timesec(t)

        assert result == 5 * 3600

    def test_timesec_from_seconds_quantity(self):
        """Test converting seconds Quantity."""
        t = 123.5 * u.s
        result = atime.timesec(t)

        assert_allclose(result, 123.5)

    def test_timesec_from_float(self):
        """Test converting plain float (assumed seconds)."""
        result = atime.timesec(100.5)

        assert result == 100.5

    def test_timesec_from_list(self):
        """Test converting list of quantities."""
        times = [1*u.hour, 2*u.hour, 3*u.hour]
        result = atime.timesec(times)

        assert result == [3600, 7200, 10800]


class TestTcTq:
    """Tests for tc() and tq() functions (time component conversion)."""

    def test_tc_from_quantity(self):
        """Test converting Quantity to component string."""
        t = 123456 * u.s
        result = atime.tc(t)

        assert result == '1d 10hr 17min 36.0s'

    def test_tq_from_string(self):
        """Test converting component string to Quantity."""
        result = atime.tq('12d 17hr 23min 33.1s')

        expected_seconds = 12*86400 + 17*3600 + 23*60 + 33.1
        assert_allclose(result.to(u.s).value, expected_seconds, rtol=1e-10)

    def test_tc_tq_roundtrip(self):
        """Test round-trip conversion."""
        original = 123456 * u.s
        string = atime.tc(original)
        result = atime.tq(string)

        assert_allclose(result.to(u.s).value, original.value, rtol=1e-10)

    def test_tq_tc_roundtrip(self):
        """Test reverse round-trip conversion."""
        original = '12d 17hr 23min 33.1s'
        quantity = atime.tq(original)
        result = atime.tc(quantity)

        assert result == original


class TestToArrayFromArray:
    """Tests for to_array() and from_array() functions."""

    def test_to_array_scalar(self):
        """Test converting scalar Time to array."""
        t = astropy.time.Time('2025-01-01T00:00:00')
        result = atime.to_array(t)

        assert isinstance(result, np.ndarray)
        assert result.shape == (1,)

    def test_to_array_multiple(self):
        """Test converting multiple times to array."""
        t = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        result = atime.to_array(t)

        assert isinstance(result, np.ndarray)
        assert result.shape == (2,)

    def test_to_array_mjd_format(self):
        """Test to_array with MJD format."""
        t = astropy.time.Time('2025-01-01T00:00:00')
        result = atime.to_array(t, format='mjd')

        assert_allclose(result[0], t.mjd, rtol=1e-10)

    def test_to_array_jd_format(self):
        """Test to_array with JD format."""
        t = astropy.time.Time('2025-01-01T00:00:00')
        result = atime.to_array(t, format='jd')

        assert_allclose(result[0], t.jd, rtol=1e-10)

    def test_from_array_scalar(self):
        """Test creating Time from scalar array."""
        mjd = np.array([60310.0])
        result = atime.from_array(mjd, format='mjd')

        assert isinstance(result, astropy.time.Time)

    def test_from_array_multiple(self):
        """Test creating Time from array."""
        mjd = np.array([60310.0, 60311.0, 60312.0])
        result = atime.from_array(mjd, format='mjd')

        assert isinstance(result, astropy.time.Time)
        assert len(result) == 3

    def test_to_array_from_array_roundtrip(self):
        """Test round-trip conversion."""
        original = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        array = atime.to_array(original)
        result = atime.from_array(array)

        assert_allclose(result.mjd, original.mjd, rtol=1e-10)

    def test_time_object_to_array_method(self):
        """Test Time.to_array() method added by module."""
        t = astropy.time.Time('2025-01-01T00:00:00')
        result = t.to_array()

        assert isinstance(result, np.ndarray)
        assert_allclose(result[0], t.mjd, rtol=1e-10)


class TestTimeSeriesFunctions:
    """Tests for TimeSeries helper functions."""

    def test_striptime(self):
        """Test removing time column from TimeSeries."""
        ts = astropy.timeseries.TimeSeries(
            time=astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00']),
            data={'value': [1, 2], 'flag': [True, False]}
        )

        result = atime.striptime(ts)

        assert 'time' not in result.colnames
        assert 'value' in result.colnames
        assert 'flag' in result.colnames

    def test_hcat(self):
        """Test horizontal concatenation of TimeSeries."""
        times = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])

        ts1 = astropy.timeseries.TimeSeries(
            time=times,
            data={'value1': [1, 2]}
        )

        ts2 = astropy.timeseries.TimeSeries(
            time=times,
            data={'value2': [3, 4]}
        )

        result = atime.hcat(ts1, ts2)

        assert 'value1' in result.colnames
        assert 'value2' in result.colnames
        assert len(result) == 2

    def test_hcat_warns_on_time_mismatch(self):
        """Test that hcat warns when times don't match."""
        ts1 = astropy.timeseries.TimeSeries(
            time=astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00']),
            data={'value1': [1, 2]}
        )

        ts2 = astropy.timeseries.TimeSeries(
            time=astropy.time.Time(['2025-01-01T00:00:00', '2025-01-03T00:00:00']),
            data={'value2': [3, 4]}
        )

        with pytest.warns(UserWarning, match="Times are not all equal"):
            result = atime.hcat(ts1, ts2)

    def test_fromtime_with_reftime(self):
        """Test fromtime with explicit reference time."""
        times = astropy.time.Time(['2025-01-01T00:00:00',
                                    '2025-01-01T01:00:00',
                                    '2025-01-01T02:00:00'])
        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'value': [1, 2, 3]}
        )

        reftime = astropy.time.Time('2025-01-01T00:00:00')
        result = atime.fromtime(ts, reftime=reftime, copy=True)

        assert 'from now' in result.colnames
        assert len(result) == 3

    def test_fromtime_with_epoch(self):
        """Test fromtime with epoch as reference."""
        times = astropy.time.Time(['2025-01-01T00:00:00',
                                    '2025-01-01T01:00:00',
                                    '2025-01-01T02:00:00'])
        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'value': [1, 2, 3]}
        )

        result = atime.fromtime(ts, reftime='epoch', copy=True)

        assert 'from epoch' in result.colnames

    def test_fromtime_removes_rows_before_reftime(self):
        """Test that fromtime removes rows before reference time."""
        times = astropy.time.Time(['2025-01-01T00:00:00',
                                    '2025-01-01T01:00:00',
                                    '2025-01-01T02:00:00',
                                    '2025-01-01T03:00:00'])
        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'value': [1, 2, 3, 4]}
        )

        reftime = astropy.time.Time('2025-01-01T01:30:00')
        result = atime.fromtime(ts, reftime=reftime, copy=True)

        # Should only have times >= reftime
        assert len(result) == 2  # Only last 2 times

    def test_timeseries_fromtime_method(self):
        """Test fromtime method added to TimeSeries."""
        times = astropy.time.Time(['2025-01-01T00:00:00',
                                    '2025-01-01T01:00:00',
                                    '2025-01-01T02:00:00'])
        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'value': [1, 2, 3]}
        )

        # Test positive offset (from start)
        result = ts.fromtime(1*u.hour)
        assert 'from start' in result.colnames[1]

    def test_timeseries_fromtime_method_negative(self):
        """Test fromtime method with negative offset (from end)."""
        times = astropy.time.Time(['2025-01-01T00:00:00',
                                    '2025-01-01T01:00:00',
                                    '2025-01-01T02:00:00'])
        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'value': [1, 2, 3]}
        )

        # Test negative offset (from end)
        result = ts.fromtime(-1*u.hour)
        assert 'from end' in result.colnames[1]

    def test_timeseries_components_method(self):
        """Test components method added to TimeSeries."""
        times = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        positions = np.array([[1, 2, 3], [4, 5, 6]]) * u.km

        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'position': positions}
        )

        ts.components('position', remove=True)

        assert 'pos x' in ts.colnames
        assert 'pos y' in ts.colnames
        assert 'pos z' in ts.colnames
        assert 'position' not in ts.colnames

    def test_timeseries_components_method_keep_original(self):
        """Test components method keeping original column."""
        times = astropy.time.Time(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
        velocities = np.array([[1, 2, 3], [4, 5, 6]]) * u.km / u.s

        ts = astropy.timeseries.TimeSeries(
            time=times,
            data={'velocity': velocities}
        )

        ts.components('velocity', remove=False)

        assert 'vel x' in ts.colnames
        assert 'vel y' in ts.colnames
        assert 'vel z' in ts.colnames
        assert 'velocity' in ts.colnames


class TestIntegration:
    """Integration tests combining multiple functions."""

    def test_create_time_series_from_relative_times(self):
        """Test creating ephemeris from relative times."""
        reftime = astropy.time.Time('2025-01-01T00:00:00')
        intervals = np.linspace(0, 60, 13) * u.minute
        times = atime.abstime(intervals, reftime)

        values = np.sin(np.linspace(0, 2*np.pi, 13))
        ts = astropy.timeseries.TimeSeries(time=times, data={'signal': values})

        assert len(ts) == 13
        assert ts.time[0].iso == '2025-01-01 00:00:00.000'

    def test_concatenate_time_objects_from_different_sources(self):
        """Test concatenating times created differently."""
        t1 = atime.abstime('2025-01-01T00:00:00')
        t2 = atime.abstime(86400, t1)  # 1 day later
        t3 = astropy.time.Time('2025-01-03T00:00:00')

        result = atime.abstime([t1, t2, t3])

        assert len(result) == 3
        assert_allclose(result[1].mjd - result[0].mjd, 1.0, rtol=1e-10)

def test_abstime():
    newyear = atime.abstime('2025-01-01T00:00:00')
    prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
    np.testing.assert_allclose(newyear.to_value('jd'), 2460676.5)
    np.testing.assert_allclose(atime.abstime(datetime.datetime(2025, 4, 17, 22, 54, 8, 684006)).to_value('jd'), \
                               2460783.454267176)
    np.testing.assert_allclose(atime.abstime(['2025-01-01T00:00:00', '2025-01-02T00:00:00', \
                                             '2025-01-03T00:00:00']).to_value('jd'), \
                               np.array([2460676.5, 2460677.5, 2460678.5]))
    np.testing.assert_allclose(atime.abstime(prop5m1h, newyear).to_value('jd'), \
                               np.array([2460676.50347222, 2460676.50694444, 2460676.51041667, \
                                         2460676.51388889, 2460676.51736111, 2460676.52083333, \
                                         2460676.52430556, 2460676.52777778, 2460676.53125, \
                                         2460676.53472222, 2460676.53819444, 2460676.54166667]))
    np.testing.assert_allclose(atime.abstime(5*u.hour, newyear).to_value('jd'), 2460676.7083333335)
    np.testing.assert_allclose(atime.abstime('12d 17hr 23min 33.1s', newyear).to_value('jd'), \
                               2460689.2246886576)
    np.testing.assert_allclose(atime.abstime([5*u.hour, '12d 17hr 23min 33.1s'], newyear).to_value('jd'), \
                               np.array([2460676.7083333335, 2460689.2246886576]))


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
