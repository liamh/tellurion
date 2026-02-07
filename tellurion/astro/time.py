"""
AstroPy definitions
"""

from typing import Final
import warnings
import datetime
import bisect
import numpy as np
import astropy.units as u
import astropy.table
import astropy.time
import astropy.timeseries
import collections.abc

################################################################################
## Time
################################################################################

def abstime(ratimes, reftime='now'):
    """Convert `ratimes`, which is a relative time (also known as
    "time delta" or "time interval"), or an absolute time, or an
    iterable of those things, into an absolute time
    (astropy.time.Time) or a list of absolute times. If `reftime` is
    not provided, it defaults to the current time.

    Parameters
    ----------
    ratimes :  str, int, float, astropy.time.Time, u.Quantity, datetime, or iterable
        The relative or absolute time(s) to convert.
    reftime : str or astropy.time.Time, optional
        The reference time for relative times. Default is ``'now'``.

    Returns
    -------
    astropy.time.Time
        The absolute time(s).

    Examples
    --------
    >>> import tellurion as tell
    >>> import astropy.units as u
    >>> import numpy as np
    >>> newyear = tell.abstime('2025-01-01T00:00:00')
    >>> prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12)
    >>> tell.abstime(['2025-01-01T00:00:00', '2025-01-02T00:00:00'])
    >>> tell.abstime(prop5m1h, newyear)
    >>> tell.abstime(5*u.hour, newyear)
    >>> tell.abstime('12d 17hr 23min 33.1s', newyear)
    >>> tell.abstime([5*u.hour, '12d 17hr 23min 33.1s'], newyear)

    Using datetime objects:

    >>> import datetime
    >>> tell.abstime(datetime.datetime(2025, 4, 17, 22, 54, 8, 684006))
    """
    if reftime=='now':
        reftime = astropy.time.Time(datetime.datetime.now(datetime.UTC), scale='utc')
    if type(ratimes)==astropy.time.Time:
        return ratimes
    if type(ratimes) in [datetime.datetime, np.datetime64]:
        return astropy.time.Time(ratimes.isoformat())
    if type(ratimes)==str:
        try:
            return (astropy.time.Time(ratimes, scale='utc'))
        except:
            try:
                return abstime(astropy.time.TimeDelta(ratimes).to_value('sec')*u.s, reftime)
            except:
                raise ValueError("Cannot interpret string as relative or absolute time")
    if type(ratimes) in [int, float, np.float64]:
        return abstime(ratimes*u.s, reftime)
    if type(ratimes) == u.Quantity:
        return reftime + ratimes
    if isinstance(ratimes, collections.abc.Iterable):
        cum = abstime(ratimes[0], reftime)
        for t1 in ratimes[1:]:
            cum = time_concat(cum, abstime(t1, reftime))
        return cum

def time_concat(time1, time2):
    def tval(time):
        if time.isscalar:
            return [time.value]
        else:
            return time.value
    return astropy.time.Time(np.concatenate([tval(time1), tval(time2)]))

def timesec(t):
    """Convert a u.Quantity to seconds as a Python float"""
    if type(t) is u.Quantity and u.get_physical_type(t) == 'time':
        pt = t.si.value.tolist() # convert to seconds and get the value_unit
    elif isinstance(t, collections.abc.Iterable):
        return [timesec(i) for i in t]
    else:
        pt = float(t) # assume seconds
    return(pt)

def tc(t):
    """Convert a u.Quantity with physical dimension time to a string
    of duration (time interval) components; this is the inverse of
    tq().

    Parameters
    ----------
    t : u.Quantity
        A quantity with time dimension.

    Returns
    -------
    str
        Duration string with components (days, hours, minutes, seconds).

    Examples
    --------
    >>> import tellurion as tell
    >>> import astropy.units as u
    >>> tell.tc(tell.tq('12d 17hr 23min 33.1s'))
    '12d 17hr 23min 33.1s'
    >>> tell. tc(123456*u.s)
    '1d 10hr 17min 36.0s'
    """
    return astropy.time.TimeDelta(t).quantity_str

def tq(compstr):
    """Convert a string of duration (time interval) components to a
    u.Quantity with physical dimension time; this is the inverse of
    tc().

    Parameters
    ----------
    compstr : str
        Duration string with components (e.g., '12d 17hr 23min 33.1s').

    Returns
    -------
    u.Quantity
        Time quantity in seconds.

    Examples
    --------
    >>> import tellurion as tell
    >>> import astropy. units as u
    >>> tell. tq('12d 17hr 23min 33.1s')
    <Quantity 1099413.1 s>
    >>> tell.tq(tell.tc(123456*u. s))
    <Quantity 123456. s>
    """
    return astropy.time.TimeDelta(compstr).to_value('sec')*u.s

################################################################################
## Time
################################################################################

prefnumabstime = 'mjd' # Preferred numerical format for absolute time
astropy.time.Time.to_array = lambda self, format=prefnumabstime: to_array(self, format)

def to_array(tms, format=prefnumabstime):
    """Convert the absolute time(s) to a Numpy array using the
    specified format to convert absolute times to a float.

    Parameters
    ----------
    tms :  astropy.time.Time
        The time or times to convert.
    format : str, optional
        The time format for conversion. Common formats include ``'mjd'``
        (Modified Julian Date), ``'jd'`` (Julian Date), ``'unix'`` (Unix
        timestamp), ``'cxcsec'`` (Chandra X-ray Center seconds), ``'gps'``
        (GPS seconds), ``'plot_date'`` (Matplotlib plot date). See the
        `astropy.time.Time formats documentation
        <https://docs.astropy.org/en/stable/time/index.html#time-format>`_
        for a complete list.  Default is ``'mjd'``.

    Returns
    -------
    numpy.ndarray
        Array of time values in the specified format.

    Examples
    --------
    >>> import tellurion as tell
    >>> import astropy.time
    >>> t = astropy.time.Time('2025-01-01T00:00:00')
    >>> tell.to_array(t, format='mjd')
    array([60310.])
    >>> tell.to_array(t, format='jd')
    array([2460310.5])
    """
    if tms.shape == ():
        mjds = np.array([tms.to_value(format=format)])
    else:
        mjds = tms.to_value(format)
    return mjds

def from_array(array, format=prefnumabstime):
    """Convert the absolute time(s) from a Numpy array using the
    specified format to convert absolute times from a float.

    Parameters
    ----------
    tms :  astropy.time.Time
        The time or times to convert.
    format : str, optional
        The time format for conversion. Common formats include ``'mjd'``
        (Modified Julian Date), ``'jd'`` (Julian Date), ``'unix'`` (Unix
        timestamp), ``'cxcsec'`` (Chandra X-ray Center seconds), ``'gps'``
        (GPS seconds), ``'plot_date'`` (Matplotlib plot date). See the
        `astropy.time.Time formats documentation
        <https://docs.astropy.org/en/stable/time/index.html#time-format>`_
        for a complete list.  Default is ``'mjd'``.

    Returns
    -------
    astropy.time.Time
        The time(s)

    Examples
    --------
    >>> import tellurion as tell
    >>> import numpy as np
    >>> tell.from_array(np.array([60310.]), format='mjd')
    <Time object: scale='utc' format='isot' value=['2024-01-01T00:00:00.000']>
    >>> tell.from_array(np.array([2460310.5]), format='jd')
    <Time object: scale='utc' format='isot' value=['2024-01-01T00:00:00.000']>
    """
    return astropy.time.Time(astropy.time.Time(array, format=format).to_value('isot'))

################################################################################
## Time series and Tables
################################################################################

def striptime(ts):
    """Remove the time column from a TimeSeries and return a plain Table"""
    return astropy.table.Table([ts[k] for k in ts.keys()[1:None]])

def hcat(ts1, ts2):
    """Concatenate timeseries by adding columns from another time series

    Warns
    -----
    If times are not all equal in both timeseries
    """
    if not(all(ts1['time'].__eq__(ts2['time']))):
        warnings.warn("Times are not all equal; using times from first set")
    return astropy.table.hstack([ts1, striptime(ts2)])

def fromtime(ts, reftime='now', label = 'from now', copy = True):
    """The time series `ts` starting at the specified reference time `reftime` (default is the current time) and new column showing the elapsed time from the reference time. Make a new series if `copy` is ``True`` (the default); otherwise, modify the original time series."""
    if reftime=='now':
        reftime = abstime(0)
    elif reftime=='epoch':
        reftime=ts.time[0]
        label='from epoch'
    if copy:
        newts = ts.copy()
    else:
        newts = ts
    rowstart = bisect.bisect_left(newts.time, reftime)
    newts.remove_rows(slice(0, rowstart))
    newcol = astropy.table.Column((newts.time - reftime).quantity_str, name = label)
    if label in newts.keys():
        newts.remove_column(label)
    newts.add_column(newcol, index=1)
    return newts

def _fromtime_method(self, reltime):
    """Extract a portion of the time series relative to its start or end time.

    If `reltime` > 0, start at that time past the start time of the table.
    If `reltime` < 0, start at abs(`reltime`) before the end time of the table.

    Parameters
    ----------
    reltime : u.Quantity
        Relative time offset. Positive values are relative to the start,
        negative values are relative to the end.

    Returns
    -------
    astropy.timeseries.TimeSeries
        Time series starting from the specified relative time with an
        elapsed time column.

    Examples
    --------
    Get time series starting 1 day from the beginning:

    >>> import astropy.units as u
    >>> simorb.ephemeris().fromtime(1*u.day)

    Get time series starting 1 hour before the end:

    >>> simorb.ephemeris().fromtime(-1*u.h)
    """
    if reltime > 0:
        return fromtime(self, abstime(reltime, self.time[0]),
                       'from start+' + reltime.to_string())
    else:
        return fromtime(self, abstime(reltime, self.time[-1]),
                       'from end' + reltime.to_string())

# Attach the method to the class
astropy.timeseries.TimeSeries.fromtime = _fromtime_method

def _components_method(self, column, remove=True):
    """Replace a 3-vector column with three scalar component columns.

    Splits a column containing 3-vectors into three separate columns
    with ' x', ' y', and ' z' suffixes so that all components may be
    viewed individually.

    Parameters
    ----------
    column : str
        Name of the column containing 3-vectors to split.
    remove : bool, optional
        If ``True``, remove the original vector column after splitting.
        Default is ``True``.

    Examples
    --------
    >>> lasthour = simorb.ephemeris().fromtime(-1*u.h)
    >>> lasthour.components('position')
    >>> lasthour.components('velocity')
    """
    # Split the 3-vector into three separate columns
    self[column[0:3]+' x'] = self[column][:, 0]
    self[column[0:3]+' y'] = self[column][:, 1]
    self[column[0:3]+' z'] = self[column][:, 2]
    # Optionally remove the original column
    if remove:
        self.remove_column(column)

astropy.timeseries.TimeSeries.components = _components_method
