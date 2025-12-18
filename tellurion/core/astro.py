"""
AstroPy definitions
"""

import warnings
import datetime
import bisect
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
import astropy.table
import astropy.time
import collections.abc

################################################################################
## Units
################################################################################

#: A revolution (full circle), useful for two-line elements (mean motion in rev/day)
u.rev = u.revolution = u.def_unit('revolution', 2*np.pi*u.radian)
u.add_enabled_units(u.rev)

#: User's preferred units
prefunits = {"time": u.second, "length": u.km, "speed": u.km/u.second,
             "angle": u.degree, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}
prefunits["posvel"] = (prefunits["length"], prefunits["speed"])
#: Unit for posvel (m, m/s)
posvelsiu = u.StructuredUnit((u.meter, u.meter/u.second))
#: Units used by Orekit
siunits = {"time": u.second, "length": u.m, "speed": u.m/u.second,
           "angle": u.radian, "angular speed": u.radian/u.second,
           "dimensionless": u.dimensionless_unscaled}
orkunits = siunits

def abstime(ratimes, reftime='now'):
    '''Convert `ratimes`, which is a relative time (also known as
    "time delta" or "time interval"), or an absolute time, or an
    iterable of those things, into an absolute time
    (astropy.time.Time) or a list of absolute times. If `reftime` is
    not provided, it defaults to the current time.

    Examples
    newyear = tell.abstime('2025-01-01T00:00:00')
    prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
    tell.abstime(['2025-01-01T00:00:00', '2025-01-02T00:00:00', '2025-01-03T00:00:00']
    tell.abstime(prop5m1h, newyear)
    tell.abstime(5*u.hour, newyear)
    tell.abstime('12d 17hr 23min 33.1s', newyear)
    tell.abstime([5*u.hour, '12d 17hr 23min 33.1s'], newyear)
    import datetime
    tell.abstime(datetime.datetime(2025, 4, 17, 22, 54, 8, 684006))

    '''
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
    '''Convert a u.Quantity to seconds as a Python float'''
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

    tc(tq('12d 17hr 23min 33.1s'))
    '12d 17hr 23min 33.1s'
    tell.tc(123456*u.s)
    '1d 10hr 17min 36.0s'
    """
    return astropy.time.TimeDelta(t).quantity_str

def tq(compstr):
    """Convert a string of duration (time interval) components to a
    u.Quantity with physical dimension time; this is the inverse of
    tc().

    tq('12d 17hr 23min 33.1s')
    <Quantity 1099413.1 s>
    tq(tc(123456*u.s))
    <Quantity 123456. s>
    """
    return astropy.time.TimeDelta(compstr).to_value('sec')*u.s

def changeunits(qsq, unitlookup=prefunits):
    '''Change the units for the quantity or structured quantity to the system of units.'''
    if type(qsq.unit) is u.StructuredUnit:
        tounits = u.StructuredUnit(tuple([unitlookup[u.get_physical_type(un)._physical_type_list[0]] \
                                          for un in qsq.unit.values()]))
    else:
        tounits = unitlookup[u.get_physical_type(qsq.unit)._physical_type_list[0]]
    return qsq.to(tounits)

def normalizeangle(angle, wrapat=u.rev/2, exclude=[]):
    '''Add or subtract multiples of full revolutions so that angle
    falls in the semi-open range [-180 degrees, +180 degrees). The cut
    point can be changed by setting wrapat differently; for example,
    for [0, 360) degrees, set to u.rev. Parts of structured quantities
    with names listed in `exclude` are not normalized. If
    exclude==True, no values are changed. Default is to exclude
    nothing.
    '''
    from tellurion.core import nquant
    if type(angle) is u.Quantity:
        if type(angle.unit) is u.StructuredUnit and exclude != True:
            return nquant.structquant([normalizeangle(kv[1], wrapat, kv[0] in exclude)
                                  for kv in nquant.namedquant(angle).items()],
                                 angle.dtype.names)
        else:
            if u.get_physical_type(angle)=='angle' and exclude != True:
                return normalizeangle(Angle(angle), wrapat)
            else:
                return angle
    elif type(angle) is Angle:
        return u.Quantity(angle.wrap_at(wrapat))
    else:
        return angle

def isupperhalfplane(angle):
    '''Angle is in the upper half plane'''
    na = normalizeangle(angle)
    return na >= 0.0 and na <= u.rev/2

################################################################################
## Time
################################################################################

prefnumabstime = 'mjd' # Preferred numerical format for absolute time
astropy.time.Time.to_array = lambda self, format=prefnumabstime: to_array(self, format)

def to_array(tms, format=prefnumabstime):
    if tms.shape == ():
        mjds = np.array([tms.to_value(format=format)])
    else:
        mjds = tms.to_value(format)
    return mjds

def from_array(tms, format=prefnumabstime):
    return astropy.time.Time(astropy.time.Time(tms, format=format).to_value('isot'))

################################################################################
## Time series and Tables
################################################################################

def striptime(ts):
    '''Remove the time column from a TimeSeries and return a plain Table'''
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
    '''The time series `ts` starting at the specified reference time `reftime` (default is the current time) and new column showing the elapsed time from the reference time. Make a new series if `copy` is `True` (the default); otherwise, modify the original time series.'''
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
