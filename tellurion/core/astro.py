"""
AstroPy definitions
"""

import warnings
import astropy.units as u
from astropy.coordinates import Angle
import astropy.table
import astropy.time
import numpy as np
import collections.abc

################################################################################
## Units
################################################################################

#: A revolution (full circle), useful for two-line elements (mean motion in rev/day)
u.rev = u.revolution = u.def_unit('revolution', 2*np.pi*u.radian)

#: User's preferred units
prefunits = {"time": u.second, "length": u.km, "speed": u.km/u.second,
             "angle": u.degree, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}
prefunits["posvel"] = (prefunits["length"], prefunits["speed"])
#: Unit for posvel (m, m/s)
posvelsiu = u.StructuredUnit((u.meter, u.meter/u.second))
#: Units used by Orekit
orkunits = {"time": u.second, "length": u.m, "speed": u.m/u.second,
             "angle": u.radian, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}

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
    """Convert a u.Quantity with physical dimension time to a string of duration (time interval) components

    tcomponents(12*u.day + 17.3*u.hour + 5*u.min + 33.1*u.s)
    '12d 17hr 23min 33.1s'
    """
    return astropy.time.TimeDelta(t).quantity_str

def tq(compstr):
    """Convert a string of duration (time interval) components to a u.Quantity with physical dimension time

    tquantity('12d 17hr 23min 33.1s')
    <Quantity 1099413.1 s>
    """
    return astropy.time.TimeDelta(compstr).to_value('sec')*u.s

# sq2split = splitsq(sq2)
# sq2sq = makesq(*sq2split)

# kep1vals = {"argper":66.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25, "sma":8000.0}
# kep1pt = {"sma":'length', "ecc":'dimensionless', "inc":'angle', "argper":'angle', "raan":'angle', "ma":'angle'}
# kep1 = makesq(kep1vals, phystype=kep1pt)

def makesq(values, names=None, units=None, phystype=None, unitlookup=prefunits):
    """Make a structured quantity from numbers or quantities

    Parameters
    ----------
    values : dict of numbers, dict of `u.Quantity`, list or array of numbers, list or array of `u.Quantity`
      If a dict is given, the keys are the names of the structure elements, and `names` is ignored
      If numbers are given for values, `units` or `phystype` must be specified
      If `u.Quanitity` are given for values, they are converted to `units` or `phystype` and `unitlookup`
      Each value can be scalars or vectors
    names : the names of the fields of the structure;
      size must match number of values (ignored if `values` is a dict), optional
    units : the units the structure; size must match number of values, optional
    phystype : the physical type (dimension) of the each argument, units are looked
      up in `unitlookup`, may specified instead of `units`, optional
    unitlookup : a dict of physical type strings and preferred units, defaults to `prefunits`, optional

    Returns
    -------
    A structured `u.Quantity`

    Examples
    --------
    >>> makesq([[1,2,3],[4,5,6]], ('pos','vel'), phystype=('length', 'speed'))
    <Quantity ([1., 2., 3.], [4., 5., 6.]) (km, km / s)>

    >>> makesq([12345.0, 45.0], ('sma','inc'), phystype=('length', 'angle'))
    <Quantity (12345., 45.) (km, deg)>

    >>> makesq([12345.0, 45.0], ('sma','inc'), ('km', 'deg'))
    <Quantity (12345., 45.) (km, deg)>

    >>> makesq([12345.0*u.km, 45.0*u.deg], ('sma','inc'), units = ('meter', 'radian'))
    <Quantity (12345000., 0.78539816) (m, rad)>

    >>> makesq([12345000*u.m, 0.125*u.rev], ('sma','inc'), phystype = ('length', 'angle'))
    <Quantity (12345., 45.) (km, deg)>

    >>> makesq({'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, units = ('meter', 'radian'))
    <Quantity (12345000., 0.78539816) (m, rad)>
    """
    def conv(val, unit):
        if type(val) is u.Quantity:
            return val.to(unit).value
        else:
            return val
    def scvec(size):
        if size==1:
            return f"f8"
        else:
            return f"({size},)f8"
    def pqlen(item):
        if type(item) is u.Quantity:
            if type(item.value) is np.ndarray:
                return len(item)
            else:
                return 1
        else:
            if type(item) is list or type(item) is tuple or type(item) is np.ndarray:
                return len(item)
            else:
                return 1
    if type(values) is dict:
        names = tuple(values.keys())
        vals = tuple(values.values())
        if type(vals[0]) is u.Quantity:  # convert units unless units=unitlookup=None
            if units==None:
                units = tuple([v.unit for v in vals])
                vals = tuple([v.value for v in vals])
            else:
                vals = [conv(v, un) for (v, un) in zip(vals, units)]
        else:
            # If values is a dict, units/phystype must also be a dict
            if units==None:
                units = tuple([unitlookup[phystype[nm]] for nm in names])
            else:
                units = tuple(units[nm] for nm in names) # units must also be a dict
    else:
        if units==None:
            if phystype==None: # values is a quantity
                units = tuple([v.unit for v in values])
                vals = tuple([v.value for v in values])
            else:
                units = tuple([unitlookup[pt] for pt in phystype])
                vals = tuple([conv(v, u) for (v, u) in zip(values, units)])
        else:
            if type(values[0]) is u.Quantity:
                vals = [conv(v, un) for (v, un) in zip(values, units)]
            else:
                vals = values
    sizes = [pqlen(v) for v in vals]
    dtype = [(n, scvec(s)) for (n, s) in zip(names, sizes)]
    npa = np.array(tuple(vals), dtype = dtype)
    return u.Quantity(npa, u.StructuredUnit(units))

def splitsq(stqu, quant=True):
    '''Make a dict of names and quantities, or values, names, and units from the structured quantity'''
    if quant:
        vnu = splitsq(stqu, False)
        qs = [v*u for (v, u) in zip(vnu[0], vnu[2])]
        return {un:val for (val, un) in zip(qs, vnu[1])}
    else:
        return (stqu.value.tolist(), stqu.dtype.names, stqu.unit.values())

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
    falls in the semi-open range [-180 degrees, +180 degrees). Parts
    of structured quantities with names listed in `exclude` are not
    normalized.
    '''
    if type(angle) is u.Quantity:
        if type(angle.unit) is u.StructuredUnit:
            return [normalizeangle(kv[1], wrapat, kv[0] in exclude) for kv in splitsq(angle).items()]
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
## Time series and Tables
################################################################################

def striptime(ts):
    '''Remove the time column from a TimeSeries and return a plain Table'''
    return astropy.table.Table([ts[k] for k in ts.keys()[1:None]])

def hcat(ts1, ts2):
    """Concatenate timeseries by adding columns

    Warns
    -----
    If times are not all equal in both timeseries
    """
    if not(all(ts1['time'].__eq__(ts2['time']))):
        warnings.warn("Times are not all equal; using times from first set")
    return astropy.table.hstack([ts1, striptime(ts2)])

def ts(data, times):
    datdict = {_eph_pos: [d[_eph_pos] for d in data], \
               _eph_vel: [d[_eph_vel] for d in data]}
    ts = TimeSeries(time=times, data=datdict)
    ts[_eph_pos].info.format = _pos_format
    ts[_eph_vel].info.format = _vel_format
    return ts
