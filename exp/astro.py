"""
AstroPy Quantity
"""

import warnings
import astropy.units as u
import astropy.table
import numpy as np
import collections.abc

################################################################################
## Units
################################################################################

"""A revolution (full circle), useful for two-line elements (mean motion in rev/day)"""
u.rev = u.revolution = u.def_unit('revolution', 2*np.pi*u.radian)

"""User's preferred units"""
prefunits = {"time": u.second, "length": u.km, "speed": u.km/u.second,
             "angle": u.degree, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}
prefunits["posvel"] = (prefunits["length"], prefunits["speed"])
posvelsiu = u.StructuredUnit((u.meter, u.meter/u.second))

def timesec(t):
    '''Convert a Quantity to seconds as a Python float'''
    if type(t) is u.Quantity and u.get_physical_type(t) == 'time':
        pt = t.si.value.tolist() # convert to seconds and get the value_unit
    elif isinstance(t, collections.abc.Iterable):
        return [timesec(i) for i in t]
    else:
        pt = float(t) # assume seconds
    return(pt)


def sq(values, names, sizes, phystype, unitlookup=prefunits):
    '''Make a structured quantity from numbers'''
    # sq([[1,2,3],[4,5,6]], ('pos','vel'), (3, 3), ('length', 'speed'))
    # sq([12345.0, 45.0], ('sma','inc'), (1,1), ('length', 'angle'))
    def scvec(size):
        if size==1:
            return f"f8"
        else:
            return f"({size},)f8"
    dtype = [(n, scvec(s)) for (n, s) in zip(names, sizes)]
    npa = np.array(tuple(values), dtype = dtype)
    units = tuple([unitlookup[pt] for pt in phystype])
    return u.Quantity(npa, u.StructuredUnit(units))
# Now have it take u.Q as input
# tl = [oel.elementval(okep1h, 'sma'), oel.elementval(okep1h, 'ecc'), oel.elementval(okep1h, 'inc')]

################################################################################
## Time series and Tables
################################################################################

def striptime(ts):
    '''Remove the time column from a TimeSeries and return a plain Table'''
    return astropy.table.Table([ts[k] for k in ts.keys()[1:None]])

def hcat(ts1, ts2):
    '''Concatenate timeseries by adding columns; warn if times are not all equal'''
    if not(all(ts1['time'].__eq__(ts2['time']))):
        warnings.warn("Times are not all equal; using times from first set")
    return astropy.table.hstack([ts1, striptime(ts2)])
