"""
Use AstroPy definitions of time and units to define Quantity with time
(TQuantity). Make TQuanity of posvel, lonlatalt; these functions are
called by user-facing modules.
"""

import astropy.units as u
import numpy as np
import copy
from . import util

################################################################################
## Units
################################################################################

"""A revolution (full circle), useful for two-line elements (mean motion in rev/day)"""
u.rev = u.revolution = u.def_unit('revolution', 2*np.pi*u.radian)

"""User's preferred units"""
prefunits = {"time": u.second, "length": u.km, "velocity": u.km/u.second,
             "angle": u.degree, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}

# Sidereal day in seconds: u.sday.to('s'), u.sday.to(u.second)

################################################################################
## TQuantity: Quantities with time
################################################################################

# The `time` argument may be anything, but is presumed to be one of:
#   None, a Time, a number, a Quantity (with a time dimension)
# How this time is used is up to the application.

class TQuantity(u.Quantity):
    """
    An AstroPy Quantity (with possible physical dimension) and an associated time stamp.
    """
    # See https://stackoverflow.com/a/28236682/238405
    def __new__(cls, value, unit=None, time=None):
        self = super().__new__(cls, value, unit)
        self.time = time
        return self
    def __init__(self, value, unit=None, time=None):
        self.time = time
    def __repr__(self):
        if hasattr(self, 'time'):
            return f"{super().__repr__()[:-1]}, time={self.time}>"
        else:
            return f"{super().__repr__()[:-1]}, time not available>"
    def __copy__(self):
        return TQuantity(self.value, unit=self.unit, time=self.time)
    def __deepcopy__(self):
        return TQuantity(copy.deepcopy(self.value), unit=copy.deepcopy(self.unit), time=copy.deepcopy(self.time))
    def convert_units(self, units):
        newobj = self.to(u.StructuredUnit(units))
        newobj.time = self.time
        return newobj

# Add a method for equality, subtraction etc. that checks time to be within some specified difference
# Apply isclose
#   astropy.time.Time('2023-09-14T08:31:00').isclose(astropy.time.Time('2023-09-14T08:31:00.00099999'), 1*u.ms)
# tquantex1 = TQuantity(15, u.m / u.s, astropy.time..Time('2023-09-14T08:30:00'))

def quant(tquant):
    "Return the AstroPy Quantity from the TQuantity"
    return u.Quantity(tquant.value, tquant.unit)

def tquant(quant, time):
    "Make a TQuantity from the AstroPy Quantity and Time"
    return TQuantity(quant.value, quant.unit, time)

################################################################################
#### LLA: Geographic coordinates, longitude, latitude, altitude
################################################################################

# Examples
#   llatime = lonlatalt([-40.0, 0.1, 250.0], nowutc())
#   llatime["lon"] # => <Quantity -40. deg>
#   llanotime = lonlatalt([-40.0, 0.1, 250.0])
#   llanotime["lon"].to(u.radian) # => <Quantity -0.6981317 rad>
#   addtime = lonlatalt(llanotime,nowutc())
#   remtime = lonlatalt(llatime)
#   changetime = lonlatalt(llatime, nowutc())
def lonlatalt(lla, time=None, length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
    if util.listnpa(lla):
        lon = lla[0]
        lat = lla[1]
        alt = lla[2]
        llatype = [('lon', 'f8'), ('lat', 'f8'), ('alt', 'f8')]
        npa = np.array((lon, lat, alt), dtype = llatype)
        strunit = u.StructuredUnit((angle_unit, angle_unit, length_unit))
        if time==None:
            ret = u.Quantity(npa, strunit)
        else:
            ret = TQuantity(npa, strunit, time)
    elif type(lla) is u.Quantity and time is not None: # Add time to an LLA that has none
        ret = TQuantity(lla.value, lla.unit, time)
    elif type(lla) is TQuantity:  # Change or remove the time
        if time is None:
            ret = u.Quantity(lla.value, lla.unit)
        else:
            ret = TQuantity(lla.value, lla.unit, time)
    else:
        ret = None
    return(ret)

def latlonalt(lla, time=None, length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
    return (lonlatalt([lla[1], lla[0], lla[2]], time, length_unit, angle_unit))

################################################################################
#### Timeseries (or ephemeris)
################################################################################

"""
Create a list of TQuantities from a column of a time series

     import astrodynamics.example as exmp
     from astrodynamics import astro
     exmp.propdemo()
     postime = astro.tscolumn(exmp.ex1.prop.ephem, "position")
     postime[2]
        <TQuantity [5241.45723203, 3974.18748273,  821.01870343] km, time=2022-06-01 12:02:30>
"""
def tscolumn(timeseries,column):
    return([tquant(row[column],row["time"]) for row in timeseries])
