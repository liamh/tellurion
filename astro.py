# Astropy definitions of time; see https://docs.astropy.org/en/stable/time/index.html
import astropy
import numpy as np
from astropy.time import Time # "apt" = astropy Time
import astropy.units as u # Define time units, e.g. nowutc() + 5*day
from astropy.units import second, minute, hour, day, Quantity, get_physical_type
# Sidereal day in seconds: u.sday.to('s'), u.sday.to(u.second)

# https://docs.astropy.org/en/stable/timeseries/times.html
# from astropy.timeseries import TimeSeries

######## TQuantity: Quantities with time
# The `time` argument may be anything, but is presumed to be one of:
#   None, a Time, a number, a Quantity (with a time dimension)
# How this time is used is up to the application.

class TQuantity(Quantity):
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
    def convert_units(self, units):
        newobj = self.to(u.StructuredUnit(units))
        newobj.time = self.time
        return newobj

# Add a method for equality, subtraction etc. that checks time to be within some specified difference
# Apply isclose
#   Time('2023-09-14T08:31:00').isclose(Time('2023-09-14T08:31:00.00099999'), 1*u.ms)
# tquantex1 = TQuantity(15, u.m / u.s, Time('2023-09-14T08:30:00'))

######## PVT: Position, velocity, and time

default_length_unit = u.km
default_velocity_unit = u.km/u.second

# Create a position-velocity as a TQuantity, call this a `pvtq`
# ex1pv = [5740.13268349499, 3314.06715, 0.0, -2.75082683526322, 4.7645718414998, 5.50165367052644]
# ex1pvtq = posvel(ex1pv, Time('2023-09-14T08:30:00'))
# ex1pvtq['p'] => <TQuantity [5740.13268349, 3314.06715   ,    0.        ] km, time not available>
# ex1pvtq.value[0] => array([5740.13268349, 3314.06715   ,    0.        ])
def posvel(pv, time=None, length_unit=default_length_unit, velocity_unit=default_velocity_unit):
    if type(pv) is list:
        if len(pv) == 6:
            pos = pv[0:3]
            vel = pv[3:6]
        elif len(pv) == 2:
            pos = pv[0]
            vel = pv[1]
    # See https://docs.astropy.org/en/stable/units/structured_units.html#example
    pvtype = [('p', '(3,)f8'), ('v', '(3,)f8')]
    pv = np.array((pos, vel), dtype = pvtype)
    pv = TQuantity(pv, u.StructuredUnit((length_unit, velocity_unit)), time)
    return(pv)

# Convert posvel units
# ex1pvtsi = ex1pvt.convert_units((u.m, u.m/u.s))
