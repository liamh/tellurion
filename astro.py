# Astropy definitions of time; see https://docs.astropy.org/en/stable/time/index.html
import astropy
import numpy as np
from astropy.time import Time # "apt" = astropy Time
import copy

#### Astropy units
import astropy.units as u # Define time units, e.g. nowutc() + 5*day
from astropy.units import second, minute, hour, day, Quantity, get_physical_type
# Sidereal day in seconds: u.sday.to('s'), u.sday.to(u.second)

# Define revolution as angle unit for two-line elements (mean motion
# in rev/day)
u.rev = u.revolution = u.def_unit('revolution', 2*np.pi*u.radian)

# Preferred units for user
prefunits = {"time": u.second, "length": u.km, "velocity": u.km/u.second,
             "angle": u.degree, "angular speed": u.radian/u.second,
             "dimensionless": u.dimensionless_unscaled}

def listnpa(thing):
    return type(thing) is list or type(thing) is np.ndarray

#### TQuantity: Quantities with time
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
#   Time('2023-09-14T08:31:00').isclose(Time('2023-09-14T08:31:00.00099999'), 1*u.ms)
# tquantex1 = TQuantity(15, u.m / u.s, Time('2023-09-14T08:30:00'))

def quant(tquant):
    return Quantity(tquant.value, tquant.unit)

def tquant(quant, time):
    return TQuantity(quant.value, quant.unit, time)

#### PVT: Position, velocity, and time

# Create a position-velocity as a TQuantity, call this a `pvtq`
# ex1pv = [5740.13268349499, 3314.06715, 0.0, -2.75082683526322, 4.7645718414998, 5.50165367052644]
# ex1pvtq = posvel(ex1pv, Time('2023-09-14T08:30:00'))
# ex1pvtq['p'] => <TQuantity [5740.13268349, 3314.06715   ,    0.        ] km, time not available>
# ex1pvtq.value[0] => array([5740.13268349, 3314.06715   ,    0.        ])
def posvel(pv, time=None, length_unit=prefunits["length"], velocity_unit=prefunits["velocity"]):
    if listnpa(pv):
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

# Scale position and velocity separately
# pv = posvel to scale
# pvscale = list or vector of length 2 to scale position, velocity
def scale_posvel(pv, pvscale):
    return posvel([pv.value[0]*pvscale[0], pv.value[1]*pvscale[1]], time=pv.time,
                  length_unit=pv.unit[0], velocity_unit=pv.unit[1])

# Convert posvel units
# ex1pvtsi = ex1pvt.convert_units((u.m, u.m/u.s))

#### LLA: Geographic coordinates, longitude, latitude, altitude
# Examples
#   llatime = lonlatalt([-40.0, 0.1, 250.0], nowutc())
#   llatime["lon"] # => <Quantity -40. deg>
#   llanotime = lonlatalt([-40.0, 0.1, 250.0])
#   llanotime["lon"].to(u.radian) # => <Quantity -0.6981317 rad>
#   addtime = lonlatalt(llanotime,nowutc())
#   remtime = lonlatalt(llatime)
#   changetime = lonlatalt(llatime, nowutc())
def lonlatalt(lla, time=None, length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
    if listnpa(lla):
        lon = lla[0]
        lat = lla[1]
        alt = lla[2]
        llatype = [('lon', 'f8'), ('lat', 'f8'), ('alt', 'f8')]
        npa = np.array((lon, lat, alt), dtype = llatype)
        strunit = u.StructuredUnit((angle_unit, angle_unit, length_unit))
        if time==None:
            ret = Quantity(npa, strunit)
        else:
            ret = TQuantity(npa, strunit, time)
    elif type(lla) is Quantity and time is not None: # Add time to an LLA that has none
        ret = TQuantity(lla.value, lla.unit, time)
    elif type(lla) is TQuantity:  # Change or remove the time
        if time is None:
            ret = Quantity(lla.value, lla.unit)
        else:
            ret = TQuantity(lla.value, lla.unit, time)
    else:
        ret = None
    return(ret)

def latlonalt(lla, time=None, length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
    return (lonlatalt([lla[1], lla[0], lla[2]], time, length_unit, angle_unit))
