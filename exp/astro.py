"""
AstroPy Quantity
"""

import astropy.units as u
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


# Convert a Quantity to seconds as a Python float
def timesec(t):
    if type(t) is u.Quantity and u.get_physical_type(t) == 'time':
        pt = t.si.value.tolist() # convert to seconds and get the value_unit
    elif isinstance(t, collections.abc.Iterable):
        return [timesec(i) for i in t]
    else:
        pt = float(t) # assume seconds
    return(pt)
