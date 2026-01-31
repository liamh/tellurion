"""Custom units

This module defines additional units beyond those in astropy.units.
"""

from typing import Final
import numpy as np
from tellurion.astro import quant
from astropy.units.core import def_unit
import astropy.units as u
from astropy.coordinates import Angle

################################################################################
##### Default and SI units
################################################################################

def gravconstunits(d):
    """Add gravitational constant units to a units dictionary"""
    d["gravconst"] = d['length']**3 / d['time']**2
    return d

# Define the unit - automatically registers in astropy.units namespace
def_unit(
    ['rev', 'revolution'],
    2 * np.pi * u.radian,
    namespace=globals(),
    doc="A revolution (full circle), useful for two-line elements (mean motion in rev/day)"
)

u.add_enabled_units(rev)
u.rev = rev

#: User's preferred units
prefunits: dict = gravconstunits({"time":   u.second, "length":  u.km, "speed": u. km/u.second,
                   "angle": u.degree, "angular speed": u.radian/u.second,
                   "dimensionless": u.dimensionless_unscaled})
""" User's preferred units, used for displaying values """

siunits: Final[dict] = gravconstunits({
    "time": u.second, "length":  u.m, "speed": u.m/u. second,
    "angle": u.radian, "angular speed":   u.radian/u. second,
    "dimensionless":  u.dimensionless_unscaled
})
""" SI units """

orkunits:  Final[dict] = siunits
""" Orekit uses SI units """

################################################################################
##### Change units
################################################################################

def changeunits(qsq, unitlookup=prefunits):
    """Change the units for the quantity or structured quantity to the system of units."""
    if type(qsq.unit) is u.StructuredUnit:
        tounits = u.StructuredUnit(tuple([unitlookup[u.get_physical_type(un)._physical_type_list[0]] \
                                          for un in qsq.unit.values()]))
    else:
        tounits = unitlookup[u.get_physical_type(qsq.unit)._physical_type_list[0]]
    return qsq.to(tounits)

################################################################################
##### Angles
################################################################################

def normalizeangle(angle, wrapat=u.rev/2, exclude=[]):
    """Add or subtract multiples of full revolutions so that angle
    falls in the semi-open range [-180 degrees, +180 degrees). The cut
    point can be changed by setting wrapat differently; for example,
    for [0, 360) degrees, set to u.rev. Parts of structured quantities
    with names listed in `exclude` are not normalized. If
    exclude==True, no values are changed. Default is to exclude
    nothing.
    """
    if type(angle) is u.Quantity:
        if type(angle.unit) is u.StructuredUnit and exclude != True:
            return quant.sq_from_dict({k: normalizeangle(v, wrapat, k in exclude)
                                       for k, v in angle.to_dict().items()},
                                      angle.isscalar)
        else:
            if u.get_physical_type(angle)=='angle' and exclude != True:
                return normalizeangle(Angle(angle), wrapat)
            else:
                return angle
    elif type(angle) is Angle:
        return u.Quantity(angle.wrap_at(wrapat))
    else:
        return angle

# Export everything
__all__ = ['rev', 'revolution', 'prefunits', 'siunits', 'orkunits', 'gravconstunits', 'changeunits', 'normalizeangle']
