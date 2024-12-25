"""
Vectors and PVT (position, velocity, time) sets in Orekit and AstroPy

No definitions for direct use
"""
import numpy as np
import astropy.units as u
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from . import astro
from . import dttm

def v3d(arg):
    # match/case will not work because `case list` causes an error
    argtype = type(arg)
    if argtype == list:
        return Vector3D(arg)
    elif argtype == np.ndarray:
        return Vector3D(arg.tolist())
    elif argtype == Vector3D:
        return arg

def orkpvt(pos, vel, time=None):
    if time==None:
        return PVCoordinates(v3d(pos), v3d(vel))
    else:
        return TimeStampedPVCoordinates(dttm.to_okad(time), v3d(pos), v3d(vel))

################################################################################
#### PVT: Position, velocity, and time
################################################################################

# Create a position-velocity as a TQuantity, call this a `pvtq`
# ex1pv = [5740.13268349499, 3314.06715, 0.0, -2.75082683526322, 4.7645718414998, 5.50165367052644]
# ex1pvtq = posvel(ex1pv, astropy.time.Time('2023-09-14T08:30:00'))
# ex1pvtq['p'] => <TQuantity [5740.13268349, 3314.06715   ,    0.        ] km, time not available>
# ex1pvtq.value[0] => array([5740.13268349, 3314.06715   ,    0.        ])
def posvel(pv, time=None, length_unit=astro.prefunits["length"], velocity_unit=astro.prefunits["velocity"]):
    if astro.listnpa(pv):
        if len(pv) == 6:
            pos = pv[0:3]
            vel = pv[3:6]
        elif len(pv) == 2:
            pos = pv[0]
            vel = pv[1]
    # See https://docs.astropy.org/en/stable/units/structured_units.html#example
    pvtype = [('p', '(3,)f8'), ('v', '(3,)f8')]
    pv = np.array((pos, vel), dtype = pvtype)
    pv = astro.TQuantity(pv, u.StructuredUnit((length_unit, velocity_unit)), time)
    return(pv)

# Scale position and velocity separately
# pv = posvel to scale
# pvscale = list or vector of length 2 to scale position, velocity
def _scale_posvel(pv, pvscale):
    return posvel([pv.value[0]*pvscale[0], pv.value[1]*pvscale[1]], time=pv.time,
                  length_unit=pv.unit[0], velocity_unit=pv.unit[1])

# Convert posvel units
# ex1pvtsi = ex1pvt.convert_units((u.m, u.m/u.s))
