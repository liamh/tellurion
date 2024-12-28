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
from . import util

# posvel preferred and SI units (for Orekit)
astro.prefunits["posvel"] = (astro.prefunits["length"], astro.prefunits["velocity"])
posvelsiu = (u.meter, u.meter/u.second)

v3dnan = Vector3D(np.nan,np.nan,np.nan) # used in Orekit when there is no velocity specified

def v3d(arg):
    # match/case will not work because `case list` causes an error
    argtype = type(arg)
    if argtype == list:
        return Vector3D(arg)
    elif argtype == np.ndarray:
        return Vector3D(arg.tolist())
    elif argtype == Vector3D:
        return arg
    elif argtype == NoneType:
        return v3dnan

def pvtork(pvc):
    '''Convert TimeStampedPVCoordinates to PVT'''
    pos = [pvc.position.x, pvc.position.y, pvc.position.z]
    if pvc.velocity.x==v3dnan.x:
        vel=None
    else:
        vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
    if type(pvc) is TimeStampedPVCoordinates:
        time=pvc.date.apt()
    else:
        time=None
    return (pos, vel, time)

def orkpvt(pos, vel, time=None):
    '''Convert PVT to TimeStampedPVCoordinates'''
    if time==None:
        return PVCoordinates(v3d(pos), v3d(vel))
    else:
        return TimeStampedPVCoordinates(dttm.to_okad(time), v3d(pos), v3d(vel))

################################################################################
#### PVT: Position, velocity, and time
################################################################################

def pvsplit(pv):
    '''Create a tuple of position and (optionally) velocity 3-vectors (or lists)'''
    match len(pv):
        case 6:            # Position and velocity
            pos = pv[0:3]
            vel = pv[3:6]
        case 2:            # Position and velocity
            pos = pv[0]
            vel = pv[1]
        case 3:            # Position only
            pos = pv[0:3]
            vel = None
        case 1:            # Position only
            pos = pv[0]
            vel = None
    return (pos,vel)

def atqpvt(pos, vel, time=None, units=astro.prefunits["posvel"]):
    '''Create a position-velocity as an astro.TQuantity'''
    if np.isnan(vel[0]):
        pvtype = [('p', np.float64)]
        npa = np.array(pos, dtype = pvtype)
        pv = astro.TQuantity(npa, u.StructuredUnit(units[0]), time)
    else:
        # See https://docs.astropy.org/en/stable/units/structured_units.html#example
        pvtype = [('p', '(3,)f8'), ('v', '(3,)f8')]
        npa = np.array((pos, vel), dtype = pvtype)
        pv = astro.TQuantity(npa, u.StructuredUnit(units), time)
    return(pv)

# Scale position and velocity separately
# pv = posvel to scale
# pvscale = list or vector of length 2 to scale position, velocity
def _scale_posvel(pv, pvscale):
    return posvel([pv.value[0]*pvscale[0], pv.value[1]*pvscale[1]], time=pv.time,
                  length_unit=pv.unit[0], velocity_unit=pv.unit[1])

# Convert posvel units
# ex1pvtsi = ex1pvt.convert_units((u.m, u.m/u.s))
