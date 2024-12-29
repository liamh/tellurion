"""
Vectors and PVT (position, velocity, time) sets in Orekit
No definitions for direct use
"""

import numpy as np
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from . import dttm

v3dnan = Vector3D(np.nan,np.nan,np.nan) # used in Orekit when there is no velocity specified

def v3d(arg):
    '''Make a Vector3D from the argument'''
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
    '''Convert TimeStampedPVCoordinates or PVCoordinates to a PVT tuple'''
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
    '''Convert PVT to TimeStampedPVCoordinates or PV to PVCoordinates'''
    if time==None:
        return PVCoordinates(v3d(pos), v3d(vel))
    else:
        return TimeStampedPVCoordinates(dttm.to_okad(time), v3d(pos), v3d(vel))
