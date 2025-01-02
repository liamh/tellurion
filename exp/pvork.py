"""
Vectors and PVT (position, velocity, time) sets in Orekit
No definitions for direct use
"""

import numpy as np
from astropy.timeseries import TimeSeries
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType

#from .
import dttm
#from .
import astro
#from .
import posvel
#from .
import force

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

def pvtork(pvc):
    '''Convert TimeStampedPVCoordinates or PVCoordinates to a PVT tuple'''
    if type(pvc) in [CartesianOrbit, Orbit]:
        return pvtork(pvc.pVCoordinates)
    else:
        pos = [pvc.position.x, pvc.position.y, pvc.position.z]
        vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
        if type(pvc) is TimeStampedPVCoordinates:
            return posvel.makepv(posvel.makepv(pos, vel, astro.posvelsiu), None), pvc.date.apt()
        else:
            return posvel.makepv(posvel.makepv(pos, vel, astro.posvelsiu), None)

def orkpvt(pv, time=None):
    '''Convert PVT to TimeStampedPVCoordinates or PV to PVCoordinates'''
    conv = pv.to(astro.posvelsiu)
    vecp = v3d(conv['p'].value)
    vecv = v3d(conv['v'].value)
    if time==None:
        return PVCoordinates(vecp, vecv)
    else:
        return TimeStampedPVCoordinates(dttm.to_okad(time), vecp, vecv)

deffe = force.setgravity(0,0)

TimeStampedPVCoordinates.cartesian = lambda self, gravity=deffe: CartesianOrbit(self, gravity['celestialframe'], gravity['earthmu'])
