"""
Vectors and PVT (position, velocity, time) sets in Orekit
No definitions for direct use
"""

import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType

#from .
import cdttm
#from .
import astro
#from .
import posvel
#from .
import ork.force as ofr


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
    '''Convert Orekit objects to a PVT tuple'''
    pos = [pvc.position.x, pvc.position.y, pvc.position.z]
    vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
    if type(pvc) is TimeStampedPVCoordinates:
        return posvel.pv(posvel.pv(pos, vel, astro.posvelsiu), None), cdttm.dttm(pvc.date)
    else:
        return posvel.pv(posvel.pv(pos, vel, astro.posvelsiu), None)

TimeStampedPVCoordinates.pvt = lambda self: pvtork(self)
Orbit.pvt = lambda self: pvtork(self.pVCoordinates)
CartesianOrbit.pvt = lambda self: pvtork(self.pVCoordinates)
TimeStampedPVCoordinates.cartesianorbit = lambda self, gravity=ofr.deffe: \
    CartesianOrbit(self, gravity['celestialframe'], gravity['earthmu'])
Vector3D.quant = lambda self, unit: u.Quantity([self.x, self.y, self.z], unit)

def orkpvt(pv, time=None):
    '''Convert PVT to TimeStampedPVCoordinates or PV to PVCoordinates'''
    conv = pv.to(astro.posvelsiu)
    vecp = v3d(conv[posvel._eph_pos].value)
    vecv = v3d(conv[posvel._eph_vel].value)
    if time==None:
        return PVCoordinates(vecp, vecv)
    else:
        return TimeStampedPVCoordinates(cdttm.okad(time), vecp, vecv)
