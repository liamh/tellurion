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

########################################
####    Vector3D                    ####
########################################

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

Vector3D.quant = lambda self, unit: u.Quantity([self.x, self.y, self.z], unit)

########################################
####  PVT tuple to and from Orekit  ####
########################################

def pvtork(pvc):
    '''Convert Orekit objects to a PVT tuple, called through methods below'''
    pos = [pvc.position.x, pvc.position.y, pvc.position.z]
    vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
    if type(pvc) is TimeStampedPVCoordinates:
        return posvel.pv(posvel.pv(pos, vel, astro.posvelsiu), None), cdttm.dttm(pvc.date)
    else:
        return posvel.pv(posvel.pv(pos, vel, astro.posvelsiu), None)

def orkpvt(pv, time=None):
    '''Convert PVT to TimeStampedPVCoordinates or PV to PVCoordinates'''
    conv = pv.to(astro.posvelsiu)
    vecp = v3d(conv[posvel._eph_pos].value)
    vecv = v3d(conv[posvel._eph_vel].value)
    if time==None:
        return PVCoordinates(vecp, vecv)
    else:
        return TimeStampedPVCoordinates(cdttm.okad(time), vecp, vecv)

########################################
####       Convert Cartesian        ####
########################################

TimeStampedPVCoordinates.pvt = lambda self: pvtork(self)
TimeStampedPVCoordinates.cartesianorbit = lambda self, fe=ofr.deffe: \
    CartesianOrbit(self, fe['celestialframe'], fe['earthmu'])
Orbit.pvt = lambda self: pvtork(self.pVCoordinates)
Orbit.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
