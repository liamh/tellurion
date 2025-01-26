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
import ork.force

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

def _pvtork(pvc, unitlookup):
    '''Convert Orekit objects to a PVT tuple, called through methods below'''
    pos = [pvc.position.x, pvc.position.y, pvc.position.z]
    vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
    pv = astro.changeunits(posvel.pv(pos, vel, astro.orkunits), unitlookup)
    if type(pvc) is TimeStampedPVCoordinates:
        return pv, cdttm.dttm(pvc.date)
    else:
        return pv
TimeStampedPVCoordinates.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self, unitlookup)
PVCoordinates.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self, unitlookup)
Orbit.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self.pVCoordinates, unitlookup)

def pvt(pv, time=None):
    '''Convert PVT or ephemeris row to TimeStampedPVCoordinates or PV to PVCoordinates'''
    if posvel.isephrow(pv):
        (pv, tpvt) = posvel.pvt(pv)
        if time==None:
            return pvt(pv, tpvt)
        elif isdttm(time):
            return pvt(pv, time)
        elif isreltime(time):
            return pvt(pv, tpvt+time)
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

TimeStampedPVCoordinates.cartesianorbit = lambda self, fe=ork.force.deffe: \
    CartesianOrbit(self, fe['celestialframe'], fe['earthmu'])
Orbit.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
