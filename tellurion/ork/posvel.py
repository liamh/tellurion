"""
Vectors and PVT (position, velocity, time) sets in Orekit
No definitions for direct use
"""

import numpy as np
import pandas as pd
import astropy.units as u
import datetime
from astropy.timeseries import TimeSeries
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
import orekit.pyhelpers as pyhelp
import org.orekit.time

from .. import astro
from .. import posvel
from . import force
from . import element

########################################
#### Convert element set to PVT     ####
########################################

def pvt(elset, time):
    '''Convert an orbital element set to Cartesian PVT'''
    return element.keplerianorbit(elset, time).pvt()

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

def okad(t):
    """ Convert time in any form to Orekit AbsoluteDate (okad), or from okad to AstroPy """
    if posvel.isdttm(t): # AstroPy
        return pyhelp.datetime_to_absolutedate(t.datetime)
    elif type(t) == np.datetime64: # NumPy
        return pyhelp.datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif type(t) == datetime.datetime:  # Python
        return pyhelp.datetime_to_absolutedate(t)
    elif type(t) is org.orekit.time.AbsoluteDate:
        return posvel.dttm(pyhelp.absolutedate_to_datetime(t))
    else:
        raise ValueException("Cannot convert value to or from Orekit AbsoluteDate")

def _pvtork(pvc, unitlookup):
    '''Convert Orekit objects to a PVT tuple, called through methods below'''
    pos = [pvc.position.x, pvc.position.y, pvc.position.z]
    vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
    pv = astro.changeunits(posvel.pv(pos, vel, astro.orkunits), unitlookup)
    if type(pvc) is TimeStampedPVCoordinates:
        return pv, okad(pvc.date)
    else:
        return pv
TimeStampedPVCoordinates.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self, unitlookup)
PVCoordinates.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self, unitlookup)
Orbit.pvt = lambda self, unitlookup=astro.prefunits: _pvtork(self.pVCoordinates, unitlookup)

def tspvc(pv, time=None):
    '''Convert PVT or ephemeris row to TimeStampedPVCoordinates or PV to PVCoordinates'''
    if posvel.isephrow(pv):
        (pv, tpvt) = posvel.pvt(pv)
        if time==None:
            return tspvc(pv, tpvt)
        elif isdttm(time):
            return tspvc(pv, time)
        elif isreltime(time):
            return tspvc(pv, tpvt+time)
    conv = pv.to(astro.posvelsiu)
    vecp = v3d(conv[posvel._eph_pos].value)
    vecv = v3d(conv[posvel._eph_vel].value)
    if time==None:
        return PVCoordinates(vecp, vecv)
    else:
        return TimeStampedPVCoordinates(okad(time), vecp, vecv)

########################################
####       Convert Cartesian        ####
########################################

TimeStampedPVCoordinates.cartesianorbit = lambda self, fe=force.deffe: \
    CartesianOrbit(self, fe['celestialframe'], fe['earthmu'])
Orbit.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
