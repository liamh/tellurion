"""Convert to and from Orekit objects"""

import numpy as np
import pandas as pd
import astropy.units as u
import datetime

import orekit_jpype.pyhelpers as pyhelp
import org.orekit.time
from org.orekit.orbits import Orbit, CartesianOrbit, KeplerianOrbit, OrbitType
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState
from org.hipparchus.geometry.euclidean.threed import Vector3D

from ..core import astro
from ..core import element
from ..core import posvel
from . import force
from . import element as oelement

###############################
####  Cartesian posvel     ####
###############################

def pvt(object, unitlookup=astro.prefunits):
    '''Make the postion, velocity, and time tuple (posvel.pv(),
    astropy.time.Time) or position and velocity from the Orekit object
    that has them defined; there is no transformation (e.g., from
    Kepler elements).
    '''
    if hasattr(object, 'pVCoordinates'):
        return pvt(object.pVCoordinates, unitlookup)
    elif hasattr(object, 'initialState'):
        return pvt(object.initialState)
    elif hasattr(object, 'position') and hasattr(object, 'velocity'):
        pos = [pvc.position.x, pvc.position.y, pvc.position.z]
        vel = [pvc.velocity.x, pvc.velocity.y, pvc.velocity.z]
        pv = astro.changeunits(posvel.pv(pos, vel, astro.orkunits), unitlookup)
        if type(pvc) is TimeStampedPVCoordinates:
            return pv, okad(pvc.date)
        else:
            return pv
    else:
        raise ValueError("Cannot convert value to position, value, and time (PVT)")

def tspvc(pv, time=None):
    '''Convert tuple (posvel.pv(), astropy.time.Time) or ephemeris row to Orekit TimeStampedPVCoordinates or posvel.pv() to PVCoordinates'''
    if posvel.isephrow(pv):
        (pv, tpvt) = posvel.pvt(pv)
        if time==None:
            return tspvc(pv, tpvt)
        elif isdttm(time):
            return tspvc(pv, time)
        elif isreltime(time):
            return tspvc(pv, tpvt+time)
    elif posvel.ispvt(pv):
        return tspvc(pv[0], pv[1])
    elif posvel.ispv(pv):
        conv = pv.to(astro.posvelsiu)
        vecp = v3d(conv[posvel._eph_pos].value)
        vecv = v3d(conv[posvel._eph_vel].value)
        if time==None:
            return PVCoordinates(vecp, vecv)
        else:
            return TimeStampedPVCoordinates(okad(time), vecp, vecv)
    else:
        raise ValueError("Cannot convert value to PVCoordinates or TimeStampedPVCoordinates")

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
        raise ValueError("Cannot convert value to or from Orekit AbsoluteDate")

def v3d(arg, unit=u.dimensionless_unscaled):
    '''Make a Vector3D from the argument; if the argument is a Vector3D, return the components as a u.Quantity'''
    # match/case will not work because `case list` causes an error
    argtype = type(arg)
    if argtype == list:
        return Vector3D(arg)
    elif argtype == np.ndarray:
        return Vector3D(arg.tolist())
    elif argtype == Vector3D:
        return u.Quantity([self.x, self.y, self.z], unit)
    else:
        raise ValueError("Cannot convert to or from Vector3D")

###############################
####  Orbital elements     ####
###############################

def kepler(object, forceenv=force.deffe): # Add prefunits
    '''The Kepler element set from the Cartesian PVT'''
    cartorb = cartesianorbit(object, forceenv)
    keporb = OrbitType.KEPLERIAN.convertType(cartorb)
    kepels = astro.makesq(oelement.elementval(keporb, element.kepeltma_names), element.kepeltma_names)
    return (kepels, object[1])

def cartesianorbit(object, forceenv=force.deffe):
    '''Find the Orekit CartesianOrbit from the object'''
    if posvel.ispvt(object) or posvel.isephrow(object):
        return cartesianorbit(tspvc(object), forceenv)
    elif type(object) is Orbit:
        return OrbitType.CARTESIAN.convertType(object)
    elif type(object) is SpacecraftState:
        return cartesianorbit(object.orbit)
    elif type(object) is BoundedPropagator:
        return cartesianorbit(object.initialState)
    elif type(object) is TimeStampedPVCoordinates:
        return CartesianOrbit(object, forceenv['celestialframe'], forceenv['earthmu'])
    else:
        raise ValueError("Cannot convert value to CartesianOrbit")
