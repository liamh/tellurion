"""Convert Cartesian state vectors to and from Orekit objects

All function accept or return Orekit objects are for internal use
and thus begins with `_`.

"""

import numpy as np
import pandas as pd
import astropy.units as u
import datetime

import orekit_jpype.pyhelpers as pyhelp
import org.orekit.time
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.propagation import Propagator, BoundedPropagator
from org.hipparchus.geometry.euclidean.threed import Vector3D

from ..core import astro
from ..core import element
from ..core import posvel
from . import force
from . import element as oelement

###############################
####  Cartesian posvel     ####
###############################

# Convert to and from Orekit representations of position and velocity,
# or position, velocity and time.
# "PV" = An object satisying ispv() (see core/posvel.py)
# "PVT" = An object satisying ispvt() (see core/posvel.py)
# _pvt(): Convert from orekit objects to PV or PVT
# _tspvc(): Convert from PV/PVT to org.orekit.utils.PVCoordinates or TimeStampedPVCoordinates

def _pvt(object, unitlookup=astro.prefunits, getpvcargs=[], additional=None):
    '''Make the postion, velocity, and time tuple (posvel.pv(),
    astropy.time.Time) or position and velocity from the Orekit object
    that has them defined; there is no transformation (e.g., from
    Kepler elements).
    '''
    if hasattr(object, 'getPosition') and hasattr(object, 'getVelocity'):
        pos = _v3d(object.getPosition(), astro.posvelsiu[0])
        vel = _v3d(object.getVelocity(), astro.posvelsiu[1])
        pv = astro.changeunits(posvel.pv(pos, vel, astro.orkunits), unitlookup)
        if hasattr(object, 'getDate'):
            return posvel.pvt(pv, _okad(object.getDate()))
        else:
            return pv
    elif hasattr(object, 'pVCoordinates'):
        return _pvt(object.pVCoordinates, unitlookup)
    elif hasattr(object, 'getPVCoordinates'):
        if len(getpvcargs)==2:
            return _pvt(object.getPVCoordinates(_okad(getpvcargs[0]), getpvcargs[1]), unitlookup)
        else:
            return _pvt(object.getPVCoordinates(), unitlookup)
    elif hasattr(object, 'initialState'):
        return _pvt(object.initialState, unitlookup)
    elif hasattr(object, 'position') and hasattr(object, 'velocity'):
        pos = _v3d(object.position, astro.posvelsiu[0])
        vel = _v3d(object.velocity, astro.posvelsiu[1])
        pv = astro.changeunits(posvel.pv(pos, vel, astro.orkunits), unitlookup)
        if type(object) is TimeStampedPVCoordinates:
            return posvel.pvt(pv, _okad(object.getDate()))
        else:
            return pv
    else:
        raise ValueError("Cannot convert value to position, value, and time (PVT)")

def _tspvc(obj, time=None):
    '''Convert tuple (posvel.pv(), astropy.time.Time) or ephemeris row to Orekit TimeStampedPVCoordinates or posvel.pv() to PVCoordinates'''
    if posvel.isephrow(obj):
        opvt = posvel.pvt(obj)
        if time==None:
            return _tspvc(opvt.pv, opvt.time)
        elif isdttm(time):
            return _tspvc(opvt.pv, time)
        elif isreltime(time):
            return _tspvc(opvt.pv, opvt.time+time)
    elif posvel.ispvt(obj):
        return _tspvc(obj.pv, obj.time)
    elif posvel.ispv(obj):
        conv = obj.to(astro.posvelsiu)
        vecp = _v3d(conv[posvel._eph_pos].value)
        vecv = _v3d(conv[posvel._eph_vel].value)
        if time==None:
            return PVCoordinates(vecp, vecv)
        else:
            return TimeStampedPVCoordinates(_okad(time), vecp, vecv)
    elif type(obj) is PVCoordinates:
        return TimeStampedPVCoordinates(_okad(time), obj)
    elif type(obj) is TimeStampedPVCoordinates:
        return obj
    else:
        raise ValueError("Cannot convert value to PVCoordinates or TimeStampedPVCoordinates")

def _okad(t):
    """ Convert time in any form to Orekit AbsoluteDate (okad), or from okad to AstroPy """
    if posvel.isdttm(t): # AstroPy
        return pyhelp.datetime_to_absolutedate(t.datetime)
    elif type(t) == np.datetime64: # NumPy
        return pyhelp.datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif type(t) == datetime.datetime:  # Python
        return pyhelp.datetime_to_absolutedate(t)
    elif type(t) is org.orekit.time.AbsoluteDate:
        return astro.abstime(pyhelp.absolutedate_to_datetime(t))
    else:
        raise ValueError("Cannot convert value to or from Orekit AbsoluteDate")

def _v3d(arg, unit=u.dimensionless_unscaled):
    '''Make a Vector3D from the argument; if the argument is a Vector3D, return the components as a u.Quantity'''
    # match/case will not work because `case list` causes an error
    argtype = type(arg)
    if argtype == list:
        return Vector3D(arg)
    elif argtype == np.ndarray:
        return Vector3D(arg.tolist())
    elif argtype == Vector3D:
        return u.Quantity([arg.getX(), arg.getY(), arg.getZ()], unit)
    else:
        raise ValueError("Cannot convert to or from Vector3D")
