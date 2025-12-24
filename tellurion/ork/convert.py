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

from tellurion.core import astro
from tellurion.core import element
from tellurion.core import posvel
from tellurion.core import pvhelper
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


def _pvt(object, unitlookup=astro.prefunits, additional=None):
    '''Make the postion, velocity (), and time tuple
    (astropy.time.Time) or position and velocity from the Orekit
    object that has them defined; there is no transformation (e.g.,
    from Kepler elements).

    '''
    if hasattr(object, 'getPosition') and hasattr(object, 'getVelocity'):
        pos = _v3d(object.getPosition(), astro.posvelsiu[0])
        vel = _v3d(object.getVelocity(), astro.posvelsiu[1])
        pv = astro.changeunits(pvhelper.cartesianpv_sep(pos, vel, astro.orkunits), unitlookup)
        if hasattr(object, 'getDate'):
            tm=_okad(object.getDate())
        elif type(additional) is astropy.time.Time:
            tm=additional
        return posvel.PositionVelocityT(time=tm, cartesian=pv)
    elif hasattr(object, 'getPVCoordinates'):
        return _pvt(object.getPVCoordinates(), unitlookup)
    elif hasattr(object, 'initialState'):
        return _pvt(object.initialState, unitlookup)
    else:
        raise ValueError("Cannot convert value to position, value, and time (PVT)")

def _tspvc(obj, time=None):
    '''Convert tuple (posvel.pv(), astropy.time.Time) or ephemeris row to Orekit TimeStampedPVCoordinates or posvel.pv() to PVCoordinates'''
    if pvhelper.isephrow(obj):
        opvt = posvel.pvt(obj) # pvhelper.cartesianpv_sep(obj['position'], obj['velocity'])
        if time==None:
            return _tspvc(opvt.pv, opvt.time)
        elif isdttm(time):
            return _tspvc(opvt.pv, time)
        elif isreltime(time):
            return _tspvc(opvt.pv, opvt.time+time)
    elif posvel.ispvtcart(obj):
        return _tspvc(obj.pv, obj.time)
    elif pvhelper.ispv(obj):
        conv = obj.to(astro.posvelsiu)
        vecp = _v3d(conv[pvhelper._eph_pos].value)
        vecv = _v3d(conv[pvhelper._eph_vel].value)
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
    if pvhelper.isdttm(t): # AstroPy
        if t.isscalar:
            return pyhelp.datetime_to_absolutedate(t.datetime)
        else:
            return [pyhelp.datetime_to_absolutedate(s.datetime) for s in t]
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
    elif argtype == u.Quantity:
        return _v3d(arg.to(astro.orkunits["length"]).value)
    elif argtype == np.ndarray:
        return Vector3D(arg.tolist())
    elif argtype == Vector3D:
        return u.Quantity([arg.getX(), arg.getY(), arg.getZ()], unit)
    else:
        raise ValueError("Cannot convert to or from Vector3D")
