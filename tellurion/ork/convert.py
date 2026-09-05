"""Convert Cartesian state vectors to and from Orekit objects

All functions accept or return Orekit objects are for internal use
and thus begins with `_`.

"""

import collections.abc
import datetime

import astropy.time
import astropy.units as u
import numpy as np
import orekit_jpype.pyhelpers as pyhelp
import pandas as pd

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.utils import TimeStampedPVCoordinates

import tellurion.astro.time as atime
from tellurion.astro import quantity_utils as quant
from tellurion.astro import units as tunits
from tellurion.core import posvel

posvelsiu = u.StructuredUnit((u.meter, u.meter / u.second))

# ------------------------------
#    Cartesian posvel
# ------------------------------

# Convert to and from Orekit representations of position and velocity,
# or position, velocity and time.
# _pvt(): Convert from orekit objects to PV or PVT
# _tspvc(): Convert from PV/PVT to org.orekit.utils.PVCoordinates or
# TimeStampedPVCoordinates


def _pvt(object, additional=None):
    """Make the postion, velocity (), and time tuple
    (astropy.time.Time) or position and velocity from the Orekit
    object or iterable of Orekit objects that has them defined; there
    is no transformation (e.g., from Kepler elements).

    """

    def orkpv(obj):
        """Extract the position and velocity from the Orekit object
        and return a tuple of them."""

        def _quant_from_v3d(arg, unit):
            return quant.make_quantity(
                np.array([arg.getX(), arg.getY(), arg.getZ()]), unit
            )

        if hasattr(obj, "getPosition") and hasattr(obj, "getVelocity"):
            pos = _quant_from_v3d(obj.getPosition(), tunits.orkunits["length"])
            vel = _quant_from_v3d(obj.getVelocity(), tunits.orkunits["speed"])
            return (pos, vel)
        elif hasattr(obj, "getPVCoordinates"):
            return orkpv(obj.getPVCoordinates())
        elif hasattr(obj, "initialState"):
            return orkpv(obj.initialState)
        else:
            raise ValueError("Cannot find position and velocity in Orekit object")

    def orktime(obj):
        if hasattr(obj, "getDate"):
            tm = _abstime_from_okad(obj.getDate())
        elif type(additional) is astropy.time.Time:
            tm = additional
        else:
            raise ValueError(
                "Cannot find time in Orekit object nor interpret `additional` as a time"
            )
        return tm

    if isinstance(object, collections.abc.Iterable):
        pvlist = [orkpv(obj) for obj in object]
        times = atime.abstime([orktime(obj) for obj in object])
        return posvel.pvtcart(pvlist, times)
    else:
        return posvel.pvtcart(orkpv(object), orktime(object))


def _pvt_from_coordinates(obj, time, frame):
    """Get PVT by calling getPVCoordinates with time and frame.

    Parameters
    ----------
    obj : Orekit object
        Object that has a getPVCoordinates method
    time : AbsoluteDate or iterable of AbsoluteDate
        Single time or sequence of times
    frame : Frame
        Reference frame for the coordinates

    Returns
    -------
    PositionVelocityT
        Position, velocity, and time data
    """

    if isinstance(time, collections.abc.Iterable):
        pv_list = [obj.getPVCoordinates(_abstime_to_okad(tm), frame) for tm in time]
        return _pvt(pv_list, additional=time)
    else:
        return _pvt(
            obj.getPVCoordinates(_abstime_to_okad(time), frame), additional=time
        )


def _tspvc(pvt):
    """Convert tuple (posvel.pv(), astropy.time.Time) or ephemeris row
    to Orekit TimeStampedPVCoordinates or posvel.pv() to
    PVCoordinates"""
    vecp = Vector3D(pvt.position_vector.si.value.tolist())
    if pvt.has_velocity:
        vecv = Vector3D(pvt.velocity_vector.si.value.tolist())
        return TimeStampedPVCoordinates(_abstime_to_okad(pvt.time), vecp, vecv)
    else:
        return TimeStampedPVCoordinates(_abstime_to_okad(pvt.time), vecp, Vector3D.ZERO)


def _abstime_from_okad(t):
    return atime.abstime(pyhelp.absolutedate_to_datetime(t))


def _abstime_to_okad(t):
    """Convert time in any form to Orekit AbsoluteDate (okad), or
    from okad to AstroPy"""
    if type(t) is astropy.time.Time:  # AstroPy
        if t.isscalar:
            return pyhelp.datetime_to_absolutedate(t.datetime)
        else:
            return [pyhelp.datetime_to_absolutedate(s.datetime) for s in t]
    elif t is np.datetime64:  # NumPy
        return pyhelp.datetime_to_absolutedate(pd.Timestamp(t).to_pydatetime())
    elif t is datetime.datetime:  # Python
        return pyhelp.datetime_to_absolutedate(t)
    else:
        raise ValueError("Cannot convert value to Orekit AbsoluteDate")
