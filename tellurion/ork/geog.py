"""Computation of GeodeticPoint and ECI computation of site vectors
and angles & range observation. These definitions are not necessary as
AstroPy does the same computation (see ../geog.py), unless the
specific earth frame used for Orekit is needed to high accuracy.

"""
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.timeseries import TimeSeries

from org.orekit.bodies import GeodeticPoint
from org.orekit.frames import TopocentricFrame
import org.orekit.models.earth as oearth

from ..core import posvel
from ..core import astro
from ..core import geog
from . import force
from . import convert

def siderealtime(time = None, location = None, forceenv=force.deffe):
    '''The Greenwich or local sidereal time(s). Default `time` is the
    current time, default `location` is the prime meridian (i.e.,
    Greenwich sidereal time).
    '''
    if not time:
        time = astro.abstime(0)
    def gst(time):
        if time.isscalar:
            vs = forceenv['gmst'].value(convert._okad(time))
        else:
            vs = [gst(tm) for tm in time]
        return coord.Longitude(vs, u.radian).to(u.deg)
    if hasattr(location, 'lon'):
        return coord.Longitude(gst(time) - location.lon)
    else:
        return gst(time)

def geodpt(earthloc):
    '''Create the Orekit GeodeticPoint from an AstroPy EarthLocation'''
    geod = earthloc.geodetic
    lat_rdn = float(geod.lat.radian)
    lon_rdn = float(geod.lon.radian)
    altitude_m = float(geod.height.si.value)
    return GeodeticPoint(lat_rdn, lon_rdn, altitude_m)
