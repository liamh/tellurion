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

def geodpt(earthloc):
    '''Create the Orekit GeodeticPoint from an AstroPy EarthLocation'''
    geod = earthloc.geodetic
    lat_rdn = float(geod.lat.radian)
    lon_rdn = float(geod.lon.radian)
    altitude_m = float(geod.height.si.value)
    return GeodeticPoint(lat_rdn, lon_rdn, altitude_m)

def _topoframe(loc, name, forceenv=force.deffe):
    return TopocentricFrame(forceenv['earth'], geodpt(loc), name)

def eciobs(loc, observation=None, name='eci obs', forceenv=force.deffe):
    '''Find the ECI position and time of the observations made from
    the location. If observation is a Time or multiple times, find the site
    vector(s). Uses Orekit.'''
    tf = _topoframe(loc, name, forceenv)
    if type(observation) is coord.sky_coordinate.SkyCoord:
        tf = TopocentricFrame(forceenv['earth'], \
                              tf.pointAtDistance(
                                  float(observation.az.radian),
                                  float(observation.alt.radian), \
                                  float(observation.distance.si.value)), \
                              name)
        dttm = observation.obstime
    elif observation==None:
        dttm = astro.abstime(0)
    else:
        dttm = observation
    if type(dttm.value) is np.ndarray:
        arr = [eciobs(loc, dt) for dt in dttm]
        datdict = {name: [ob.cartesian.xyz for ob in arr]}
        times = [ob.obstime for ob in arr]
        ts = TimeSeries(time=times, data=datdict)
        ts[name].info.format = posvel._pos_format
        return ts
    else:
        pos = convert._pvt(tf, getpvcargs=[dttm, forceenv['celestialframe']]).pv['position']
        return posvel.makept(pos, dttm)

def aer(postime, location, forceenv=force.deffe):
    '''The az-el-range for a position-time (output of `makept()`) from an observer location'''
    eci = postime.cartesian.xyz
    pt = convert._v3d(eci)
    tf = _topoframe(location, "local frame", forceenv)
    cf = forceenv['celestialframe']
    time = postime.obstime
    oktime = convert._okad(time)
    azm = tf.getAzimuth(pt, cf, oktime)
    elv = tf.getElevation(pt, cf, oktime)
    rng = tf.getRange(pt, cf, oktime)*astro.orkunits["length"] # Can't convert units .to(astro.prefunits["length"])
    return geog.azelrange(coord.Angle(azm, u.radian), coord.Angle(elv, u.radian), rng, location, time)
