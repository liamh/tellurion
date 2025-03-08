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

from .. import posvel
from . import force
from . import posvel as oposvel

def geodpt(earthloc):
    '''Create the Orekit GeodeticPoint from an AstroPy EarthLocation'''
    geod = earthloc.geodetic
    lat_rdn = float(geod.lat.radian)
    lon_rdn = float(geod.lon.radian)
    altitude_m = float(geod.height.si.value)
    return GeodeticPoint(lat_rdn, lon_rdn, altitude_m)

def eciobs(loc, observation=None, name='eci obs', forceenv=force.deffe):
    '''Find the ECI position and time of the observations made from
    the location. If observation is a Time or multiple times, find the site
    vector(s). Uses Orekit.'''
    tf = TopocentricFrame(forceenv['earth'], geodpt(loc), name)
    if type(observation) is coord.sky_coordinate.SkyCoord:
        tf = TopocentricFrame(forceenv['earth'], \
                              tf.pointAtDistance(
                                  float(observation.az.radian),
                                  float(observation.alt.radian), \
                                  float(observation.distance.si.value)), \
                              name)
        dttm = observation.obstime
    elif observation==None:
        dttm = posvel.nowutc()
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
        posv3d = tf.getPVCoordinates(oposvel.okad(dttm), forceenv['celestialframe']).getPosition()
        pos = posvel.makepos(posv3d.quant(u.m))
        return posvel.makept(pos, dttm)
