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

import tellurion.astro.time as atime
from tellurion.astro import units
from tellurion.astro import quantity_utils as quant
from tellurion.core import posvel
from tellurion.core import geog
from tellurion.core import obs
from tellurion.ork import force
from tellurion.ork import convert

def siderealtime(time = None, location = None, forceenv=force.deffe):
    '''The Greenwich or local sidereal time(s). Default `time` is the
    current time, default `location` is the prime meridian (i.e.,
    Greenwich sidereal time).
    '''
    if not time:
        time = atime.abstime(0)
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

def _topoframe(loc, name="", forceenv=force.deffe):
    return TopocentricFrame(forceenv['earth'], geodpt(loc), name)

def sitevec(loc, times, name="sitevec", unitlookup=units.prefunits, forceenv=force.deffe):
    '''The site vector for an earth location at the specified times.'''
    if times.isscalar:
        tf = _topoframe(loc, name, forceenv)
        tspvc = tf.getPVCoordinates(convert._okad(times), forceenv['celestialframe'])
        ret = convert._pvt(tspvc, unitlookup)
    else:
        pvts = [sitevec(loc, t, name, unitlookup, forceenv) for t in times]
        ret = posvel.PositionVelocityT(time=times, cartesian=\
                                       quant.vstack(tuple([pvt.cartesian for pvt in pvts])))
    if name:
        ret.name = name
    return ret

def eciaer(observation, name="", unitlookup=units.prefunits, forceenv=force.deffe):
    '''Find the ECI position and time of the az-el-range observations
    made from the location. If observation is a Time or multiple
    times, find the site vector(s). Uses Orekit.
    '''
    tf = _topoframe(observation.loc, name, forceenv)
    aerork = quant.change_units(observation.obs, units.orkunits)
    tf = TopocentricFrame(forceenv['earth'], \
                          tf.pointAtDistance(
                              float(aerork['azim'].value), \
                              float(aerork['elev'].value), \
                              float(aerork['range'].value)), \
                          name)
    tspvc = tf.getPVCoordinates(convert._okad(observation.time), forceenv['celestialframe'])
    pvt = convert._pvt(tspvc, unitlookup)
    return pvt.position

def aereci(postime, location, forceenv=force.deffe):
    '''The EarthObservation (az-el-range) for a PositionT from an observer location'''
    tf = _topoframe(location, "local frame", forceenv)
    cf = forceenv['celestialframe']
    def aer(postime):
        if postime.isscalar:
            pt = convert._v3d(postime.cartesian)
            oktime = convert._okad(postime.time)
            azm = tf.getAzimuth(pt, cf, oktime)
            elv = tf.getElevation(pt, cf, oktime)
            rng = tf.getRange(pt, cf, oktime)
            return [azm, elv, rng]
        else:
            return [aer(pt) for pt in postime]
    data = np.array(aer(postime))
    if data.ndim == 1:
        azm = coord.Angle(data[0], u.radian)
        elv = coord.Angle(data[1], u.radian)
        rng = data[2]*u.m
    else:
        azm = coord.Angle(data[:,0], u.radian)
        elv = coord.Angle(data[:,1], u.radian)
        rng = data[:,2]*u.m
    return obs.azelrange(azm, elv, rng, location, postime.time)
