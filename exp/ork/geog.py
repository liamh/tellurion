import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries

from org.orekit.bodies import GeodeticPoint
from org.orekit.frames import TopocentricFrame
import org.orekit.models.earth as oearth

import geog
import posvel
import cdttm
import ork.force as ofr

def geodpt(earthloc):
    '''Create the Orekit GeodeticPoint from an AstroPy EarthLocation'''
    geod = earthloc.geodetic
    lat_rdn = float(u.Quantity(geod.lat).to(u.radian).value)
    lon_rdn = float(u.Quantity(geod.lon).to(u.radian).value)
    altitude_m = float(u.Quantity(geod.height).to(u.m).value)
    return GeodeticPoint(lat_rdn, lon_rdn, altitude_m)

def sitevec(loc, dttm, name='sitevec', forceenv=ofr.deffe):
    '''Find the site vector (in the GCRS frame) of the EarthLocation
    using Orekit. If `dttm` has multiple times, then a time series is
    created of all the site vectors. If it is `None`, then the current
    time is used.
    '''
    if dttm==None:
        dttm = posvel.nowutc()
    if type(dttm.value) is np.ndarray:
        arr = [sitevec(loc, dt) for dt in dttm]
        datdict = {name: [x[0] for x in arr]}
        times = [x[1] for x in arr]
        ts = TimeSeries(time=times, data=datdict)
        ts[name].info.format = posvel._pos_format
        return ts
    topoframe = TopocentricFrame(forceenv['earth'], geodpt(loc), name)
    posv3d = topoframe.getPVCoordinates(cdttm.okad(dttm), forceenv['celestialframe']).getPosition()
    pos = posvel.makepos(posv3d.quant(u.m))
    return (pos, dttm)
