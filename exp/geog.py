# Site vector calculation
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.timeseries import TimeSeries
from astropy.time import Time
import cdttm
import posvel
import geonames

# List of sites: astropy.coordinates.EarthLocation.get_site_names()
def sitevec(loc, dttm=None, name='sitevec'):
    '''Find the site vector (in the GCRS frame) of the EarthLocation
    using AstroPy. If `dttm` has multiple times, then a time series is
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
    itrs = coord.ITRS(coord.CartesianRepresentation(x=loc.x, y=loc.y, z=loc.z), obstime=dttm)
    gcrs = itrs.transform_to(coord.GCRS(obstime=dttm)).cartesian
    return (posvel.makepos(gcrs), dttm)

def earthloc(lon, lat, elevation=None):
    '''Make an earth location; if elev=None (default), look it up.'''
    longi = coord.Angle(lon)
    latit = coord.Angle(lat)
    if elevation is None:
        elevation = geonames.elev(latit.value, longi.value)*u.m
    return coord.EarthLocation.from_geodetic(lon=longi, lat=latit, height=elevation.si)

nullisland = earthloc(0*u.deg, 0*u.deg, 0*u.m)

def siderealtime(time = None, location = nullisland):
    '''
    The sidereal time of the location; `time=None` (default) gives
    the current time, default `location` gives GST
    '''
    if time==None:
        time = posvel.nowutc()
    obstm = Time(time, location = location)
    return obstm.sidereal_time('mean')
