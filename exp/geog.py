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

def earthloc(lon, lat, elevation=None):
    '''Make an earth location; if elev=None (default), look it up.'''
    longi = coord.Angle(lon)
    latit = coord.Angle(lat)
    if elevation is None:
        elevation = geonames.elev(latit.value, longi.value)*u.m
    return coord.EarthLocation.from_geodetic(lon=longi, lat=latit, height=elevation.si)

nullisland = earthloc(0*u.deg, 0*u.deg, 0*u.m)

def eciobs(loc, observation=None, name='eci obs'):
    '''Find the ECI position and time of the observations made from
    the location. If observation is a Time or multiple times, find the site
    vector(s). Uses AstroPy.'''
    if observation==None:
        observation = posvel.nowutc()
    if type(observation) is Time:
        if type(observation.value) is np.ndarray:
            arr = [eciobs(loc, obs) for obs in observation]
            datdict = {name: [posvel.makepos(ob) for ob in arr]}
            times = [ob.obstime for ob in arr]
            ts = TimeSeries(time=times, data=datdict)
            ts[name].info.format = posvel._pos_format
            return ts
        else:
            itrs = coord.ITRS(coord.CartesianRepresentation(x=loc.x, y=loc.y, z=loc.z), \
                              obstime=observation)
            gcrs = itrs.transform_to(coord.GCRS(obstime=observation))
            return coord.SkyCoord(gcrs)
    elif type(observation) is coord.sky_coordinate.SkyCoord:
        return observation.transform_to(coord.GCRS)
    else:
        return None

def azelrange(az, el, rang, obsloc, obstime):
    return coord.SkyCoord(coord.AltAz(az=az, alt=el, distance=rang, location=obsloc, obstime=obstime))

def siderealtime(time = None, location = nullisland):
    '''
    The sidereal time of the location; `time=None` (default) gives
    the current time, default `location` gives GST
    '''
    if time==None:
        time = posvel.nowutc()
    obstm = Time(time, location = location)
    return obstm.sidereal_time('mean')
