"""Unused conversions and observation formats. Where AstroPy
functionality duplicates Orekit's, we usually prefer Orekit for
consistency.

"""

from astropy.timeseries import TimeSeries
from astropy.time import Time

def eciobs(loc, observation=None, name='eci obs'):
    '''Find the ECI position and time of the observations made from
    the location. If observation is a Time or multiple times, find the site
    vector(s). Uses AstroPy.'''
    if observation==None:
        observation = astro.abstime(0)
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
    '''Create an azimuth, elevation, and range observation'''
    return coord.SkyCoord(coord.AltAz(az=az, alt=el, distance=rang, location=obsloc, obstime=obstime))

def radecrange(ra, dec, rang, obsloc, obstime, frame='gcrs'):
    '''Create an right ascension, declination, and range observation'''
    # This doesn't seem to offset the origin to the observer location
    return coord.SkyCoord(ra=ra, dec=dec, distance=rang, obstime=obstime, frame=frame, \
                          obsgeoloc=posvel.makepos(obsloc))

def siderealtime(time = None, location = nullisland):
    '''
    The sidereal time of the location; `time=None` (default) gives
    the current time, default `location` gives GST
    '''
    if time==None:
        time = astro.abstime(0)
    obstm = Time(time, location = location)
    return obstm.sidereal_time('mean')
