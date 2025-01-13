# Site vector calculation
import numpy as np
import astropy.coordinates as coord
from astropy.timeseries import TimeSeries
import cdttm
import posvel

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
