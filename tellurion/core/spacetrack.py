import collections
import itertools
import astropy.units as u
import spacetrack
from tellurion.core import astro
from tellurion.core import nquant
from tellurion.core import posvel


# Set stclient using registered space-track.org username and password
# Place these lines with correct username and password in ~/.ipython/profile_default/startup/50-spacetrack.py
#   import spacetrack
#   stclient = spacetrack.SpaceTrackClient(identity="myemail@example.com", password="mypw")

# Example with Sentinel 3A
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sentst = isssent['SENTINEL 3A']
# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE OSCULATING:
# sent_badpvt = tork.cartesian(tell.kepler(sentst.els, sentst.t))
# A better choice would be to use Orekit to propagate/convert; see ork/tle.py for `sent_goodpvt`.
# sentgen = tork.SGP4gen(isssent['SENTINEL 3A'], 1*u.day, {'altitude': 125.0*u.km, 'eclipse': True, 'visibility': []})
# tell.posdiff(sentgen['pvt0'].pv, sent_badpvt.pv)
# <Quantity 38.40421256 km>

MeanElementSetT = collections.namedtuple('MeanElementSetT', 'els t tle model scdata')

def mestrt(mest):
    '''Create a readable (only strings, numbers, dict) tuple from the
    `MeanElementSetT`. To save a mean element set from spacetrack and
    later recreate it without access to `space-track.org`, use this
    function to create the readable tuple, save to a Python source
    file and call `makemest()` on it.
    '''
    return(nquant.namedquant(mest.els, False, True), \
           mest.t.to_string(), \
           mest.tle, \
           mest.model, \
           mest.scdata)

def makemest(elstuples, timestr, tle, model, scdata):
    '''Create a MeanElementSetT from readable arguments'''
    return MeanElementSetT(nquant.namedquant(*elstuples), \
                           astro.abstime(timestr), tle, model, scdata)

def satdata(stdict):
    epoch = astro.abstime(stdict['EPOCH']) + float(stdict['EPOCH_MICROSECONDS'])*u.microsecond
    orbels = {'sma': float(stdict['SEMIMAJOR_AXIS'])*u.km,
              'ecc': float(stdict['ECCENTRICITY'])*u.dimensionless_unscaled,
              'inc': float(stdict['INCLINATION'])*u.deg,
              'raan': float(stdict['RA_OF_ASC_NODE'])*u.deg,
              'argper': float(stdict['ARG_OF_PERICENTER'])*u.deg,
              'ma': float(stdict['MEAN_ANOMALY'])*u.deg,
              'memo': float(stdict['MEAN_MOTION'])*u.rev/u.day,
              'memod': float(stdict['MEAN_MOTION_DOT'])*u.rev/(u.day*u.day),
              'memodd': float(stdict['MEAN_MOTION_DDOT'])*u.rev/(u.day*u.day*u.day),
              'period': float(stdict['PERIOD'])*u.min,
              'peralt': float(stdict['PERIGEE'])*u.km,
              'apoalt': float(stdict['APOGEE'])*u.km,
              'B': 12.7416*float(stdict['BSTAR'])*u.m*u.m/u.kg}
    return (orbels, epoch)

def stscdata(stdata):
    return {'name': stdata['OBJECT_NAME'],
            'type': stdata['OBJECT_TYPE'],
            'catid': int(stdata['OBJECT_NUMBER']),
            'intldes': stdata['OBJECT_ID']}

def spacetrack_latest(stclient, satnums):
    stdata = stclient.tle_latest(norad_cat_id=satnums, ordinal=1)
    sattle = (stclient.tle_latest(norad_cat_id=satnums, ordinal=1, format='tle')).splitlines()
    ret = [MeanElementSetT(*satdata(std), tle, 'SGP4', stscdata(std)) \
           for (std, tle) in zip(stdata, itertools.batched(sattle, 2))]
    if type(satnums) is int:
        return ret[0]
    else:
        return dict(zip([el.scdata['name'] for el in ret], ret))
