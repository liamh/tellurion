import os
import astropy.units as u
import spacetrack
from . import astro
from . import posvel

# Set stclient using registered space-track.org username and password
# Place these lines with correct username and password in ~/.ipython/profile_default/startup/50-spacetrack.py
#   import spacetrack
#   stclient = spacetrack.SpaceTrackClient(identity="myemail@example.com", password="mypw")

# Example with Sentinel 3A
# sentst = tell.spacetrack_latest(stclient, 41335)
# sentkep = tell.makesq(sentst[0][0])  # Keplerian elements
# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE OSCULATING:
# sent_badpvt = tork.keplerianorbit(sentkep, sentst[0][1]).pvt()
# A better choice would be to use Orekit to propagate/convert; see ork/tle.py for `sent_goodpvt`.
# tell.posdiff(sent_goodpvt[0], sent_badpvt[0])
# <Quantity 38.40421256 km>

def satdata(stdict):
    epoch = posvel.dttm(stdict['EPOCH']) + float(stdict['EPOCH_MICROSECONDS'])*u.microsecond
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

def stmetadata(stdata):
    return {'name': stdata['OBJECT_NAME'],
            'type': stdata['OBJECT_TYPE'],
            'catid': int(stdata['OBJECT_NUMBER']),
            'intldes': stdata['OBJECT_ID']}

def spacetrack_latest(stclient, satnums):
    stdata = stclient.tle_latest(norad_cat_id=satnums, ordinal=1)
    ret = [(satdata(std), stmetadata(std)) for std in stdata]
    if type(satnums) is int:
        return ret[0]
    else:
        return ret
