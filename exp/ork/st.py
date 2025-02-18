import os
import astropy.units as u
import spacetrack
import astro
import posvel
import org.orekit.propagation.analytical.tle as tle

# Set stclient using registered space-track.org username and password
# stclient = spacetrack.SpaceTrackClient(identity="myemail@example.com", password="mypw")

# from experimental import *
# from ork.st import *
# isssent = spacetrack_latest(stclient, [25544, 41335])
# (iss, sentinel3a) = (isssent[0], isssent[1])
# isskep = astro.makesq(iss[0][0])
# NOT AN ACCURATE COMPUTATION OF CARTESIAN POSITION, IT ASSUMES KEPLER ELEMENTS ARE OSCULATING:
# iss_badpvt = ork.element.keplerianorbit(isskep, iss[0][1]).pvt()
#
# A better choice would be to use Orekit to propagate/convert
# isstle = tle_latest(stclient, 25544)

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

# Uses Orekit
def tle_latest(stclient, satnum):
    sattle = stclient.tle_latest(norad_cat_id=satnum, ordinal=1, format='tle')
    sat2lines = sattle.split("\n")
    return tle.TLE(sat2lines[0], sat2lines[1])
