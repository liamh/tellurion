"""
Position, velocity and time sets in AstroPy
"""
import datetime
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.timeseries import TimeSeries
import astropy.table.row
import astropy.coordinates as coord
# from . import astro
import astro
import cdttm

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]

_pos_format = '10.3f'
_vel_format = '10.6f'

##################################################
####   Tests for posvel and related types     ####
##################################################

def isq3vec(obj, physdim):
    '''Is a 3-vector u.Quantity with the specified physical dimension'''
    return type(obj) is u.Quantity \
        and u.get_physical_type(obj) == u.get_physical_type(physdim) \
        and np.size(obj)==3

def ispv(obj):
    return type(obj) == u.Quantity and obj.dtype.names is not None \
        and _eph_pos in obj.dtype.names and _eph_vel in obj.dtype.names \
        and isq3vec(obj[_eph_pos],'length') and isq3vec(obj[_eph_vel],'speed')

def isdttm(obj):
    return type(obj) is astropy.time.Time

def isreltime(obj):
    return type(obj) is u.Quantity and u.get_physical_type(obj) == u.get_physical_type('time')

def ispvt(obj):
    return type(obj) == tuple and len(obj) == 2 \
        and ispv(obj[0]) and isdttm(obj[1])

def isephem(ts):
    '''Argument is an ephemeris table'''
    return type(ts) is TimeSeries \
        and all([k in ts.keys() for k in _ephemeris_columns])

def isephrow(row):
    '''Argument is a row of an ephemeris table'''
    return type(row) is astropy.table.row.Row \
        and all([row.keys().__contains__(k) for k in _ephemeris_columns])

##################################################
####   Make posvel and related types          ####
##################################################

def pv(position, velocity, unitlookup=astro.prefunits):
    '''Make a posvel from separate position and velocity; if argument `position` is a posvel, then convert units'''
    if ispv(position):
        return(position.to(unitlookup['length']))
    else:
        return astro.makesq((position, velocity), (_eph_pos, _eph_vel), \
                            phystype=('length','speed'), unitlookup=unitlookup)

##################################################
#### Dates, times, and PVT                    ####
##################################################

def nowutc():
    """
    The time now (in UTC) as dttm.
    """
    return astropy.time.Time(datetime.datetime.now(datetime.UTC), scale='utc')

# To add timezone to datetime
#import pytz
#def utcdt(datetime):
#    return pytz.utc.localize(datetime)

# A PVT consists of a tuple a posvel (as defined by ispv()) and an astropy.time.Time
def pvt(obj, item=None):
    '''Return a tuple of posvel and time from a variety of sources.
    '''
    if isephrow(obj):
        # PVT from an ephemeris row
        pos = obj[_eph_pos]
        vel = obj[_eph_vel]
        return (pv(pos, vel), obj[_eph_time])
    elif type(obj) is TimeSeries:
        # Select a row from an ephemeris by index, absolute time, or relative time
        if item == None:   # Return the last row
            item = -1
        try:
            row = obj[item]
        except:
            try:
                row = obj.loc[cdttm.dttm(item)]
            except:
                row = obj.loc[obj[0]['time'] + item]
        return pvt(row)
    elif ispv(obj):
        if item==None:
            # Add the current time to the PV
            return (obj, nowutc())
        else:
            # Add the specified time to the PV
            return (obj, cdttm.dttm(item))
    elif ispvt(obj):
        if isdttm(item):
            # Replace the timestamp in the PVT
            return (obj[0], item)
        else:
            # Displace the timestamp in the PVT by the given relative time
            return (obj[0], obj[1] + item)
    elif type(obj) is tuple:
        # Create a PVT from the three P, V, T
        return (pv(obj[0], obj[1]), obj[2])
    else:
        raise ValueError('Cannot make a PVT from this object')

def makepos(pos, unit=astro.prefunits['length']):
    '''Create a position vector or convert units'''
    if type(pos) is coord.SkyCoord:
        return makepos(pos.cartesian, unit)
    if type(pos) is coord.representation.cartesian.CartesianRepresentation:
        return makepos(pos.xyz, unit)
    if type(pos) is tuple:
        return u.Quantity(pos, unit)
    if u.get_physical_type(unit) == 'length':
        if isq3vec(pos, 'length'):
            return pos.to(unit)
        elif type(pos) is u.Quantity:
            raise ValueError('Argument does not represent a position 3-vector')
        else:
            return pos*unit
    else:
        raise ValueError('Unit does not represent a position')

def makept(pos, dttm, frame='gcrs'):
    '''Make a SkyCoord GCRS location'''
    return coord.SkyCoord(coord.CartesianRepresentation(pos), obstime=dttm, frame=frame)

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    '''Magnitude of the difference of two vectors'''
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])
