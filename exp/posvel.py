"""
Position, velocity sets in AstroPy
"""
import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries
import astropy.table.row
# from . import astro
import astro
import dttm

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]

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

# This will be needed for makept
def makepos(pos, unit=astro.prefunits['length']):
    '''Create a position vector or convert units'''
    if u.get_physical_type(unit) == 'length':
        if isq3vec(pos, 'length'):
            return pos.to(unit)
        elif type(pos) is u.Quantity:
            raise ValueError('Argument does not represent a position 3-vector')
        else:
            return pos*unit
    else:
        raise ValueError('Unit does not represent a position')

def makepv(pos, vel, units=astro.prefunits["posvel"]):
    '''Make a posvel from separate position and velocity; if argument `pos` is a posvel, then convert units'''
    # See https://docs.astropy.org/en/stable/units/structured_units.html#example
    if ispv(pos):
        return(pos.to(units))
    if isq3vec(pos, 'length') and isq3vec(vel, 'speed'):
        return(makepv(pos.value, vel.value, u.StructuredUnit((pos.unit, vel.unit))).to(units))
    else:
        pvtype = [(_eph_pos, '(3,)f8'), (_eph_vel, '(3,)f8')]
        npa = np.array((pos, vel), dtype = pvtype)
        pv = u.Quantity(npa, u.StructuredUnit(units))
        return(pv)

# A PVT consists of a tuple a posvel (as defined by ispv()) and an astropy.time.Time
def makepvt(obj):
    '''Return a tuple of posvel and time, from either an ephrow or a (pos, vel, time) tuple.'''
    if isephrow(obj):
        pos = obj[_eph_pos]
        vel = obj[_eph_vel]
        return (makepv(pos, vel), obj[_eph_time])
    elif ispv(obj):
        return (obj, dttm.nowutc())
    elif type(obj) is tuple and len(obj)==3:
        return (makepv(obj[0], obj[1]), obj[2])
    else:
        raise ValueError('Cannot make a PVT from this object')
