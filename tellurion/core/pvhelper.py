import numpy as np
import astropy.units as u
import astropy.time
from tellurion.core import astro
from tellurion.core import nquant

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]
_ephemeris_columns_pos_xyz = [_eph_time, 'px','py','pz']
_ephemeris_columns_pos_only = [_eph_time, _eph_pos]

# These should be conditional on the units used
_pos_format = '10.3f'
_vel_format = '10.6f'

##################################################
####   Utility functions for coordinates      ####
##################################################

# pos1 = cartesianpv(np.array([1,2,3]), ['m'])
# pv1 = cartesianpv(np.array([1,2,3,4,5,6]), ['m', 'm/s'])
# pos2 = cartesianpv(np.array([[1,2,3], [-1,-2,-3], [-3,-2,-1]]), ['m'])
# pv2 = cartesianpv([[1,2,3,4,5,6], [-1,-2,-3,-4,-5,-6]], ['m', 'm/s'])
def cartesianpv(array, units=['length', 'speed'], unitlookup=astro.prefunits):
    '''Make a Cartesian position or position-velocity structured
    quantity from the array (or something that can be converted to an
    array with np.array). May be a scalar (1d array) or vector (2d
    array) structure.
    '''
    arr = np.array(array)
    if arr.shape[arr.ndim-1] == 6:
        nsu = [(_eph_pos, 3, units[0]), (_eph_vel, 3, units[1])]
        velp = True
    elif arr.shape[arr.ndim-1] == 3:
        nsu = [(_eph_pos, 3, units[0])]
        velp = False
    else:
        raise ValueError(f"Incorrect shape of array to make posvel {arr.shape}")
    return (nquant.sq(arr, nsu, unitlookup=unitlookup), velp)


def sphericalpv(sph_position, sph_velocity=None, labels=['rtasc','decl','distance'], unitlookup=astro.prefunits):
    '''Make a spherical coordinate set for position and velocity

    Args:
        sph_position: [cylang, sphang, distance] where each can be:
            - Scalar number or u.Quantity
            - List of numbers or u.Quantity
            - 1D or 2D np.ndarray
        sph_velocity: Time derivatives of sph_position (optional, same format)
        labels: Field names for the coordinates
        unitlookup: Unit lookup dictionary

    Returns:
        Structured quantity with spherical coordinates.
    '''
    if sph_velocity:
        # Include velocity components
        labels_r = [sym+'_r' for sym in labels]
        return nquant.structquant(
            sph_position + sph_velocity,
            labels + labels_r,
            units=['angle', 'angle', 'length',
                   'angular speed', 'angular speed', 'speed'],
            unitlookup=unitlookup
        )
    else:
        # Position only
        return nquant.structquant(
            sph_position,
            labels,
            units=['angle', 'angle', 'length'],
            unitlookup=unitlookup
        )

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
    """Object is a relative time: is a u.Quantity with physical type 'time'"""
    return type(obj) is u.Quantity and u.get_physical_type(obj) == u.get_physical_type('time')

def isephem(ts):
    '''Argument is an ephemeris table'''
    return type(ts) is TimeSeries \
        and all([k in ts.keys() for k in _ephemeris_columns])

def isephrow(row):
    '''Argument is a row of an ephemeris table'''
    return type(row) is astropy.table.row.Row \
        and all([row.keys().__contains__(k) for k in _ephemeris_columns])

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    '''Magnitude of the difference of two vectors'''
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])

def posdiff(a, b):
    return np.linalg.norm(a[_eph_pos]-b[_eph_pos])
