import numbers
import numpy as np
import astropy.units as u
import astropy.time
from tellurion.astro import units
from tellurion.astro import quantity_utils as quant

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = 'time'
_eph_pos = 'position'
_eph_vel = 'velocity'
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]
_ephemeris_columns_pos_xyz = [_eph_time, 'px','py','pz']
_ephemeris_columns_pos_only = [_eph_time, _eph_pos]
_posvel_pt = {_eph_pos: 'length', _eph_vel: 'speed'}

# These should be conditional on the units used
_pos_format = '10.3f'
_vel_format = '10.6f'

##################################################
####   Utility functions for coordinates      ####
##################################################

def cartesianpv(posvel, isscalar, given=units.prefunits, convert=units.prefunits):
    """Create a u.Quantity position or position-velocity scalar or
    array. Inputs `posvel` may be a tuple of (position, velocity) or
    (position), with values that are u.Quantity, np.array with final
    dimension 6 (for position and velocity) or 3 (for position)."""
    scalar = isscalar
    pvunit = _posvel_pt
    if type(posvel) is tuple:
        if len(posvel) == 2:
            pv = {_eph_pos: posvel[0], _eph_vel: posvel[1]}
        else:
            pv = {_eph_pos: posvel[0]}
    elif type(posvel) is np.ndarray:
        if posvel.shape[posvel.ndim-1] == 6:
            if posvel.ndim == 2:
                pv = {_eph_pos: posvel[:,0:3], _eph_vel: posvel[:,3:6]}
                scalar = False
            else:
                pv = {_eph_pos: posvel[0:3], _eph_vel: posvel[3:6]}
                scalar = True
        elif posvel.shape[posvel.ndim-1] == 3:
            pvunit = _posvel_pt[_eph_pos]
            if posvel.ndim == 2:
                pv = posvel[:,0:3]
                scalar = False
            else:
                pv = posvel[0:3]
                scalar = True
    else:
        pv = posvel
    if convert and given != convert:
        return quant.change_units(quant.make_quantity(pv, pvunit, scalar, given), convert)
    else:
        return quant.make_quantity(pv, pvunit, scalar, given)

sphpospt = ['angle', 'angle', 'length']
sphvelpt = ['angular speed', 'angular speed', 'speed']

def sphericalpv(sph_position, sph_velocity=None, labels=['rtasc','decl','distance'], unitlookup=units.prefunits):
    """Make a spherical coordinate set for position and velocity

    Args:
        sph_position: [cylang, sphang, distance] where each can be:
            - Scalar number or u.Quantity
            - List of numbers or u.Quantity
            - np.array
        sph_velocity: Time derivatives of sph_position (optional, same format)
        labels: Field names for the coordinates
        unitlookup: Unit lookup dictionary

    Returns:
        Structured quantity with spherical coordinates.
    """
    isscalar = not(hasattr(sph_position[0], '__getitem__'))
    posdict = dict(zip(labels, sph_position))
    if sph_velocity:
        # Include velocity components
        labels_r = [sym+'_r' for sym in labels]
        veldict = dict(zip(labels_r, sph_velocity))
        unt = dict(zip(labels+labels_r, sphpospt+sphvelpt))
        return quant.make_quantity(posdict | veldict, unt, isscalar, unitlookup)
    else:
        # Position only
        unt = dict(zip(labels, sphpospt))
        return quant.make_quantity(posdict, unt, isscalar, unitlookup)

##################################################
####   Tests for posvel and related types     ####
##################################################

def isq3vec(obj, physdim):
    """Is a 3-vector u.Quantity with the specified physical dimension"""
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
    """Argument is an ephemeris table"""
    return type(ts) is TimeSeries \
        and all([k in ts.keys() for k in _ephemeris_columns])

def isephrow(row):
    """Argument is a row of an ephemeris table"""
    return type(row) is astropy.table.row.Row \
        and all([row.keys().__contains__(k) for k in _ephemeris_columns])

##################################################
#### Compare positions, velocities            ####
##################################################

def magdiff(a, b):
    """Magnitude of the difference of two vectors"""
    return u.Quantity([np.linalg.norm(ai - bi) for (ai, bi) in zip(a, b)])

def posdiff(a, b):
    return np.linalg.norm(a[_eph_pos]-b[_eph_pos])
