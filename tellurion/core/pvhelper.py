import numpy as np
import astropy.units as u
import astropy.time
from tellurion.astro import units
from tellurion.astro import quant

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

def ncartesianpv(pv, given=units.prefunits, convert=units.prefunits):
    """Create a u.Quantity position or position-velocity array or tuple of arrays."""
    if type(pv) is tuple:
        tup = pv
    else:
        arr = np.array(pv)
        if arr.shape[arr.ndim-1] == 6:
            if arr.ndim == 2:
                tup = (arr[:,0:3], arr[:,3:6])
            else:
                tup = (arr[0:3], arr[3:6])
    pos = tup[0]*given['length']
    if len(tup) == 2:
        vel = tup[1]*given['speed']
        ou = quant.hstack((pos.structure(_eph_pos, pos.ndim==1), \
                           vel.structure(_eph_vel, vel.ndim==1)))
        if convert:
            return quant.nchangeunits(ou, convert)
        else:
            return ou
    else:
        if convert:
            return quant.nchangeunits(pos, convert)
        else:
            return pos

# pos1 = cartesianpv(np.array([1,2,3]), ['m'])
# pv1 = cartesianpv(np.array([1,2,3,4,5,6]), ['m', 'm/s'])
# pos2 = cartesianpv(np.array([[1,2,3], [-1,-2,-3], [-3,-2,-1]]), ['m'])
# pv2 = cartesianpv([[1,2,3,4,5,6], [-1,-2,-3,-4,-5,-6]], ['m', 'm/s'])
def cartesianpv(array, units=('length', 'speed'), unitlookup=units.prefunits):
    """Make a Cartesian position or position-velocity structured
    quantity from the array (or something that can be converted to an
    array with np.array). May be a scalar (1d array) or vector (2d
    array) structure.

    If units is not None, it should be a tuple or list of length two
    with units for length and speed.

    """
    if type(array) is tuple:
        arr = np.array((np.hstack(array)))
    else:
        arr = np.array(array)
    if arr.shape[arr.ndim-1] == 6:
        nsu = [(_eph_pos, 3, units[0]), (_eph_vel, 3, units[1])]
        return quant.sq(arr, nsu, unitlookup=unitlookup)
    elif arr.shape[arr.ndim-1] == 3:
        nsu = [(_eph_pos, 3, units[0])]
        sq = quant.sq(arr, nsu, unitlookup=unitlookup)
        return sq[_eph_pos]
    else:
        raise ValueError(f"Incorrect shape of array to make posvel {arr.shape}")

def cartesianpv_sep(position, velocity):
    """Make a Cartesian position-velocity from separate position and velocity quantities."""
    if type(position) is u.Quantity and type(velocity) is u.Quantity:
      posdec = position.decompose()
      veldec = velocity.decompose()
      dval = [(_eph_pos, '<f8', (3,)), (_eph_vel, '<f8', (3,))]
      return quant.compose_sq(np.hstack((posdec[0], veldec[0])), dval, (posdec[2], veldec[2]))
    else:
        raise ValueError(f"Position and velocity {position, velocity} must be of type u.Quantity")

def sphericalpv(sph_position, sph_velocity=None, labels=['rtasc','decl','distance'], unitlookup=units.prefunits):
    """Make a spherical coordinate set for position and velocity

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
    """
    if sph_velocity:
        # Include velocity components
        posdict = {label: u.Quantity(value) for label, value in zip(labels, sph_position)}
        labels_r = [sym+'_r' for sym in labels]
        veldict = {label: u.Quantity(value) for label, value in zip(labels_r, sph_velocity)}
        return quant.sq_from_dict(posdict | veldict)
    else:
        # Position only
        posdict = {label: u.Quantity(value) for label, value in zip(labels, sph_position)}
        return quant.sq_from_dict(posdict)

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
