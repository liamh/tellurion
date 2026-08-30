import numpy as np

from tellurion.astro import quantity_utils as quant, units

##################################################
####   Constants used to define field names   ####
##################################################

_eph_time = "time"
_eph_pos = "position"
_eph_vel = "velocity"
_ephemeris_columns = [_eph_time, _eph_pos, _eph_vel]
_ephemeris_columns_pos_only = [_eph_time, _eph_pos]
_posvel_pt = {_eph_pos: "length", _eph_vel: "speed"}
_sphpospt = ["angle", "angle", "length"]
_sphvelpt = ["angular speed", "angular speed", "speed"]



# These should be conditional on the units used
_pos_format = "10.3f"
_vel_format = "10.6f"

##################################################
####   Utility functions for coordinates      ####
##################################################

def _cartesianpv(posvel, isscalar, given=units.prefunits, convert=units.prefunits):
    """Create a u.Quantity position or position-velocity scalar or array."""
    scalar = isscalar
    pvunit = _posvel_pt

    if type(posvel) is tuple:
        scalar = len(posvel[0].shape)==1
        if len(posvel) == 2:
            pv = {_eph_pos: posvel[0], _eph_vel: posvel[1]}
        else:
            pv = {_eph_pos: posvel[0]}
    else:
        # Convert to numpy array (handles list, ndarray, and other sequences)
        arr = np.asarray(posvel)
        width = arr.shape[-1] if arr.ndim > 0 else 0

        if width >= 6:
            if arr.ndim == 2:
                pv = {_eph_pos: arr[:,0:3], _eph_vel: arr[:,3:6]}
                scalar = False
            else:
                pv = {_eph_pos: arr[0:3], _eph_vel: arr[3:6]}
                scalar = True
        elif width == 3 or width == 4:
            pvunit = _posvel_pt[_eph_pos]
            pv = {_eph_pos: arr}
            scalar = arr.ndim == 1

    return quant.change_units(quant.make_quantity(pv, pvunit, scalar, given), convert)


def _sphericalpv(sph_position, sph_velocity=None, labels=["rtasc","decl","distance"],
                 unitlookup=units.prefunits):
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
    isscalar = not(hasattr(sph_position[0], "__getitem__"))
    posdict = dict(zip(labels, sph_position))
    if sph_velocity:
        # Include velocity components
        labels_r = [sym+"_r" for sym in labels]
        veldict = dict(zip(labels_r, sph_velocity))
        unt = dict(zip(labels+labels_r, _sphpospt+_sphvelpt))
        return quant.make_quantity(posdict | veldict, unt, isscalar, unitlookup)
    else:
        # Position only
        unt = dict(zip(labels, _sphpospt))
        return quant.make_quantity(posdict, unt, isscalar, unitlookup)
