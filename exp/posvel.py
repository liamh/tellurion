"""
Position, velocity sets in AstroPy
"""
import numpy as np
import astropy.units as u
# from . import astro
import astro

def isq3vec(obj, physdim):
    '''Is a 3-vector u.Quantity with the specified physical dimension'''
    return type(obj) is u.Quantity \
        and u.get_physical_type(obj) == u.get_physical_type(physdim) \
        and np.size(obj)==3

def ispv(obj):
    return type(obj) == u.Quantity and obj.dtype.names is not None \
        and 'p' in obj.dtype.names and 'v' in obj.dtype.names \
        and isq3vec(obj['p'],'length') and isq3vec(obj['v'],'speed')

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
        pvtype = [('p', '(3,)f8'), ('v', '(3,)f8')]
        npa = np.array((pos, vel), dtype = pvtype)
        pv = u.Quantity(npa, u.StructuredUnit(units))
        return(pv)
