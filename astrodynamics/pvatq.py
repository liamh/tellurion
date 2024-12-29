"""
Vectors and PVT (position, velocity, time) sets in AstroPy

No definitions for direct use
"""
import numpy as np
import astropy.units as u
from . import astro
from . import util

# posvel preferred and SI units (for Orekit)
astro.prefunits["posvel"] = (astro.prefunits["length"], astro.prefunits["velocity"])
posvelsiu = (u.meter, u.meter/u.second)

################################################################################
#### PVT: Position, velocity, and time
################################################################################

def atqpvt(fromthing, time=None, units=astro.prefunits["posvel"]):
    '''Create astro.TQuantity of posvel from many possible different inputs'''

    ### Subfunctions needed to sort out the input "fromthing"
    def isnppvt(obj):
        '''Object is a length-2 tuple with posvel as numpy array and time as datetime64'''
        return(type(obj) is tuple and len(obj)==2 \
               and type(obj[0]) is np.ndarray and type(obj[1]) is np.datetime64)

    def isatqpvt(obj):
        '''Object is an atq position possibly with time'''
        isatq = type(obj) == astro.TQuantity or type(obj) == u.Quantity
        if not isatq:
            return False
        return 'p' in obj.dtype.names and 'v' in obj.dtype.names

    atqargunits = None

    # Convert all possible inputs into standardized pos, vel, time
    if isnppvt(fromthing): # An np tuple (pv, time)
        ftis = "np"
        (pv, time) = fromthing
        pos=pv[0:3]
        vel=pv[3:6]
    elif util.listnpa(fromthing): # Either python list or numpy array without time
        (pos, vel) = pvsplit(fromthing)
    elif isatqpvt(fromthing):  # Already an atq, but may need unit conversion
        atqargunits = fromthing.convert_units(units)
        time = fromthing.time.datetime
    elif type(fromthing) == astropy.table.row.Row:
        pos = fromthing['position'].convert_units(units[0])
        vel = fromthing['velocity'].convert_units(units[1])
        time = fromthing['time']
    else:
        warnings.warn(f"Cannot convert object `{fromthing}` to a PVT")
        return None

    # Set the TQuantity variable
    if atqargunits is None:
        return(atqptpvt(pos, vel, time, units=units))
    else:
        return(atqargunits)

def pvsplit(pv):
    '''Create a tuple of position and (optionally) velocity 3-vectors (or lists)'''
    match len(pv):
        case 6:            # Position and velocity
            pos = pv[0:3]
            vel = pv[3:6]
        case 2:            # Position and velocity
            pos = pv[0]
            vel = pv[1]
        case 3:            # Position only
            pos = pv[0:3]
            vel = None
        case 1:            # Position only
            pos = pv[0]
            vel = None
    return (pos,vel)

def atqptpvt(pos, vel, time=None, units=astro.prefunits["posvel"]):
    '''Create a PT or PVT as an astro.TQuantity'''
    if np.isnan(vel[0]):
        pvtype = [('p', np.float64)]
        npa = np.array(pos, dtype = pvtype)
        pv = astro.TQuantity(npa, u.StructuredUnit(units[0]), time)
    else:
        # See https://docs.astropy.org/en/stable/units/structured_units.html#example
        pvtype = [('p', '(3,)f8'), ('v', '(3,)f8')]
        npa = np.array((pos, vel), dtype = pvtype)
        pv = astro.TQuantity(npa, u.StructuredUnit(units), time)
    return(pv)

# Scale position and velocity separately
# pv = posvel to scale
# pvscale = list or vector of length 2 to scale position, velocity
def _scale_posvel(pv, pvscale):
    return posvel([pv.value[0]*pvscale[0], pv.value[1]*pvscale[1]], time=pv.time,
                  length_unit=pv.unit[0], velocity_unit=pv.unit[1])

# Convert posvel units
# ex1pvtsi = ex1pvt.convert_units((u.m, u.m/u.s))
