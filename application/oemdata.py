import numpy as np
import astropy.units as u
from oem import OrbitEphemerisMessage
import tellurion as tell

oemunits = {"time": u.second, "length": u.km, "speed": u.km/u.second,
            "angle": u.degree, "angular speed": u.radian/u.second,
            "dimensionless": u.dimensionless_unscaled}

# simorb = oempvt("../../../gmat/EphemerisFile1.oem")
def oempvt(file):
    '''Read the CCSDS Orbital Ephemeris Message data file and create a PositionVelocityT from the data.'''
    ephemeris = OrbitEphemerisMessage.open(file)
    array = np.vstack(tuple([np.hstack((row.position, row.velocity)) \
                             for row in ephemeris.states]))
    pv = tell.ncartesianpv(array, oemunits)
    times = tell.abstime([s.epoch for s in ephemeris.states])
    return tell.PositionVelocityT(time=times, cartesian=pv)
