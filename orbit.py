from orkinit import *
from dttm import *
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory
# ECIS
from ecis import *

# Orekit configuration
okc = {'cartesian': OrbitType.CARTESIAN,
       'meananom': PositionAngleType.MEAN,
       'trueanom': PositionAngleType.TRUE}

# Orbital environment constants
envct = {
        'earthframe': FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        'celestialframe': FramesFactory.getEME2000(),
        'earthrad': Constants.IERS2010_EARTH_EQUATORIAL_RADIUS,
        'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY,
        'earthJ2': -Constants.IERS2010_EARTH_C20,
        'earthmu': Constants.IERS2010_EARTH_MU,
        'earthflat': Constants.IERS2010_EARTH_FLATTENING,
        'sun': CelestialBodyFactory.getSun(),
}

######## Convert between Python and Orekit PVT

# PVT as Orekit arrays and AbsoluteDate
# pos, vel, dttm: Python position 3-vector, velocity 3-vector, and datetime
def pvtok(pos, vel, dttm):
    return(TimeStampedPVCoordinates(datetime_to_absolutedate(dttm), Vector3D(pos), Vector3D(vel)))

# PVT as Python arrays and datetime
TimeStampedPVCoordinates.pyrep = lambda self: [[self.position.x, self.position.y, self.position.z],
                                               [self.velocity.x, self.velocity.y, self.velocity.z],
                                               absolutedate_to_datetime(self.date)]

######## Make orbits

# Make a Cartesian orbit from Orekit PVT
TimeStampedPVCoordinates.cartesian = lambda self: CartesianOrbit(self, envct['celestialframe'], envct['earthmu'])

# Make a Kepler orbital element set - this should be
def kepler(epoch, sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, timeelt_type):
    return(KeplerianOrbit(sma, # Semimajor Axis (m)
                          ecc,    # Eccentricity
                          radians(inc_deg),  # Inclination (deg)
                          radians(argper_deg),   # Perigee argument (deg)
                          radians(raan_deg),   # Right ascension of ascending node (degrees)
                          radians(timeelt_deg),  # Time element (deg)
                          timeelt_type,  # Sets which type of anomaly we use
                          envct['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                          epoch,   # Sets the date of the orbital parameters
                          envct['earthmu']))   # Sets the central attraction coefficient (m³/s²)

######## Convert orbits

# Convert to the requested orbit type "cart", "kep"
def convert(tree, orbtype):
    if "cart" in tree:
        orbit = tree.cart
    elif "kep" in tree:
        orbit = tree.kep
    elif "circ" in tree:
        orbit = tree.circ
    elif "equi" in tree:
        orbit = tree.equi
    elif "pvt" in tree:
        tree.update(cart = tree.pvt.cartesian())
        orbit = tree.cart
    else:
        raise ValueError("Nothing to convert")
    match orbtype:
        case "cart":
            ret = OrbitType.CARTESIAN.convertType(orbit)
            tree.update(cart = ret)
        case "kep":
            ret = OrbitType.KEPLERIAN.convertType(orbit)
            tree.update(kep = ret)
        case "circ":
            ret = OrbitType.CIRCULAR.convertType(orbit)
            tree.update(circ = ret)
        case "equi":
            ret = OrbitType.EQUINOCTIAL.convertType(orbit)
            tree.update(equi = ret)
        case "pvt":
            ret = orbit.pVCoordinates
            tree.update(pvt = ret)
        case _:
            raise ValueError("Type \"" + orbtype + "\" unknown")
    return(ret)

Ecis.cartesian = lambda tree: convert(tree, "cart")
Ecis.kepler = lambda tree: convert(tree, "kep")
Ecis.circular = lambda tree: convert(tree, "circ")
Ecis.equinoctial = lambda tree: convert(tree, "equi")
Ecis.pvt = lambda tree: convert(tree, "pvt")

######## Properties of orbits

# The position-velocity-time for the state
# The inverse of cartesian()
def pvt(orbit):
    return(orbit.pVCoordinates)

# The geocentric distance of the orbit
def posmag(orbit):
    return(orbit.pVCoordinates.position.norm)

def period(orbit):
    return(orbit.getKeplerianPeriod())

## Make a computation tree from PVT

def cspvt(pos, vel, dttm):
    return(newtree('pvt', pvtok(pos, vel, dttm)))

## Example orbit
ex1 = cspvt([5740132.68349499, 3314067.15, 0.0],
            [-2750.82683526322, 4764.5718414998, 5501.65367052644],
            datetime(2022, 6, 1, 12, 0, 0))

# ex1.cartesian()

# Need to build and convert a Kepler
