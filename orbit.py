######### Representation of state as Cartesian, Kepleran, Circular, Equinoctial

## A position-velocity-time PVT is a Cartesian state, but
## Orekit has a separate class of CartesianOrbit that can be
## propagated, so these are have different representation.

from orkinit import *
from dttm import *
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory
# ECIS - exploratory computation in stages
from ecis import *

# Orekit configuration
okc = {'cartesian': OrbitType.CARTESIAN}

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

######## Snaglab (Python) and Orekit representation of a postion-velocity-time (PVT)

# PVT as Orekit arrays and AbsoluteDate
# These are made with the `ork` variable set to a TimeStampedPVCoordinates instance
class PVT:
    position: list # Position 3-vector in meters
    velocity: list # Velocity 3-vector in meters/seconds
    dttm: datetime # Epoch time
    ork: TimeStampedPVCoordinates

    def __init__(self, position, velocity, dttm, ork=None):
        self.position = position
        self.velocity = velocity
        self.dttm = dttm
        if ork is None:
            self.ork = TimeStampedPVCoordinates(datetime_to_absolutedate(self.dttm),
                                                Vector3D(self.position),
                                                Vector3D(self.velocity))
        else:
            self.ork = ork
    def __repr__(self):
        return f"<PVT position: {self.position} (m) velocity:{self.velocity} (m/s) epoch {self.dttm} (UTC)>"

# .snl(): Convert PVT from Orekit to Python
TimeStampedPVCoordinates.snl \
    = lambda self: PVT([self.position.x, self.position.y, self.position.z],
                       [self.velocity.x, self.velocity.y, self.velocity.z],
                       absolutedate_to_datetime(self.date),
                       self)

######## Make orbits

# Make a Cartesian orbit from Orekit PVT
TimeStampedPVCoordinates.cartesian = lambda self: CartesianOrbit(self, envct['celestialframe'], envct['earthmu'])

# Make a Kepler orbital element set
def kepler(epoch, sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, mean_timeelt):
    if mean_timeelt:
        timeelt_type = PositionAngleType.MEAN
    else:
        timeelt_type = PositionAngleType.TRUE
    return(KeplerianOrbit(sma, # Semimajor Axis (m)
                          ecc,    # Eccentricity
                          radians(inc_deg),  # Inclination (deg)
                          radians(argper_deg),   # Perigee argument (deg)
                          radians(raan_deg),   # Right ascension of ascending node (degrees)
                          radians(timeelt_deg),  # Time element (deg)
                          timeelt_type,  # Sets which type of anomaly we use
                          envct['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                          datetime_to_absolutedate(epoch),   # Sets the date of the orbital parameters
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
        tree.update(cart = tree.pvt.ork.cartesian())
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
            ret = orbit.posveltime()
            tree.update(pvt = ret)
        case _:
            raise ValueError("Type \"" + orbtype + "\" unknown")
    return(ret)

Ecis.cartesian = lambda tree: convert(tree, "cart")
Ecis.kepler = lambda tree: convert(tree, "kep")
Ecis.circular = lambda tree: convert(tree, "circ")
Ecis.equinoctial = lambda tree: convert(tree, "equi")
Ecis.posveltime = lambda tree: convert(tree, "pvt")

######## Properties of orbits

# The position-velocity-time for the state
# The inverse of .cartesian()
Orbit.posveltime = lambda orbit: orbit.pVCoordinates.snl()

# The geocentric distance of the orbit
def posmag(orbit):
    return(orbit.pVCoordinates.position.norm)

def period(orbit):
    return(orbit.getKeplerianPeriod())

## Make a computation tree from position, velocity, and datetime

def new_posveltime(pos, vel, dttm):
    ret = newtree('pvt', PVT(pos, vel, dttm)) # Create the tree and set the first component to the PVT
    ret.cartesian() # Convert the PVT to the Orekit CartesianOrbit and save that as the next component
    return(ret)

## Make a computation tree from Kepler elements and datetime

def new_kepler(sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, mean_timeelt, epoch):
    ret = newtree('kep', kepler(epoch, sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, mean_timeelt))
    convert(ret, "pvt")
    return(ret)

## Example orbit
ex1 = new_posveltime([5740132.68349499, 3314067.15, 0.0],
            [-2750.82683526322, 4764.5718414998, 5501.65367052644],
            datetime(2022, 6, 1, 12, 0, 0))
ex1.pvtorbrec = ex1.cart.posveltime() # The PVT recalculated from the Cartesian orbit
ex1.pvtorkrec = ex1.pvt.ork.snl() # The PVT recalculated from the Orekit representation
# ex1.keys()
# ex1.cartesian()

# Need to build and convert a Kepler
ex2 = new_kepler(8.0e6, 0.1, 42.0, 217.4, -90.0, 7.25, True, datetime(2023, 9, 14, 8, 30, 0))
