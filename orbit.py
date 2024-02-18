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
# pvt: Orekit PVT
def pvtpy(pvt):
    jp = [pvt.position.x, pvt.position.y, pvt.position.z]
    jv = [pvt.velocity.x, pvt.velocity.y, pvt.velocity.z]
    epoch = absolutedate_to_datetime(pvt.date)
    return([jp, jv, epoch])

######## Make orbits

# Make a Cartesian orbit from the Orekit PVT
# The inverse of pvt()
def cartesian(pvt):
    return(CartesianOrbit(pvt, envct['celestialframe'], envct['earthmu']))

# Make a Kepler orbital element set
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

## Example orbit
# ex1p = [5740132.68349499, 3314067.15, 0.0]
# ex1v = [-2750.82683526322, 4764.5718414998, 5501.65367052644]
# ex1t = datetime(2022, 6, 1, 12, 0, 0)
# ex1pvt = pvtok(ex1p,ex1v,ex1t)
# ex1orb = cartesian(ex1pvt)

######## Convert orbits

# Convert to the requested orbit type "cart", "kep"
def convert(orbit, orbtype):
    match orbtype:
        case "cart":
            return(OrbitType.CARTESIAN.convertType(orbit))
        case "kep":
            return(OrbitType.KEPLERIAN.convertType(orbit))
        case "circ":
            return(OrbitType.CIRCULAR.convertType(orbit))
        case "equi":
            return(OrbitType.EQUINOCTIAL.convertType(orbit))
        case _:
            raise ValueError("Type \"" + orbtype + "\" unknown")

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
