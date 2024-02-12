from orkinit import *
from dttm import *
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit

okc['meananom']= PositionAngleType.MEAN
okc['trueanom']= PositionAngleType.TRUE

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

## TODO
# To be consistent, datm() should generate a Python datetime so that pvtok can convert.
# ex1p = [5740132.68349499, 3314067.15, 0.0]
# ex1v = [-2750.82683526322, 4764.5718414998, 5501.65367052644]
# ex1t = datm(2022, 6, 1, 12, 0, 00.000)

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
