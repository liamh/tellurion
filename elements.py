from propagate import *
from org.orekit.orbits import KeplerianOrbit, PositionAngleType

okc['meananom']= PositionAngleType.MEAN
okc['trueanom']= PositionAngleType.TRUE

# Make a Kepler orbital element set
def kepler_oes(epoch, sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, timeelt_type):
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

def period(orbit):
    return(convert(orbit,"kep").getKeplerianPeriod())
