from setup import *

meananom = PositionAngleType.MEAN
trueanom = PositionAngleType.TRUE

# Make a Kepler orbital element set
def kepler_oes(epoch, sma, ecc, inc_deg, raan_deg, argper_deg, timeelt_deg, timeelt_type):
    return(KeplerianOrbit(sma, # Semimajor Axis (m)
                          ecc,    # Eccentricity
                          radians(inc_deg),  # Inclination (deg)
                          radians(argper_deg),   # Perigee argument (deg)
                          radians(raan_deg),   # Right ascension of ascending node (degrees)
                          radians(timeelt_deg),  # Time element (deg)
                          timeelt_type,  # Sets which type of anomaly we use
                          FramesFactory.getEME2000(), # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                          epoch,   # Sets the date of the orbital parameters
                          earthmu))   # Sets the central attraction coefficient (m³/s²)

# Transform to Kepler orbital elements
def orbkep(orbit):
    return(OrbitType.KEPLERIAN.convertType(orbit))

# Transform to Cartesian orbit
def cartorb_conv(orbit):
    return(OrbitType.CARTESIAN.convertType(orbit))

def period(orbit):
    return(orbkep(orbit).getKeplerianPeriod())
