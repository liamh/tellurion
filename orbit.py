######### Representation of state as Cartesian, Kepleran, Circular, Equinoctial

## A position-velocity-time PVT is a Cartesian state, but
## Orekit has a separate class of CartesianOrbit that can be
## propagated, so these are have different representation.

from force import *
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
# ECIS - exploratory computation in stages
from ecis import *

# Orekit configuration
okc = {'cartesian': OrbitType.CARTESIAN}

######## Snaglab (Python) and Orekit representation of a postion-velocity-time (PVT)

# PVT as Orekit arrays and AbsoluteDate
# These are made with the `ork` variable set to a TimeStampedPVCoordinates instance
class PVT:
    position: list # Position 3-vector in meters
    velocity: list # Velocity 3-vector in meters/seconds
    dttm: Time # Epoch time
    ork: TimeStampedPVCoordinates

    def __init__(self, position, velocity, dttm, ork=None):
        self.position = position
        self.velocity = velocity
        self.dttm = dttm
        if ork is None:
            self.ork = TimeStampedPVCoordinates(self.dttm.okad(),
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
                       apt(self.date),
                       self)

Quantity.snl \
    = lambda self, date: PVT(ex1pvq.si['p'][0].value.tolist(),
                             ex1pvq.si['v'][0].value.tolist(),
                             date)

######## Make orbits

# Make a Cartesian orbit from Orekit PVT
TimeStampedPVCoordinates.cartesian = lambda self, gravity: CartesianOrbit(self, gravity['celestialframe'], gravity['earthmu'])

# Make a Kepler orbital element set
# `oes` is the orbital element set as a dictionary with the following elements
# zapo_m and zper_m or sma_m and ecc
# 'inc_deg' (default 0)
# 'argper_deg' (default -90.0)    Perigee argument (deg)
# 'raan_deg' (default 0.0)        Right ascension of ascending node (degrees)
# 'timeelt_deg' (default 0.0)     Time element (deg)
# 'mean_timeelt' True or False (default)  Whether time element is true or mean anomaly
def kepler(oes, epoch, constants):
    zapo = oes.get('zapo_m')
    zper = oes.get('zper_m')
    sma = oes.get('sma_m')
    ecc = oes.get('ecc')
    if zapo is not None and zper is not None:
        rapo = constants['earthrad'] + zapo
        rper = constants['earthrad'] + zper
        sma = (rapo+rper)/2.0
        ecc = (rapo-rper)/(2.0*sma)
    if oes.get('mean_timeelt', False):
        timeelt_type = PositionAngleType.MEAN
    else:
        timeelt_type = PositionAngleType.TRUE
    return(KeplerianOrbit(sma, # Semimajor Axis (m)
                          ecc,    # Eccentricity
                          radians(oes.get('inc_deg', 0.0)),  # Inclination (deg)
                          radians(oes.get('argper_deg', -90.0)),   # Perigee argument (deg)
                          radians(oes.get('raan_deg', 0.0)),   # Right ascension of ascending node (degrees)
                          radians(oes.get('timeelt_deg', 0.0)),  # Time element (deg)
                          timeelt_type,  # Sets which type of anomaly we use (true
                          constants['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                          epoch.okad(),   # Sets the date of the orbital parameters
                          constants['earthmu']))   # Sets the central attraction coefficient (m³/s²)

######## Convert orbits

# Convert to the requested orbit type "cart", "kep"
def convert(tree, orbtype):
    [orbit, parent] = thingofclass(tree, Orbit)
    if orbit is None:
        [orbit, parent] = thingofclass(tree, PVT)
        parent.update(cart = parent.pvt.ork.cartesian(tree.default))
        orbit = parent.cart
    match orbtype:
        case "cart":
            ret = OrbitType.CARTESIAN.convertType(orbit)
            parent.update(cart = ret)
        case "kep":
            ret = OrbitType.KEPLERIAN.convertType(orbit)
            parent.update(kep = ret)
        case "circ":
            ret = OrbitType.CIRCULAR.convertType(orbit)
            parent.update(circ = ret)
        case "equi":
            ret = OrbitType.EQUINOCTIAL.convertType(orbit)
            parent.update(equi = ret)
        case "pvt":
            ret = orbit.posveltime()
            parent.update(pvt = ret)
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

def new_posveltime(pv, dttm, constants):
#    if type(pv) is Quantity
    if type(pv) is list:
        if len(pv) == 6:
            pos = pv[0:3]
            vel = pv[3:6]
        elif len(posvel):
            pos = pv[0]
            vel = pv[1]
    ret = newtree('pvt', PVT(pos, vel, dttm), constants) # Create the tree and set the first component to the PVT
    ret.cartesian() # Convert the PVT to the Orekit CartesianOrbit and save that as the next component
    return(ret)

## Make a computation tree from Kepler elements and datetime

def new_kepler(oes, epoch, constants):
    ret = newtree('kep', kepler(oes, epoch, constants), constants)
    convert(ret, "pvt")
    return(ret)

## Example orbit
ex1pvt = posvel([5740.13268349499, 3314.06715, 0.0],
                [-2.75082683526322, 4.7645718414998, 5.50165367052644]).snl(Time('2022-06-01T12:00:00.000000'))
ex1 = new_posveltime([5740132.68349499, 3314067.15, 0.0,
                      -2750.82683526322, 4764.5718414998, 5501.65367052644],
                     Time('2022-06-01T12:00:00.000000'),
                     setgravity(0,0))
ex1.pvtorbrec = ex1.cart.posveltime() # The PVT recalculated from the Cartesian orbit
ex1.pvtorkrec = ex1.pvt.ork.snl() # The PVT recalculated from the Orekit representation
# ex1.keys()
# ex1.cartesian()

# Build and convert a Kepler
ex2 = new_kepler({'sma_m': 8.0e6, 'ecc': 0.1, 'inc_deg':42.0,
                  'raan_deg':217.4, 'argper_deg':-90.0,
                  'timeelt_deg':7.25, 'mean_timeelt':True},
                 Time('2023-09-14T08:30:00'),
                 setgravity(0, 0))
