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
import collections.abc
# Representation using AstroPy
from astro import *
# ECIS - exploratory computation in stages
import ecis
from cartprodparam import ensurelist
ecis.ecisdefault='forceenv' # Forces and other environmental constants

# Orekit configuration
okc = {'cartesian': OrbitType.CARTESIAN}

######## Snaglab (Python) and Orekit representation of a postion-velocity-time (PVT)

# PVT as
#  .ork:  TimeStampedPVCoordinates (Orekit)
#  .pvtq: TQuantity (AstroPy)
class PVT:
    pvtq: TQuantity
    ork: TimeStampedPVCoordinates

    def __init__(self, fromthing, time=None, ork=None,
                 units=[prefunits["length"], prefunits["velocity"]]):
        if type(fromthing)==list or type(fromthing) == np.ndarray:
            self.pvtq = posvel([fromthing[0:3], fromthing[3:6]], time,
                               length_unit=units[0], velocity_unit=units[1])
            if ork is None:
                psi = self.pvtq.si['p'].value.tolist()
                vsi = self.pvtq.si['v'].value.tolist()
                self.ork = TimeStampedPVCoordinates(to_okad(time), Vector3D(psi), Vector3D(vsi))
            else:
                self.ork = ork
        elif type(fromthing) == TimeStampedPVCoordinates:
            self.__init__([fromthing.position.x, fromthing.position.y, fromthing.position.z,
                           fromthing.velocity.x, fromthing.velocity.y, fromthing.velocity.z],
                          fromthing.date.apt(),
                          fromthing,
                          [u.meter, u.meter/u.second])
            self.convert_units(units[0],units[1])
        elif type(fromthing) == TQuantity:
            self.pvtq = fromthing
            psi = fromthing.si['p'].value.tolist()
            vsi = fromthing.si['v'].value.tolist()
            self.ork = TimeStampedPVCoordinates(datetime_to_absolutedate(fromthing.time.datetime),
                                                Vector3D(psi), Vector3D(vsi))
        elif type(fromthing) == astropy.table.row.Row:
            self.__init__(posvel([fromthing['position'], fromthing['velocity']], fromthing['time']))
        elif type(fromthing) == tuple and len(fromthing) == 3 and type(fromthing[0]) == Time \
             and type(fromthing[1]) == TQuantity and type(fromthing[2]) == TQuantity:
            self.__init__(posvel([fromthing[1], fromthing[2]], fromthing[0]))
        elif type(fromthing) == Orbit: # The inverse of .cartesian()
            self.__init__(fromthing.pVCoordinates)
    def __repr__(self):
        return f"<PVT position: {self.pvtq['p'].value.tolist()} ({self.pvtq.unit[0].to_string()}) velocity:{self.pvtq['v'].value.tolist()} ({self.pvtq.unit[1].to_string()}) epoch {self.pvtq.time} (UTC)>"
    def convert_units(self, length_unit=prefunits["length"], velocity_unit=prefunits["velocity"]):
        self.pvtq = self.pvtq.convert_units((length_unit,velocity_unit))
        return self
    def makenp(self):
        return (np.concatenate((self.pvtq.value[0], self.pvtq.value[1])), self.pvtq.time.datetime64)

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
                          to_okad(epoch),   # Sets the date of the orbital parameters
                          constants['earthmu']))   # Sets the central attraction coefficient (m³/s²)


######## Properties of orbits

# The geocentric distance of the orbit
def posmag(orbit):
    return(orbit.pVCoordinates.position.norm)

def period(orbit):
    return(orbit.getKeplerianPeriod())

## Make a computation tree from position, velocity, and datetime

def new_cart(pv, time, constants):
    pvt = PVT(pv, time)
    ret = ecis.newtree('pvt', pvt, constants) # Create the tree and set the first component to the PVT
    ret.cartesian() # Convert the PVT to the Orekit CartesianOrbit and save that as the next component
    ret.kepler()
    return(ret)

## Make a computation tree from Kepler elements and datetime

def new_kepler(oes, epoch, constants):
    ret = ecis.newtree('kep', kepler(oes, epoch, constants), constants)
    ret.cartesian()
    return(ret)

################################################################################
## Conversion
################################################################################

# Convert to the requested orbit type "cart", "kep", "circ", "equi"
def convert(tree, orbtype, searchtype = Orbit):
    [orbit, parent] = ecis.thingofclass(tree, searchtype)
    if orbit is None:
        [orbit, parent] = ecis.thingofclass(tree, PVT)
        parent.update(cart = parent.pvt.ork.cartesian(tree.forceenv))
        orbit = parent.cart
    match orbtype:
        case "cart":
            ret = OrbitType.CARTESIAN.convertType(orbit)
            if parent is not None:
                parent.update(cart = ret)
        case "kep":
            # Must cast; see https://forum.orekit.org/t/convert-orbit-to-keplerian/1441/2
            ret = KeplerianOrbit.cast_(OrbitType.KEPLERIAN.convertType(orbit))
            if parent is not None:
                parent.update(kep = ret)
        case "circ":
            ret = CircularOrbit.cast_(OrbitType.CIRCULAR.convertType(orbit))
            if parent is not None:
                parent.update(circ = ret)
        case "equi":
            ret = EquinoctialOrbit.cast_(OrbitType.EQUINOCTIAL.convertType(orbit))
            if parent is not None:
                parent.update(equi = ret)
        case _:
            raise ValueError("Type \"" + orbtype + "\" unknown")
    return(ret)

ecis.Ecis.cartesian = lambda tree: convert(tree, "cart")
ecis.Ecis.kepler = lambda tree: convert(tree, "kep")
ecis.Ecis.circular = lambda tree: convert(tree, "circ")
ecis.Ecis.equinoctial = lambda tree: convert(tree, "equi")

################################################################################
## Get orbital elements
################################################################################

# Both `orbit` and `element` can be a lists e.g. ["sma", "ecc"]
# orbit:   orbit (of any type)
# element: the orbital element desired, see list in elvals
def elementval (orbit, element):
    if isinstance(orbit, collections.abc.Iterable):
        return [elementval(orb, element) for orb in orbit]
    else:
        if type(orbit) is KeplerianOrbit:
            orbkep = orbit
        else:
            orbkep = convert(orbit,"kep")
        if isinstance(element, list):
            return [elget(orbkep, el) for el in element]
        else:
            return elget(orbkep, element)

elkeys = ["name", "description", "phystype", "orkunit", "getter"]
elvals = [["sma", "semimajor axis", "length", u.meter, KeplerianOrbit.getA],
          ["ecc", "eccentricity", "dimensionless", None, KeplerianOrbit.getE],
          ["inc", "inclination", "angle", u.radian, KeplerianOrbit.getI],
          ["argper", "argument of perigee", "angle", u.radian, KeplerianOrbit.getPerigeeArgument],
          ["raan", "right ascension of the ascending node", "angle", u.radian,
           KeplerianOrbit.getRightAscensionOfAscendingNode],
          ["ta", "true anomaly", "angle", u.radian, KeplerianOrbit.getTrueAnomaly],
          ["ma", "mean anomaly", "angle", u.radian, KeplerianOrbit.getMeanAnomaly],
          ["memo", "mean motion", "angular speed", u.radian/u.second, KeplerianOrbit.getKeplerianMeanMotion],
          ["period", "orbital period", "time", u.second, KeplerianOrbit.getKeplerianPeriod]]
eldict = dict(zip([ev[0] for ev in elvals], [dict(zip(elkeys,ev)) for ev in elvals]))

def elget(orbkep, el):
    lookup = eldict[el]
    orkval = lookup["getter"](orbkep)
    orkunit = lookup["orkunit"]
    if orkunit is None:
        return Quantity(orkval)
    else:
        return Quantity(orkval, orkunit).to(prefunits[lookup["phystype"]])


################################################################################
## Examples
################################################################################

ex1 = new_cart([5740.13268349499, 3314.06715, 0.0,
                -2.75082683526322, 4.7645718414998, 5.50165367052644],
               Time('2022-06-01T12:00:00.000000'),
               setgravity(0,0))
# All these are the same as ex1.pvt: PVT(ex1.cart), PVT(ex1.pvt.pvtq), PVT(ex1.pvt.ork), PVT(*ex1.pvt.makenp())
# In [3]: ex1.pvt
# Out[3]: <PVT position: [5740.13268349499, 3314.06715, 0.0] (km) velocity:[-2.75082683526322, 4.7645718414998, 5.50165367052644] (km/s) epoch 2022-06-01T12:00:00.000 (UTC)>
# In [4]: ex1.pvt.pvtq
# Out[4]: <TQuantity ([5740.13268349, 3314.06715   ,    0.        ], [-2.75082684,  4.76457184,  5.50165367]) (km, km / s), time=2022-06-01T12:00:00.000>
# In [8]: ex1.pvt.ork
# Out[8]: <TimeStampedPVCoordinates: {2022-06-01T12:00:00.000, P(5740132.68349499, 3314067.15, 0.0), V(-2750.82683526322, 4764.5718414998, 5501.65367052644), A(0.0, 0.0, 0.0)}>
# elementval(ex1.cart, "sma")
# Out[20]: <Quantity 6672.37441023 km>

# ex1.cartesian()

# Build and convert a Kepler
ex2 = new_kepler({'sma_m': 8.0e6, 'ecc': 0.1, 'inc_deg':42.0,
                  'raan_deg':217.4, 'argper_deg':-90.0,
                  'timeelt_deg':7.25, 'mean_timeelt':True},
                 Time('2023-09-14T08:30:00'),
                 setgravity(0, 0))
