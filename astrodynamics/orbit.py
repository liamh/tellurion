""""
Representation of state as Cartesian, Kepleran, Circular, Equinoctial

Exported: PVT, posmag, new_cart, elementval, kepler, new_kepler
Examples: ex1, ex2
"""

import astropy
import astropy.time # For ex1, ex2
import astropy.units as u
import orekit.pyhelpers as pyhelp
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
from org.orekit.utils import Constants
# Orbital elements and PVT
from org.orekit.utils import TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
import collections.abc
#from astropy.table.row import *

from . import astro  # AstroPy and TQuantity
from . import posvel # Vectors and PVT (position, velocity, time) sets in Orekit and AstroPy
from . import dttm   # Dates and times and conversions in various packages: Orekit, AstroPy, NumPy, Python
from . import force
from . import ecis
from . import orbit
from . import util

# ECIS - exploratory computation in stages
ecis.ecisdefault='forceenv' # Forces and other environmental constants

# Orekit configuration
_okc = {'cartesian': OrbitType.CARTESIAN}

######## Snaglab (Python) and Orekit representation of a postion-velocity-time (PVT)

# PVT as
#  .ork:  TimeStampedPVCoordinates (Orekit)
#  .pvtq: astro.TQuantity (AstroPy)
class PVT:
    """A position-value-time, represented in two ways,
       pvtq as a astro.TQuantity (based on AstroPy's u.Quantity)
       ork  as an Orekit TimeStampedPVCoordinates

    They can be made from a
    * list or numpy ndarray of 6 components (position and velocity), and an absolute time specification
    * astro.TQuantity
    * TimeStampedPVCoordinates
    * astropy.table.row.Row
    * tuple in order (astropy.time.Time, u.Quantity (position), u.Quantity (velocity))
    * Orbit

    If the units are not in the input, they may be specified in the
    units argument, which defaults to `[prefunits["length"],
    prefunits["velocity"]]`
    """
    pvtq: astro.TQuantity
    ork: TimeStampedPVCoordinates

    def __init__(self, fromthing, time=None, ork=None,
                 units=[astro.prefunits["length"], astro.prefunits["velocity"]]):
        if util.listnpa(fromthing):
            self.pvtq = posvel.posvel([fromthing[0:3], fromthing[3:6]], time,
                               length_unit=units[0], velocity_unit=units[1])
            if ork is None:
                psi = self.pvtq.si['p'].value.tolist()
                vsi = self.pvtq.si['v'].value.tolist()
                self.ork = posvel.orkpvt(psi,vsi,time)
            else:
                self.ork = ork
        elif type(fromthing) == TimeStampedPVCoordinates:
            self.__init__([fromthing.position.x, fromthing.position.y, fromthing.position.z,
                           fromthing.velocity.x, fromthing.velocity.y, fromthing.velocity.z],
                          fromthing.date.apt(),
                          fromthing,
                          [u.meter, u.meter/u.second])
            self.convert_units(units[0],units[1])
        elif type(fromthing) == astro.TQuantity:
            self.pvtq = fromthing
            psi = fromthing.si['p'].value.tolist()
            vsi = fromthing.si['v'].value.tolist()
            self.ork = posvel.orkpvt(psi, vsi, fromthing.time.datetime)
        elif type(fromthing) == astropy.table.row.Row:
            self.__init__(posvel.posvel([fromthing['position'], fromthing['velocity']], fromthing['time']))
        elif type(fromthing) == tuple and len(fromthing) == 3 and type(fromthing[0]) == astropy.time.Time \
             and type(fromthing[1]) == astro.TQuantity and type(fromthing[2]) == astro.TQuantity:
            self.__init__(astro.posvel([fromthing[1], fromthing[2]], fromthing[0]))
        elif type(fromthing) == Orbit: # The inverse of .cartesian()
            self.__init__(fromthing.pVCoordinates)
    def __repr__(self):
        return f"<PVT position: {self.pvtq['p'].value.tolist()} ({self.pvtq.unit[0].to_string()}) velocity:{self.pvtq['v'].value.tolist()} ({self.pvtq.unit[1].to_string()}) epoch {self.pvtq.time} (UTC)>"
    def scale(self, pvscale):
        # Multiple the position by a scalar (pvscale[0]) and velocity by another scalar (pvscale[1])
        return PVT(posvel._scale_posvel(self.pvtq, pvscale))
    def convert_units(self, length_unit=astro.prefunits["length"], velocity_unit=astro.prefunits["velocity"]):
        self.pvtq = self.pvtq.convert_units((length_unit,velocity_unit))
        return self
    def makenp(self):
        # Make a NumPy object
        return np.concatenate((self.pvtq.value[0], self.pvtq.value[1])), self.pvtq.time.datetime64
    def lla(self, forceenv):
        # Convert to geographic coordinates
        # ex1.pvt.lla(ex1.forceenv)
        return LLA(llafrompt(self.ork.getPosition(), self.ork.getDate(), forceenv),
                   self.pvtq.time)

######## Properties of orbits

def posmag(orbit):
    """The geocentric distance of the orbit"""
    return(orbit.pVCoordinates.position.norm)

######## Make orbits

# Make a Cartesian orbit from Orekit PVT
TimeStampedPVCoordinates.cartesian = lambda self, gravity: CartesianOrbit(self, gravity['celestialframe'], gravity['earthmu'])

def new_cart(pv, time, constants):
    """
    Create a new computation tree based on the Cartesian
    position-velocity and time, and constants (output from
    `setgravity()`) to be used in the computation.
    """
    pvt = PVT(pv, time)
    ret = ecis.newtree('pvt', pvt, constants) # Create the tree and set the first component to the PVT
    ret.cartesian() # Convert the PVT to the Orekit CartesianOrbit and save that as the next component
    ret.kepler()
    return(ret)

################################################################################
## Conversion
################################################################################

# Convert to the requested orbit type "cart", "kep", "circ", "equi"
def _convert(tree, orbtype, searchtype = Orbit):
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

ecis.Ecis.cartesian = lambda tree: _convert(tree, "cart")
ecis.Ecis.kepler = lambda tree: _convert(tree, "kep")
ecis.Ecis.circular = lambda tree: _convert(tree, "circ")
ecis.Ecis.equinoctial = lambda tree: _convert(tree, "equi")

################################################################################
## Orbital elements
################################################################################

def elementval (orbit, element, earthrad=None):
    """
    Compute the orbital element from the orbit
    Arguments
      orbit:    orbit (of any type), may be a list
      element:  the orbital element desired, see list eldict.keys(); may be a list, e.g. ["sma", "ecc"]
      earthrad: the radius of the earth, necessary to provide for altitudes of perigee and apogee
    """
    if isinstance(orbit, collections.abc.Iterable):
        return [elementval(orb, element) for orb in orbit]
    else:
        if type(orbit) is KeplerianOrbit:
            orbkep = orbit
        else:
            orbkep = _convert(orbit,"kep")
        if isinstance(element, list):
            return [_elget(orbkep, el, earthrad) for el in element]
        else:
            return _elget(orbkep, element, earthrad)

_elkeys = ["name", "description", "phystype", "orkunit", "getter"]
_elvals = [["sma", "semimajor axis", "length", u.meter, KeplerianOrbit.getA],
          ["ecc", "eccentricity", "dimensionless", u.dimensionless_unscaled, KeplerianOrbit.getE],
          ["inc", "inclination", "angle", u.radian, KeplerianOrbit.getI],
          ["argper", "argument of perigee", "angle", u.radian, KeplerianOrbit.getPerigeeArgument],
          ["radper", "radius of perigee", "length", u.meter,
           lambda kep: kep.a*(1.0-kep.e)],
          ["radapo", "radius of apogee", "length", u.meter,
           lambda kep: kep.a*(1.0+kep.e)],
          ["altper", "altitude of perigee", "length", u.meter,
           lambda kep, earthrad: kep.a*(1.0-kep.e)-earthrad],
          ["altapo", "altitude of apogee", "length", u.meter,
           lambda kep, earthrad: kep.a*(1.0+kep.e)-earthrad],
          ["raan", "right ascension of the ascending node", "angle", u.radian,
           KeplerianOrbit.getRightAscensionOfAscendingNode],
          ["ta", "true anomaly", "angle", u.radian, KeplerianOrbit.getTrueAnomaly],
          ["ma", "mean anomaly", "angle", u.radian, KeplerianOrbit.getMeanAnomaly],
          ["memo", "mean motion", "angular speed", u.radian/u.second, KeplerianOrbit.getKeplerianMeanMotion],
          ["period", "orbital period", "time", u.second, KeplerianOrbit.getKeplerianPeriod]]
_eldict = dict(zip([ev[0] for ev in _elvals], [dict(zip(_elkeys,ev)) for ev in _elvals]))

# Get altitude of perigee/apogee by subtracting ex2.prop.forceenv["earthrad"]
def _elget(orbkep, el, earthrad=None):
    lookup = _eldict[el]
    getter = lookup["getter"]
    if '__code__' in dir(getter) and len(getter.__code__.co_varnames) > 1:
        orkval = getter(orbkep, earthrad)
    else:
        orkval = getter(orbkep)
    orkunit = lookup["orkunit"]
    return u.Quantity(orkval, orkunit).to(astro.prefunits[lookup["phystype"]])

def _elmake(el, value):
    return u.Quantity(value, astro.prefunits[_eldict[el]["phystype"]])

def _orkkep(el, oes):
    return float(oes[el].to(_eldict[el]['orkunit']).value)

def _makekep(elvald):
    return {e: _elmake(e, elvald[e]) for e in elvald}

def kepler(oes, epoch, constants):
    foes = _makekep({"ecc":0.0, "inc":38.0, "raan":0.0, "argper":-90.0, "ma":0}) | _makekep(oes)
    if 'ma' in foes:
        timeelt_type = PositionAngleType.MEAN
        timeelt = _orkkep('ma',foes)
    else:
        timeelt_type = PositionAngleType.TRUE
        timeelt = _orkkep('ta',foes)
    return KeplerianOrbit(_orkkep('sma',foes),
                          _orkkep('ecc',foes),
                          _orkkep('inc',foes),
                          _orkkep('argper',foes),
                          _orkkep('raan',foes),
                          timeelt,
                          timeelt_type,  # Sets which type of anomaly we use (true
                          constants['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                          dttm.to_okad(epoch),   # Sets the date of the orbital parameters
                          constants['earthmu'])   # Sets the central attraction coefficient (m³/s²)

def new_kepler(oes, epoch, constants):
    """
    Create a new computation tree based on the Keplerian
    (classical) orbital elements and time, and constants (output from
    `setgravity()`) to be used in the computation.
    """
    ret = ecis.newtree('kep', kepler(oes, epoch, constants), constants)
    ret.cartesian()
    return(ret)
