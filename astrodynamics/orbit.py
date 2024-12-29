""""
Representation of state as Cartesian, Kepleran, Circular, Equinoctial

Exported: PVT, posmag, new_cart, elementval, kepler, new_kepler
Examples: ex1, ex2
"""

import astropy
import astropy.time # For ex1, ex2
import astropy.units as u
import warnings
import numpy as np
import orekit.pyhelpers as pyhelp
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
from org.orekit.utils import Constants
# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
import collections.abc
#from astropy.table.row import *

from . import astro  # AstroPy and TQuantity
from . import pvatq # Vectors and PVT (position, velocity, time) sets in AstroPy
from . import pvork  # Vectors and PVT (position, velocity, time) sets in Orekit
from . import dttm   # Dates and times and conversions in various packages: Orekit, AstroPy, NumPy, Python
from . import force
from . import frames
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
#  .atq: astro.TQuantity (AstroPy)
class PVT:
    """A position-value-time, represented in two ways,
       atq as a astro.TQuantity (based on AstroPy's u.Quantity)
       ork as an Orekit TimeStampedPVCoordinates

    They can be made from a
    * list or numpy ndarray of 6 components (position and velocity), and an absolute time specification
    * astro.TQuantity
    * TimeStampedPVCoordinates
    * astropy.table.row.Row
    * tuple in order (astropy.time.Time, u.Quantity (position), u.Quantity (velocity))
    * Orbit

    If the units are not in the input, they may be specified as a tuple in the
    units argument, which defaults to `(prefunits["length"],
    prefunits["velocity"])`
    """
    atq: astro.TQuantity
    ork: TimeStampedPVCoordinates

    def __init__(self, fromthing, time=None, ork=None, units=astro.prefunits["posvel"]):

        if (type(fromthing) is TimeStampedPVCoordinates or type(fromthing) is PVCoordinates) \
           or type(fromthing)==Orbit:
            # fromthing is an Orekit object
            if type(fromthing)==Orbit:
                ft=fromthing.pVCoordinates
            else:
                ft=fromthing
            (pos, vel, time) = pvork.pvtork(ft)
            atq = pvatq.atqptpvt(pos, vel, time, pvatq.posvelsiu).convert_units(units)
            self.ork = fromthing
        else: # fromthing is not an Orekit object
            atq = pvatq.atqpvt(fromthing, time, units=units)
            self.ork = pvork.orkpvt(atq.si['p'].value.tolist(), \
                                    atq.si['v'].value.tolist(), time)
        self.atq = atq

    def __repr__(self):
        if 'v' in self.atq.dtype.names:
            return f"<PVT position: {self.atq['p'].value.tolist()} ({self.atq.unit[0].to_string()}) " \
                f"velocity:{self.atq['v'].value.tolist()} ({self.atq.unit[1].to_string()}) " \
                f"epoch {self.atq.time} (UTC)>"
        else:
            return f"<PT position: {self.atq['p'].value.tolist()} ({self.atq.unit[0].to_string()}) " \
                f"epoch {self.atq.time} (UTC)>"
    def scale(self, pvscale):
        # Multiple the position by a scalar (pvscale[0]) and velocity by another scalar (pvscale[1])
        return PVT(pvatq._scale_posvel(self.atq, pvscale))
    def convert_units(self, units=astro.prefunits["posvel"]):
        self.atq = self.atq.convert_units(units)
        return self
    def makenp(self):
        # Make a NumPy object
        return np.concatenate((self.atq.value[0], self.atq.value[1])), self.atq.time.datetime64
    def lla(self, forceenv=force.setgravity(0,0)):
        '''Convert to geographic coordinates'''
        return frames.llafrompt(self.ork.getPosition(), self.ork.getDate(), forceenv)
        #return frames.LLA(forceenv['earth'].transform(self.ork.getPosition(), \
        #                                              forceenv['celestialframe'], self.ork.getDate()))

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
    # ret.kepler() # This gets preferred for propagation, don't want that
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
