""""
Representation of state as Cartesian, Kepleran, Circular, Equinoctial

Exported: PVT, posmag, new_cart, elementval, kepler, new_kepler
Examples: ex1, ex2
"""

from force import *
from frames import *
import orbit
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import CircularOrbit
from org.orekit.orbits import EquinoctialOrbit
# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
import collections.abc
from astro import *
#from astropy.table.row import *
# ECIS - exploratory computation in stages
import ecis
from cartprodparam import ensurelist
ecis.ecisdefault='forceenv' # Forces and other environmental constants

# Orekit configuration
_okc = {'cartesian': OrbitType.CARTESIAN}

######## Snaglab (Python) and Orekit representation of a postion-velocity-time (PVT)

# PVT as
#  .ork:  TimeStampedPVCoordinates (Orekit)
#  .pvtq: TQuantity (AstroPy)
class PVT:
    """A position-value-time, represented in two ways,
       pvtq as a TQuantity (based on AstroPy's Quantity)
       ork  as an Orekit TimeStampedPVCoordinates

    They can be made from a
    * list or numpy ndarray of 6 components (position and velocity), and an absolute time specification
    * TQuantity
    * TimeStampedPVCoordinates
    * astropy.table.row.Row
    * tuple in order (astropy.time.Time, Quantity (position), Quantity (velocity))
    * Orbit

    If the units are not in the input, they may be specified in the
    units argument, which defaults to [prefunits["length"],
    prefunits["velocity"]]
    """
    pvtq: TQuantity
    ork: TimeStampedPVCoordinates

    def __init__(self, fromthing, time=None, ork=None,
                 units=[prefunits["length"], prefunits["velocity"]]):
        if listnpa(fromthing):
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
    def scale(self, pvscale):
        # Multiple the position by a scalar (pvscale[0]) and velocity by another scalar (pvscale[1])
        return PVT(scale_posvel(self.pvtq, pvscale))
    def convert_units(self, length_unit=prefunits["length"], velocity_unit=prefunits["velocity"]):
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
    return Quantity(orkval, orkunit).to(prefunits[lookup["phystype"]])

def _elmake(el, value):
    return Quantity(value, prefunits[_eldict[el]["phystype"]])

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
                          to_okad(epoch),   # Sets the date of the orbital parameters
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
# ex2 = new_kepler({'sma_m': 8.0e6, 'ecc': 0.1, 'inc_deg':42.0,
#                   'raan_deg':217.4, 'argper_deg':-90.0,
#                   'timeelt_deg':7.25, 'mean_timeelt':True},
#                  Time('2023-09-14T08:30:00'),
#                  setgravity(0, 0))

ex2 = new_kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "raan":217.4, "ma":7.25},
                 Time('2023-09-14T08:30:00'),
                 setgravity(0, 0))
