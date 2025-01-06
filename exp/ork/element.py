""""
Orbital elements in Orekit
"""

import collections.abc
import astropy.units as u
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.orbits import KeplerianOrbit, PositionAngleType

import astro
import cdttm
import posvel
import ork.force as ofr
import ork.posvel as opv

def elementval (orbit, elt, earthrad=None):
    """
    Compute the orbital element from the orbit
    Arguments
      orbit:    orbit (of any type), may be a list
      elt:  the orbital element desired, see list eldict.keys(); may be a list, e.g. ["sma", "ecc"]
      earthrad: the radius of the earth, necessary to provide for altitudes of perigee and apogee
    """
    if isinstance(orbit, collections.abc.Iterable):
        return [elementval(orb, elt) for orb in orbit]
    else:
        if type(orbit) is KeplerianOrbit:
            orbkep = orbit
        else:
            orbkep = _convert(orbit,"kep")
        if isinstance(elt, list):
            return [_elget(orbkep, el, earthrad) for el in elt]
        else:
            return _elget(orbkep, elt, earthrad)

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

def kepler(oes, epoch=None, constants=ofr.deffe):
    '''Find the org.orekit.orbits.KeplerianOrbit from the input'''
    if posvel.ispvt(oes) and epoch is None:
        return opv.orkpvt(*oes).keplerianorbit()
    if posvel.ispv(oes):
        return opv.orkpvt(oes, epoch).keplerianorbit()
    if posvel.isephrow(oes):
        return kepler(posvel.pvt(oes))
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
                          cdttm.okad(epoch),   # Sets the date of the orbital parameters
                          constants['earthmu'])   # Sets the central attraction coefficient (m³/s²)

# Convert from KeplerianOrbit
KeplerianOrbit.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
KeplerianOrbit.pvt = lambda self: opv.pvtork(self.pVCoordinates)

# Convert to KeplerianOrbit
TimeStampedPVCoordinates.keplerianorbit = lambda self, gravity=ofr.deffe: self.cartesianorbit(gravity).keplerianorbit()
CartesianOrbit.keplerianorbit = lambda self: KeplerianOrbit.cast_(OrbitType.KEPLERIAN.convertType(self))
Orbit.keplerianorbit = lambda self: KeplerianOrbit.cast_(OrbitType.KEPLERIAN.convertType(self))
