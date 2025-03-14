"""
Orbital elements in Orekit
"""

import collections.abc
import astropy.units as u
from astropy.timeseries import TimeSeries
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType, CircularOrbit
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.orekit.orbits import KeplerianOrbit, PositionAngleType

from ..core import astro
from ..core import posvel
from ..core import element
from . import force
from . import posvel as oposvel

########################################
####  Convert PVT to Kepler elset   ####
########################################

def kepler(pvt, units=(astro.prefunits['length'], astro.prefunits['angle'])):
    '''Convert Cartesian PVT to a Kepler orbital element set'''
    return oposvel.tspvc(*pvt).kepler(units)

########################################
####    Element values              ####
########################################

def tselements(ephem, elements):
    '''Make a time series of selected orbital elements'''
    return TimeSeries(time=ephem.time,
                      data=[dict(zip(elements, elementval(ephrow, elements)))
                            for ephrow in ephem])

def elementval (orbit, elt, earthrad=force.deffe["earthrad"]):
    """
    Compute the orbital element from the orbit
    Arguments
      orbit:    orbit (of any type), may be a list
      elt:  the orbital element desired, see list eldict.keys(); may be a list, e.g. ["sma", "ecc"]
      earthrad: the radius of the earth, necessary to provide for altitudes of perigee and apogee
    """
    if posvel.isephrow(orbit):
        return elementval(oposvel.tspvc(orbit), elt, earthrad)
    elif isinstance(orbit, collections.abc.Iterable):
        return [elementval(orb, elt) for orb in orbit]
    else:
        orbkep = orbit.keplerianorbit()
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

########################################
####    Make Kepler element set     ####
########################################

def keplerianorbit(oes, epoch, units=(astro.prefunits['length'], astro.prefunits['angle']), fe=force.deffe):
    '''Make a org.orekit.orbits.KeplerianOrbit from orbital elements as a u.Quantity or Dict'''
    if type(oes) is dict:
        oes = element.kepler(oes, None, units)

    oessi = oes.si.value
# Future: convert zp/za to a/e
#    oessidict = astro.splitsq(oessi)

    if 'ma' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.MEAN, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              oposvel.okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'])   # Sets the central attraction coefficient (m³/s²)
    elif 'ta' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.TRUE, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              oposvel.okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'])   # Sets the central attraction coefficient (m³/s²)
    else:
        raise ValueError('Time element (ma or ta) required in element set')
    return

def _kepler(orkobj, units=(astro.prefunits['length'], astro.prefunits['angle'])):
    '''Make a kepler u.Quantity from the Orekit object.'''
    return element.kepler({'sma': elementval(orkobj, 'sma'),
                           'ecc': elementval(orkobj, 'ecc'),
                           'inc': elementval(orkobj, 'inc'),
                           'argper': elementval(orkobj, 'argper'),
                           'raan': elementval(orkobj, 'raan'),
                           'ma': elementval(orkobj, 'ma')},
                          oposvel.okad(orkobj.getDate()),
                          units)

########################################
####    Convert element types       ####
########################################


# Convert to KeplerianOrbit (Orekit)
#TimeStampedPVCoordinates.keplerianorbit = lambda self, gravity=force.deffe: self.cartesianorbit(gravity).keplerianorbit()
# Orbit.keplerianorbit = lambda self: KeplerianOrbit.cast_(OrbitType.KEPLERIAN.convertType(self))

# Convert to kepler (u.Quantity)
#TimeStampedPVCoordinates.kepler = \
#    lambda self, units=(astro.prefunits['length'], astro.prefunits['angle']): \
#        self.cartesianorbit(force.deffe).kepler(units)
#Orbit.kepler = lambda self, units=(astro.prefunits['length'], astro.prefunits['angle']): _kepler(self, units)

# Circular orbit
# TimeStampedPVCoordinates.circularorbit = \
#    lambda self, gravity=force.deffe: CircularOrbit(self, gravity['celestialframe'], gravity['earthmu'])
#Orbit.circularorbit = lambda self: CircularOrbit(self)
