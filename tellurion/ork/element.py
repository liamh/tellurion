"""
Orbital elements in Orekit
"""

import collections.abc
import astropy.units as u
from astropy.timeseries import TimeSeries
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType, CircularOrbit, KeplerianOrbit, PositionAngleType
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates

from ..core import astro
from ..core import posvel
from ..core import element
from . import force
from . import convert

########################################
####    Element values              ####
########################################

def tselements(ephem, elements):
    '''Make a time series of selected orbital elements'''
    return TimeSeries(time=ephem.time,
                      data=[dict(zip(elements, elementval(ephrow, elements)))
                            for ephrow in ephem])

def elementval (orbstate, elt, earthrad=force.deffe["earthrad"]):
    """
    Compute the orbital element from the orbital state
    Arguments
      orbstate:    Representation of orbital state in any form
      elt:  the orbital element desired, see list eldict.keys(); may be a list, e.g. ["sma", "ecc"]
      earthrad: the radius of the earth, necessary to provide for altitudes of perigee and apogee
    """
    if isinstance(elt, list):
        return [elementval(orbstate, el, earthrad) for el in elt]
    else:
        lookup = _eldict[elt]
        getter = lookup["getter"]
        orbkep = _keplerianorbit(orbstate)
        if '__code__' in dir(getter) and len(getter.__code__.co_varnames) > 1:
            orkval = getter(orbkep, earthrad)
        else:
            orkval = getter(orbkep)
        orkunit = lookup["orkunit"]
        return u.Quantity(orkval, orkunit).to(astro.prefunits[lookup["phystype"]])

_elkeys = ["name", "description", "phystype", "orkunit", "getter"]
_elvals = [["sma", "semimajor axis", "length", u.meter, KeplerianOrbit.getA],
          ["ecc", "eccentricity", "dimensionless", u.dimensionless_unscaled, KeplerianOrbit.getE],
          ["inc", "inclination", "angle", u.radian, KeplerianOrbit.getI],
          ["argper", "argument of perigee", "angle", u.radian, KeplerianOrbit.getPerigeeArgument],
          ["radper", "radius of perigee", "length", u.meter,
           lambda kep: kep.getA()*(1.0-kep.getE())],
          ["radapo", "radius of apogee", "length", u.meter,

           lambda kep: kep.getA()*(1.0+kep.getE())],
          ["altper", "altitude of perigee", "length", u.meter,
           lambda kep, earthrad: kep.getA()*(1.0-kep.getE())-earthrad],
          ["altapo", "altitude of apogee", "length", u.meter,
           lambda kep, earthrad: kep.getA()*(1.0+kep.getE())-earthrad],
          ["raan", "right ascension of the ascending node", "angle", u.radian,
           KeplerianOrbit.getRightAscensionOfAscendingNode],
          ["ta", "true anomaly", "angle", u.radian, KeplerianOrbit.getTrueAnomaly],
          ["ma", "mean anomaly", "angle", u.radian, KeplerianOrbit.getMeanAnomaly],
          ["memo", "mean motion", "angular speed", u.radian/u.second, KeplerianOrbit.getKeplerianMeanMotion],
          ["period", "orbital period", "time", u.second, KeplerianOrbit.getKeplerianPeriod]]
_eldict = dict(zip([ev[0] for ev in _elvals], [dict(zip(_elkeys,ev)) for ev in _elvals]))

def _elmake(el, value):
    return u.Quantity(value, astro.prefunits[_eldict[el]["phystype"]])

def _orkkep(el, oes):
    return float(oes[el].to(_eldict[el]['orkunit']).value)

def _makekep(elvald):
    return {e: _elmake(e, elvald[e]) for e in elvald}

########################################
####    Make Kepler element set     ####
########################################

def _keplerianorbit(oes, units=(astro.prefunits['length'], astro.prefunits['angle']),
                    forceenv=force.deffe):
    '''Make a org.orekit.orbits.KeplerianOrbit from anything'''
    if element.iskepels(oes):
        return _keporb_from_components(*oes, units, forceenv)
    elif type(oes) is KeplerianOrbit:
        return oes
    elif posvel.ispvt(oes):
        co = CartesianOrbit(convert._tspvc(*oes),
                            forceenv['celestialframe'], forceenv['earthmu'])
    elif posvel.isephrow(oes):
        co = CartesianOrbit(convert._tspvc(oes),
                            forceenv['celestialframe'], forceenv['earthmu'])
    else:
        raise ValueError('Cannot transform to Keplerian elements')
    return OrbitType.KEPLERIAN.convertType(co)

def _keporb_from_components(oes, epoch, units=(astro.prefunits['length'], astro.prefunits['angle']), fe=force.deffe):
    '''Make a org.orekit.orbits.KeplerianOrbit from orbital elements as a u.Quantity or Dict'''
    oessi = oes.si.value
    # Future: convert zp/za to a/e
    #    oessidict = astro.splitsq(oessi)

    if 'ma' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.MEAN, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              convert._okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'])   # Sets the central attraction coefficient (m³/s²)
    elif 'ta' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.TRUE, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              convert._okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'])   # Sets the central attraction coefficient (m³/s²)
    else:
        raise ValueError('Time element (ma or ta) required in element set')

###############################
####  Transformations      ####
###############################

# Transformations between Cartesian state vector and Kepler elements

def kepler(object, forceenv=force.deffe, mean_time_element=True): # Add prefunits
    '''The Kepler element set from the Cartesian PVT or equivalent'''
    co = CartesianOrbit(convert._tspvc(*object),
                        forceenv['celestialframe'], forceenv['earthmu'])
    ko = OrbitType.KEPLERIAN.convertType(co)
    if mean_time_element:
        elnames = element.kepeltma_names
    else:
        elnames = element.kepeltta_names
    kepels = astro.makesq(elementval(ko, elnames), elnames)
    return (kepels, object[1])

def cartesian(object, dttm=None):
    '''Make a Cartesian PVT from the object, or an ephemeris
    generator, which has an initial state. If the object is an
    elements set without a datetime, it must be supplied in `dttm`.'''
    if element.iskepels(object):
        if type(object) is tuple:
            return convert._pvt(_keplerianorbit(object))
        else:
            return convert._pvt(_keplerianorbit(object, dttm))
    else:
        raise ValueError('Can only transform Kepler element sets')
