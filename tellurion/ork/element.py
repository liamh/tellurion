"""
Orbital elements in Orekit
"""

import numpy as np
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

def allplane(oes, forceenv=force.deffe, units=astro.prefunits):
    '''Generate all plane pairs (sma, ecc), (radper, radapo), (altper, altapo) from the first or last pairs; additionally, the mean motion can be substituted for semimajor axis in the first pair.

    Example 1, convert from altitudes of perigee and apogee to semimajor axis and eccentricity
    byalts = tell.kepler({"altper":160*u.km, "altapo":20250*u.km, \
                          "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg}, \
                          tell.dttm('2022-02-15T08:30:00'))
    smaecc = tork.allplane(byalts)
    smaecc[0]['sma'] # <Quantity 16583.13646 km>
    smaecc[0]['ecc'] # <Quantity 0.60573583>
    tell.iskepels(byalts) # False
    tell.iskepels(smaecc) # True

    Example 2, convert from semimajor axis and eccentricity to altitudes of perigee and apogee
    bysmaecc = tell.kepler({"sma":8000, "ecc":0.1, \
                          "inc":45, "argper": 120.0, "raan": 80.0, "ma": 0.0}, \
                          tell.dttm('2025-02-01T12:30:00'))
    alts = tork.allplane(bysmaecc)
    alts[0]['altper'] # <Quantity 821.86354 km>
    alts[0]['altapo'] # <Quantity 2421.86354 km>

    Example 3, define a geosynchronous orbit by mean motion
    geo = tork.allplane(tell.kepler({"memo":1.0*u.rev/u.sday, "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg}))

    Example 4, define a geosynchronous transfer orbit
    gto = tork.allplane(tell.kepler({"altper": 350*u.km, "altapo": tork.smamemo(1.0,True), \
                          "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg}))
    '''
    if type(oes) is tuple:
        (oesd, oest) = oes
    else:
        oesd = oes
        oest = None
    names = oesd.dtype.names
    if ('sma' in names or 'memo' in names) and 'ecc' in names:
        if 'memo' in names:
            sma = smamemo(oesd['memo'], forceenv, units)
            new = {'sma': sma, 'radper': sma*(1-oesd['ecc']), 'radapo': sma*(1+oesd['ecc'])}
        else:
            sma = oesd['sma']
            memo = np.sqrt((forceenv['earthmu']*u.m**3/u.s**2)/sma**3)
            new = {'memo': memo, 'radper': sma*(1-oesd['ecc']), 'radapo': sma*(1+oesd['ecc'])}
        new['altper'] = new['radper'] - forceenv['earthrad']*u.m
        new['altapo'] = new['radapo'] - forceenv['earthrad']*u.m
    elif 'altper' in names and 'altapo' in names:
        new = {'radper': oesd['altper'] + forceenv['earthrad']*u.m, \
               'radapo': oesd['altapo'] + forceenv['earthrad']*u.m}
        new['sma'] = (new['radapo']+new['radper'])/2
        new['ecc'] = (new['radapo']-new['radper'])/(new['radapo']+new['radper'])
    else:
        raise ValueError('Plane must be defined by either (sma, ecc) or (altper, altapo)')
    all = astro.splitsq(oesd) | new
    if oest is None:
        return astro.makesq(all)
    else:
        return (astro.makesq(all), oest)

def smamemo(meanmotion, altitude=False, forceenv=force.deffe, units=astro.prefunits):
    '''Find the semimajor axis from the mean motion; if mean motion is
    a number, units are presumed to be revolutions/sidereal day.

    Example of geosynchronous satellite semimajor axis
      tork.smamemo(1.0)
      <Quantity 42164.1696233 km>
    '''
    if type(meanmotion) is u.Quantity:
        memod = meanmotion
    else:
        memod = meanmotion*1.0*u.rev/u.sday
    sma1 = np.cbrt((forceenv['earthmu']*u.m**3/u.s**2)/memod**2)
    sma = sma1.decompose().to(units['length'], equivalencies=u.dimensionless_angles())
    if altitude:
        return sma-forceenv['earthrad']*u.m
    else:
        return sma
