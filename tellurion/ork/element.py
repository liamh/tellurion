"""
Orbital elements in Orekit
"""

import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType, CircularOrbit, KeplerianOrbit, PositionAngleType
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates

from ..core import astro
from ..core import nquant
from ..core import posvel
from ..core import element
from . import force
from . import convert

###########################################
#### Element values from orbital state ####
###########################################

_eldict = element.sfdict([["sma", "semimajor axis", "length", u.meter, KeplerianOrbit.getA],
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
                          ["period", "orbital period", "time", u.second, KeplerianOrbit.getKeplerianPeriod]])

_elphystype = {key: value['phystype'] for key, value in _eldict.items()}

def elementval (orbstate, elt, earthrad=force.deffe["earthrad"].si.value):
    """
    Compute the orbital element from the orbital state
    Arguments
      orbstate:    Representation of orbital state in any form
      elt:  the orbital element desired, see list eldict.keys(); may be a list, e.g. ["sma", "ecc"]
      earthrad: the radius of the earth, necessary to provide for altitudes of perigee and apogee
    """
    return element.statefnval(_keplerianorbit(orbstate), elt, _eldict, earthrad)

def tselements(ephem, elements):
    '''Make a time series of selected orbital elements'''
    return TimeSeries(time=ephem.time,
                      data=[dict(zip(elements, elementval(ephrow, elements)))
                            for ephrow in ephem.pvt()])

########################################
####    Make Kepler element set     ####
########################################

def _keplerianorbit(oes, units=(astro.prefunits['length'], astro.prefunits['angle']),
                    forceenv=force.deffe):
    '''Make a org.orekit.orbits.KeplerianOrbit from anything'''
    if element.iskepels(oes):
        return _keporb_from_components(oes.els, oes.t, units, forceenv)
    elif type(oes) is KeplerianOrbit:
        return oes
    elif posvel.ispvtcart(oes):
        co = CartesianOrbit(convert._tspvc(oes),
                            forceenv['celestialframe'], forceenv['earthmu'].si.value)
        return OrbitType.KEPLERIAN.convertType(co)
    elif hasattr(oes,'pvt'):
        _keplerianorbit(oes.pvt(), units=units, forceenv=forceenv)
    else:
        raise ValueError('Cannot transform to Keplerian elements')

def _keporb_from_components(oes, epoch, units=(astro.prefunits['length'], astro.prefunits['angle']), fe=force.deffe):
    '''Make a org.orekit.orbits.KeplerianOrbit from orbital elements as a u.Quantity or Dict'''
    oessi = oes.si.value
    if 'ma' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.MEAN, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              convert._okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'].si.value)   # Sets the central attraction coefficient (m³/s²)
    elif 'ta' in oessi.dtype.names:
        return KeplerianOrbit(float(oessi['sma']), float(oessi['ecc']), float(oessi['inc']), \
                              float(oessi['argper']), float(oessi['raan']), \
                              float(oessi['ma']), PositionAngleType.TRUE, \
                              fe['celestialframe'], # The frame in which the parameters are defined (must be a pseudo-inertial frame)
                              convert._okad(epoch),   # Sets the date of the orbital parameters
                              fe['earthmu'].si.value)   # Sets the central attraction coefficient (m³/s²)
    else:
        raise ValueError('Time element (ma or ta) required in element set')

###############################
####  Transformations      ####
###############################

# Transformations between Cartesian state vector and Kepler elements

def kepler(object, forceenv=force.deffe, mean_time_element=True): # Add prefunits
    '''The Kepler element set from the Cartesian PVT or equivalent'''
    co = CartesianOrbit(convert._tspvc(object),
                        forceenv['celestialframe'], forceenv['earthmu'].si.value)
    ko = OrbitType.KEPLERIAN.convertType(co)
    if mean_time_element:
        elnames = element.kepeltma_names
    else:
        elnames = element.kepeltta_names
    if hasattr(object, 't'):
        dttm = object.t
    elif hasattr(object, 'time'):
        dttm = object.time
    kepels = dict(zip(elnames, elementval(ko, elnames)))
    return element.kepler(kepels, dttm)

def pvt(object, dttm=None):
    '''Make a Cartesian PVT from the object, or an ephemeris
    generator, which has an initial state. If the object is an
    elements set without a datetime, it must be supplied in `dttm`.'''
    if element.iskepels(object, True):
        return convert._pvt(_keplerianorbit(object))
    elif element.iskepels(object, False):
        return convert._pvt(_keplerianorbit(object, dttm))
    else:
        raise ValueError('Can only transform Kepler element sets')

def allplane(oesdict, forceenv=force.deffe, unitlookup=astro.prefunits):
    '''Generate all plane pairs (sma, ecc), (radper, radapo), (altper, altapo) from the first or last pairs; additionally, the mean motion can be substituted for semimajor axis in the first pair.

    Example 1, convert from altitudes of perigee and apogee to semimajor axis and eccentricity
    byalts = tell.kepler({"altper":160*u.km, "altapo":20250*u.km, \
                          "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg}, \
                          tell.abstime('2022-02-15T08:30:00'))
    smaecc = tork.allplane(byalts)
    smaecc[0]['sma'] # <Quantity 16583.13646 km>
    smaecc[0]['ecc'] # <Quantity 0.60573583>
    tell.iskepels(byalts) # False
    tell.iskepels(smaecc) # True

    Example 2, convert from semimajor axis and eccentricity to altitudes of perigee and apogee
    bysmaecc = {"sma":8000, "ecc":0.1, "inc":45, "argper": 120.0, "raan": 80.0, "ma": 0.0}
    alts = tork.allplane(bysmaecc)
    alts[0]['altper'] # <Quantity 821.86354 km>
    alts[0]['altapo'] # <Quantity 2421.86354 km>

    Example 3, define a geosynchronous orbit by mean motion
    geo = tork.allplane({"memo":1.0*u.rev/u.sday, "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg})

    Example 4, define a geosynchronous transfer orbit
    gto = tork.allplane({"altper": 350*u.km, "altapo": tork.sma(1.0,True), \
                          "ecc":0.0*u.dimensionless_unscaled, \
                          "inc":0.0*u.deg, "argper": 120.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    '''
    names = oesdict.keys()
    if ('sma' in names or 'memo' in names) and 'ecc' in names:  # OR PERIOD IN NAMES
        if 'memo' in names:
            smav = sma(oesdict['memo'], False, forceenv, unitlookup)
            new = {'sma': smav, 'radper': smav*(1-oesdict['ecc']), 'radapo': smav*(1+oesdict['ecc'])}
        else:
            smav = oesdict['sma']
            memo = np.sqrt(forceenv['earthmu'].to(unitlookup['gravconst'])/smav**3)\
                     .to(unitlookup['angular speed'], equivalencies=u.dimensionless_angles())
            new = {'memo': memo, 'radper': smav*(1-oesdict['ecc']), 'radapo': smav*(1+oesdict['ecc'])}
        new['altper'] = new['radper'] - forceenv['earthrad']
        new['altapo'] = new['radapo'] - forceenv['earthrad']
    elif 'altper' in names and 'altapo' in names:
        new = {'radper': oesdict['altper'] + forceenv['earthrad'], \
               'radapo': oesdict['altapo'] + forceenv['earthrad']}
        new['sma'] = (new['radapo']+new['radper'])/2
        new['ecc'] = (new['radapo']-new['radper'])/(new['radapo']+new['radper'])
    else:
        raise ValueError('Plane must be defined by either (sma, ecc) or (altper, altapo)')
    return oesdict | new

def sma(input, altitude=False, forceenv=force.deffe, unitlookup=astro.prefunits):
    '''Find the semimajor axis from the mean motion, orbital perioid,
    altitude, or specific energy; if `input` is a number, it is
    assumed to be a mean motion in revolutions/sidereal day.

    Example of geosynchronous satellite semimajor axis
      tork.sma(1.0)
      <Quantity 42164.1696233 km>
    Example of orbital period
      tork.sma(10000*u.s)
      <Quantity 10032.11910363 km>
    Example of specific energy
      tork.sma(-20*(u.km/u.s)**2)

    '''
    mu = forceenv['earthmu']
    if type(input) is u.Quantity:
        quant = input
    else:
        quant = input*1.0*u.rev/u.sday
    pdim=u.get_physical_type(quant)
    if pdim == 'angular speed':
        s1 = np.cbrt(mu/quant**2)
        s = s1.to(unitlookup['length'], equivalencies=u.dimensionless_angles())
        if altitude:
            return s-forceenv['earthrad']
        else:
            return s
    elif pdim == 'time':
        return sma(u.rev/quant, altitude, forceenv, unitlookup)
    elif pdim == 'specific energy':
        return (-mu/(2*quant)).to(unitlookup['length'])
    elif pdim == 'length':
        return quant + forceenv['earthrad']
    else:
        raise ValueError('Cannot convert quantity to semimajor axis')
