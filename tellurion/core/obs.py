"""Observations from a fixed location on earth"""

import collections
import dataclasses
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time
from tellurion.core import astro
from tellurion.core import nquant
#from tellurion.core import element

obsdict1 = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed', 'rtasc': 'angle', 'decl': 'angle'}

EarthObservationT = collections.namedtuple('EarthObservationT', 'loc obs time')

def azelrange(azim, elev, range, obsloc, obstime, unitlookup=astro.prefunits):
    '''Create an azimuth, elevation, and range observation'''
    def sqaer(azim, elev, range, obstime):
        if obstime.isscalar and azim.isscalar and elev.isscalar:
            if range == None:
                sq = nquant.structquant({'azim' : azim, 'elev' : elev}, obsdict1)
            else:
                sq = nquant.structquant({'azim' : azim, 'elev' : elev, 'range' : range}, obsdict1)
        else:
            if range == None:
                sq = nquant.vstack(tuple([sqaer(a, e, t) \
                                          for a, e, t in zip(azim,elev,obstime)]))
            else:
                sq = nquant.vstack(tuple([sqaer(a, e, r, t) \
                                          for a, e, r, t in zip(azim,elev,range,obstime)]))
        return sq
    return EarthObservationT(obsloc, astro.changeunits(sqaer(azim, elev, range, obstime)), obstime)
