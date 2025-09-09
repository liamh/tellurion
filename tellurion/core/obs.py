"""Observations from a fixed location on earth"""

import collections
import dataclasses
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time
from tellurion.core import util
from tellurion.core import astro
from tellurion.core import element

obsdict1 = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed', 'rtasc': 'angle', 'decl': 'angle'}

#obsdict = element.sfdict([["azim", "azimuth", "angle", u.radian, fnazim],
#                          ["elev", "elevation", "angle", u.radian, fnelev]])

EarthObservationT = collections.namedtuple('EarthObservationT', 'loc obs time')

def azelrange(azim, elev, range, obsloc, obstime):
    '''Create an azimuth, elevation, and range observation'''
    if range == None:
        sq = astro.makesq({'azim' : azim, 'elev' : elev}, obsdict1)
    else:
        sq = astro.makesq({'azim' : azim, 'elev' : elev, 'range' : range}, obsdict1)
    return EarthObservationT(obsloc, sq, obstime)
