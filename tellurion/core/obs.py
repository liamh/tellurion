"""Observations from a fixed location on earth"""

import collections
import dataclasses
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time
from tellurion.core import util
from tellurion.core import astro

@dataclasses.dataclass
class EarthObservation(util.QuantT):
    '''An observation or time series of observations from the earth'''
    location: coord.EarthLocation
    # observation: u.Quantity
    # '''The orbital state vector as an astropy.units 6-vector with structured quantity of physical dimension length, speed'''
    # time: astropy.time.Time
    # '''The date and time of the state'''

    def __init__(self, q, time, location):
        super().__init__(self)
        self.location = location

    def __getitem__(self, index):
        qt = super().__getitem__(self, index)
        rqt = util.QuantT(qt.q, qt.time, self.location[index])



obsdict = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed', 'rtasc': 'angle', 'decl': 'angle'}


def azelrange(azim, elev, range, obsloc, obstime):
    '''Create an azimuth, elevation, and range observation'''
    if range == None:
        sq = astro.makesq({'azim' : azim, 'elev' : elev}, obsdict)
    else:
        sq = astro.makesq({'azim' : azim, 'elev' : elev, 'range' : range}, obsdict)
    return EarthObservation(obsloc, sq, obstime)
