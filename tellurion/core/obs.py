"""Observations from a fixed location on earth"""

import collections
import dataclasses
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time

from tellurion.astro import units
from tellurion.astro import quant

# obsdict1 = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed', 'rtasc': 'angle', 'decl': 'angle'}

EarthObservationT = collections.namedtuple('EarthObservationT', 'loc obs time')

def azelrange(azim, elev, rnge, obsloc, obstime, unitlookup=units.prefunits):
    """Create an azimuth, elevation, and range observation; the values
    for azim, elev, rnge must be u.Quantity."""
    def sqaer(azim, elev, rnge, obstime):
        if obstime.isscalar and azim.isscalar and elev.isscalar:
            if rnge == None:
                sq = quant.sq_from_dict({'azim' : azim, 'elev' : elev}, True)
            else:
                sq = quant.sq_from_dict({'azim' : azim, 'elev' : elev, 'range' : rnge}, True)
        else:
            if rnge == None:
                sq = quant.vstack(tuple([sqaer(a, e, t) \
                                          for a, e, t in zip(azim,elev,obstime)]))
            else:
                sq = quant.vstack(tuple([sqaer(a, e, r, t) \
                                          for a, e, r, t in zip(azim,elev,rnge,obstime)]))
        return sq
    return EarthObservationT(obsloc, units.changeunits(sqaer(azim, elev, rnge, obstime)), obstime)
