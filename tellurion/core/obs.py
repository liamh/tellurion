"""Observations from a fixed location on earth"""

import collections
import dataclasses
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import astropy.time

from tellurion.astro import units
from tellurion.astro import quant

obsdict1 = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed', 'rtasc': 'angle', 'decl': 'angle'}

EarthObservationT = collections.namedtuple('EarthObservationT', 'loc obs time')

def azelrange(azim, elev, rnge, obsloc, obstime, unitlookup=units.prefunits, convunits=units.prefunits):
    """Create an azimuth, elevation, and range observation; the values
    for azim, elev, rnge must be u.Quantity."""
    issc = obstime.isscalar and azim.isscalar and elev.isscalar
    if rnge == None:
        sq = quant.make_quantity({'azim' : azim, 'elev' : elev}, obsdict1, True, unitlookup)
    else:
        sq = quant.make_quantity({'azim' : azim, 'elev' : elev, 'range' : rnge}, obsdict1, True, unitlookup)
    sqc = quant.changeunits(sq, convunits)
    return EarthObservationT(obsloc, sqc, obstime)
