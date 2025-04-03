# TLEs and their propagation
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
# from org.orekit.frames import FramesFactory # Will need to get the correct frame for TLEs; TEME?
from org.orekit.utils import PVCoordinatesProvider
from ..core import astro
from . import force
from . import convert

# See core/spacetrack.py for defining stclient

def tle_latest(stclient, satnum):
    '''Find the latest tle from the space-track.org. The first
    argument should define should be a
    spacetrack.base.SpaceTrackClient using a valid identity and
    password.'''
    sattle = stclient.tle_latest(norad_cat_id=satnum, ordinal=1, format='tle')
    sat2lines = sattle.split("\n")
    return TLE(sat2lines[0], sat2lines[1])

# senttle = tork.tle_latest(stclient, 41335)
# sent_goodpvt = tork.tleprop(senttle)
def tleprop(tle, proptime=0.0, forceenv=force.deffe):
    '''Propagate the two-line elements using SGP4'''
    propagator = TLEPropagator.selectExtrapolator(tle)
    return convert._pvt(propagator.getPVCoordinates(tle.getDate().shiftedBy(astro.timesec(proptime)), \
                                                    forceenv['celestialframe']))
