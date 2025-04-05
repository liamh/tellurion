# TLEs and their propagation
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
# from org.orekit.frames import FramesFactory # Will need to get the correct frame for TLEs; TEME?
from org.orekit.utils import PVCoordinatesProvider
from ..core import astro
from . import force
from . import convert

# See core/spacetrack.py for defining stclient

# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sent_goodpvt = tork.tleprop(isssent[1][1])
def tleprop(tlestr, proptime=0.0, forceenv=force.deffe):
    '''Propagate the two-line elements using SGP4'''
    tle = TLE(*tlestr)
    propagator = TLEPropagator.selectExtrapolator(tle)
    return convert._pvt(propagator.getPVCoordinates(tle.getDate().shiftedBy(astro.timesec(proptime)), \
                                                    forceenv['celestialframe']))
