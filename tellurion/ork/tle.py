# TLEs and their propagation
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
# from org.orekit.frames import FramesFactory # Will need to get the correct frame for TLEs; TEME?
from org.orekit.utils import PVCoordinatesProvider
from ..core import astro
from . import force
from . import convert

# See core/spacetrack.py for defining stclient

# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sent_goodpvt = tork.SGP4prop(isssent['SENTINEL 3A'])

def SGP4prop(meanels, proptime=0.0, forceenv=force.deffe):
    '''Propagate mean elements using SGP4'''
    if hasattr(meanels, 'model') and meanels.model == 'SGP4':
        tle = TLE(*meanels.tle)
        propagator = TLEPropagator.selectExtrapolator(tle)
        return convert._pvt(propagator.getPVCoordinates(tle.getDate().shiftedBy(astro.timesec(proptime)), \
                                                        forceenv['celestialframe']))
    else:
        raise ValueError("Can only propagte SGP4 mean elements with SGP4")
