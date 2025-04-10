# TLEs and their propagation
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
# from org.orekit.frames import FramesFactory # Will need to get the correct frame for TLEs; TEME?
from org.orekit.utils import PVCoordinatesProvider
from ..core import astro
from . import force
from . import convert
from . import niprop
from . import eclipse

# See core/spacetrack.py for defining stclient

# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sent_goodpvt = tork.SGP4prop(isssent['SENTINEL 3A'])

def SGP4prop(meanels, proptime=0.0, forceenv=force.deffe, events=niprop.defev):
    '''Propagate mean elements using SGP4'''
    if hasattr(meanels, 'model') and meanels.model == 'SGP4':
        ret = {}
        tle = TLE(*meanels.tle)
        propto = tle.getDate().shiftedBy(astro.timesec(proptime))
        propagator = TLEPropagator.selectExtrapolator(tle)
        if events['eclipse']:
            (ret['umbradet'], logger_umb) = eclipse._make_eclipsedet(propagator, forceenv, True)
            (ret['penumbradet'], logger_pen) = eclipse._make_eclipsedet(propagator, forceenv, False)
        tspvc = propagator.getPVCoordinates(propto, forceenv['celestialframe'])
        ret['pvt'] = convert._pvt(tspvc)
        if events['eclipse']:
            ret['sun transition'] = eclipse._eclipse_transition_table(logger_umb, logger_pen)
        return ret
    else:
        raise ValueError("Can only propagate SGP4 mean elements with SGP4")
