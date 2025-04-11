# TLEs and their propagation
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
# from org.orekit.frames import FramesFactory # Will need to get the correct frame for TLEs; TEME?
from org.orekit.utils import PVCoordinatesProvider
from ..core import astro
from . import force
from . import convert
from . import niprop
from . import eclipse

# import astropy.units as u
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sentgen = tork.SGP4gen(isssent['SENTINEL 3A'], 1*u.day, {'altitude': 125.0*u.km, 'eclipse': True, 'visibility': []})
# tork.propagate(sentgen, np.linspace(5.0*u.minute, 60.0*u.minute, 12), True)

def SGP4gen(meanels, proptime, events=niprop.defev, forceenv=force.deffe):
    '''Propagate mean elements using SGP4'''
    if hasattr(meanels, 'model') and meanels.model == 'SGP4':
        ret = {}
        tle = TLE(*meanels.tle)
        ret['epoch'] = tle.getDate()
        propagator = TLEPropagator.selectExtrapolator(tle)
        ret['propfn'] = lambda propto: propagator.getPVCoordinates(propto, forceenv['celestialframe'])
        ret['pvt0'] = convert._pvt(ret['propfn'](ret['epoch']))
        if events['eclipse']:
            (ret['umbradet'], logger_umb) = eclipse._make_eclipsedet(propagator, forceenv, True)
            (ret['penumbradet'], logger_pen) = eclipse._make_eclipsedet(propagator, forceenv, False)
            ret['propfn'](ret['epoch'].shiftedBy(astro.timesec(proptime)))
            ret['sun transition'] = eclipse._eclipse_transition_table(logger_umb, logger_pen)
        return ret
    else:
        raise ValueError("Can only propagate SGP4 mean elements with SGP4")
