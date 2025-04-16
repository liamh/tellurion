"""Internal-only functions for computing eclipsing and transitions in and out of eclipses"""

import operator
import numpy as np
import astropy.units as u
import astropy.time
from org.orekit.propagation.events import EclipseDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent
from ..core import posvel
from ..core import astro
from . import convert

def _make_eclipsedet(propagator, forceenv, umbra, occluder='earth'):
    '''Make an eclipse detector for either umbra (`umbra=True`) or penumbra (`umbra=False`) and add it to the `propagator`.'''
    logger = EventsLogger()
    if umbra:
        # Not necessary to have .withUmbra(), it is already set that way
        eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, forceenv[occluder]).withUmbra()
        handled = eclipsedet.withHandler(ContinueOnEvent())
    else:
        # Necessary to have withPenumbra(), as it is not changed in the instance
        eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, forceenv[occluder]).withPenumbra()
        handled = eclipsedet.withHandler(ContinueOnEvent())
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return (eclipsedet, logger)

def _eclipse_transition_table(logger_umb, logger_pen, reftime='epoch'):
    umbra = _eclipse_transitions(logger_umb, True)
    penumbra = _eclipse_transitions(logger_pen, False)
    suntr = sorted(umbra + penumbra, key=operator.itemgetter(2))
    if len(suntr) > 0:
        pvs = u.Quantity([m[0] for m in suntr])
        times = astropy.time.Time([m[2] for m in suntr])
        txyz = posvel.posxyz(posvel.tsephem(pvs, times))
        txyz['suntrans'] = [m[1] for m in suntr]
        elapsed = [dt.quantity_str for dt in np.diff(txyz['time'])]
        elapsed.insert(0,'')
        txyz.add_column(elapsed, index=1, name='elapsed')
        if reftime is not None:
            astro.fromtime(txyz, reftime=reftime, copy=False)
        return txyz
    return None

def _eclipse_transitions(logger, umbra):
    '''Find the transitions in and out of eclipse'''
    loggedevents = logger.getLoggedEvents()
    def suntrans(ev):
        if umbra:
            if ev.isIncreasing():
                return 'up' # Transition from umbra to penumbra
            else:
                return 'pu' # Transition from penumbra to umbra
        else:
            if ev.isIncreasing():
                return 'ps' # Transition from penumbra to full sunlight
            else:
                return 'sp' # Transition from full sunlight to penumbra
    def pvet(ev):
        '''A 3-tuple of posvel, sun transition (2-character string with prior and posterior sun state), and time.'''
        pvt = convert._pvt(ev.getState().getPVCoordinates())
        st = suntrans(ev)
        return (pvt[0], st, pvt[1])
    return [pvet(ev) for ev in loggedevents]

def _solar_illumination_state(umbra_detector, penumbra_detector, spacecraft_state):
    '''A single character, one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial eclipse, or 's' (full sun).'''
    # https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EclipseDetector.html
    # g: Compute the value of the switching function. This
    # function becomes negative when entering the region of
    # shadow and positive when exiting.
    umbsl = umbra_detector.g(spacecraft_state)
    pensl = penumbra_detector.g(spacecraft_state)
    if umbsl < 0.0 and pensl < 0.0:
        return 'u'
    elif umbsl*pensl < 0.0:
        return 'p'
    else:
        return 's'
