"""Internal-only functions for computing eclipsing and transitions in and out of eclipses"""

from org.orekit.propagation.events import EclipseDetector
from . import event

_sun_states = ['u', 'p', 's'] # Umbra, penumbra, sun

def _make_eclipsedet(propagator, forceenv, umbra, occluder='earth'):
    '''Make an eclipse detector for either umbra (`umbra=True`) or penumbra (`umbra=False`) and add it to the `propagator`.'''
    if umbra:
        detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                           forceenv[occluder]).withUmbra()
    else:
        detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                           forceenv[occluder]).withPenumbra()
    return (detector, event._make_evdet(propagator, detector, True))

def _eclipse_transition_table(logger_umb, logger_pen, reftime='epoch'):
    return event._event_transition_table([logger_umb, logger_pen], 'suntrans', _sun_states, reftime)

def _solar_illumination_state(umbra_detector, penumbra_detector, spacecraft_state):
    '''A single character, one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial eclipse, or 's' (full sun).'''
    # https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EclipseDetector.html
    # g: Compute the value of the switching function. This
    # function becomes negative when entering the region of
    # shadow and positive when exiting.
    umbsl = umbra_detector.g(spacecraft_state)
    pensl = penumbra_detector.g(spacecraft_state)
    return event._label_three_state([umbsl, pensl], _sun_states)
