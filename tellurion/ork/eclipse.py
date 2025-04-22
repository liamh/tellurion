"""The `'eclipse'` dictionary entry in the `events` argument to
`prepare()` is a list of two booleans, indicating whether to detect
umbra, penumbra or both.  A transitions table will be created by
`prepare()` in `['sun transition']`, and subsequent propagation will
add a column `sunlight` that shows the sunlight state (`u`= umbra,
`p`=penumbra, `s`=full sun).
"""

from org.orekit.propagation.events import EclipseDetector
from ..core import astro
from . import event

def _sun_states(umbpen):
    sun_states = 'ups' # Umbra, penumbra, sun
    if umbpen[0] and umbpen[1]:
        return sun_states
    elif umbpen[0]:
        return sun_states[0] + sun_states[2]
    else:
        return sun_states[1] + sun_states[2]

def _remove_none(lst):
    return list(filter(lambda x: x is not None, lst))

def _umbra_penumbra(umbra_penumbra, propagator, gendict, proptime, forceenv, reftime, occluder='earth'):
    '''Define eclipse detectors for umbra, penumbra, or both. The
    first argument `umbra_penumbra` should be a two-element Boolean
    list defining which events to include.'''
    if umbra_penumbra[0]:
        (gendict['umbradet'], logger_umb) = _make_eclipsedet(propagator, forceenv, True, occluder)
    else:
        logger_umb = None
    if umbra_penumbra[1]:
        (gendict['penumbradet'], logger_pen) = _make_eclipsedet(propagator, forceenv, False, occluder)
    else:
        logger_pen = None
    gendict['propfn'](gendict['epoch'].shiftedBy(astro.timesec(proptime)))
    aloggers = _remove_none([logger_umb, logger_pen])
    gendict['sun transition'] = _eclipse_transition_table(aloggers, _sun_states([logger_umb, logger_pen]), reftime)
    return gendict

def _make_eclipsedet(propagator, forceenv, umbra, occluder='earth'):
    '''Make an eclipse detector for either umbra (`umbra=True`) or penumbra (`umbra=False`) and add it to the `propagator`.'''
    if umbra:
        detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                           forceenv[occluder]).withUmbra()
    else:
        detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                           forceenv[occluder]).withPenumbra()
    return (detector, event._make_evdet(propagator, detector, True))

def _eclipse_transition_table(loggers, state_chars, reftime='epoch'):
    return event._event_transition_table(loggers, 'suntrans', state_chars, reftime)

def _solar_illumination_state(umbra_detector, penumbra_detector, spacecraft_state):
    '''A single character, one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial eclipse, or 's' (full sun).'''
    # https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EclipseDetector.html
    # g: Compute the value of the switching function. This
    # function becomes negative when entering the region of
    # shadow and positive when exiting.
    if umbra_detector:
        umbsl = umbra_detector.g(spacecraft_state)
    else:
        umbsl = None
    if penumbra_detector:
        pensl = penumbra_detector.g(spacecraft_state)
    else:
        pensl = None
    return event._label_positive_count(_remove_none([umbsl, pensl]), _sun_states([umbsl, pensl]))
