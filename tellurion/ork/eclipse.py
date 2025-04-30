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

def _do_eclipse(umbpen):
    return type(umbpen) is list and len(umbpen)==2 and any(umbpen)

def _sun_states(umbpen):
    sun_states = 'ups' # Umbra, penumbra, sun
    if umbpen[0] and umbpen[1]:
        return sun_states
    elif umbpen[0]:
        return sun_states[0] + sun_states[2]
    else:
        return sun_states[1] + sun_states[2]

def _umbra_penumbra(umbra_penumbra, propagator, gendict, proptime, forceenv, reftime, occluder='earth', \
                    output='et'):
    '''Define eclipse detectors for umbra, penumbra, or both. The
    first argument `umbra_penumbra` should be a two-element Boolean
    list defining which events to include.'''
    if _do_eclipse(umbra_penumbra):
        if umbra_penumbra[0]:
            (detector_umb, logger_umb) = _make_eclipsedet(propagator, forceenv, True, occluder)
        else:
            logger_umb = None
            detector_umb = None
        if umbra_penumbra[1]:
            (detector_pen, logger_pen) = _make_eclipsedet(propagator, forceenv, False, occluder)
        else:
            logger_pen = None
            detector_pen = None
        gendict['propfn'](gendict['epoch'].shiftedBy(astro.timesec(proptime)))
        loggers = [logger_umb, logger_pen]
        gendict['eclipsedet'] = [detector_umb, detector_pen]
        gendict['sun transition'] = _eclipse_transition_table(loggers, _sun_states(loggers), reftime, \
                                                              output=output)
    else:
        gendict['eclipsedet'] = [None, None]
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

def _eclipse_transition_table(loggers, state_chars, reftime='epoch', output='et'):
    pvt = event._event_transition_table(loggers, 'suntrans', state_chars, reftime)
    if output=='pvt':
        return pvt
    else:
        return pvt.ephemeris()

def _solar_illumination_state(detectors, spacecraft_state):
    '''A single character, one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial eclipse, or 's' (full sun).'''
    # https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EclipseDetector.html
    # g: Compute the value of the switching function. This
    # function becomes negative when entering the region of
    # shadow and positive when exiting.
    if _do_eclipse(detectors):
        if detectors[0]:
            umbsl = detectors[0].g(spacecraft_state)
        else:
            umbsl = None
        if detectors[1]:
            pensl = detectors[1].g(spacecraft_state)
        else:
            pensl = None
        return ('sunlight', event._label_positive_count([umbsl, pensl], _sun_states([umbsl, pensl])))
    else:
        return ()
