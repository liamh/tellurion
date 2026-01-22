"""The `'eclipse'` dictionary entry in the `events` argument to
`prepare()` is a list of two booleans, indicating whether to detect
umbra, penumbra or both.  A transitions table will be created by
`prepare()` in `['sun transition']`, and subsequent propagation will
add a column `sunlight` that shows the sunlight state (`u`= umbra,
`p`=penumbra, `s`=full sun).
"""

from org.orekit.propagation.events import EclipseDetector
from . import util

############################################
### Required by the propagator (prop.py) ###
############################################

def _mkdetlog(events, propagator, forceenv, occluder='earth'):
    '''Define eclipse detectors for umbra, penumbra, or both, and add
    them to the propagator. The first argument `umbra_penumbra` should
    be a two-element Boolean list defining which events to include.
    '''
    def mkdet(propagator, forceenv, umbra, occluder='earth'):
        '''Make an eclipse detector for either umbra (`umbra=True`) or
        penumbra (`umbra=False`) and add it to the `propagator`.'''
        if umbra:
            detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                               forceenv[occluder]).withUmbra()
        else:
            detector = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, \
                                               forceenv[occluder]).withPenumbra()
        return (detector, util._make_evdet(propagator, detector, True))

    umbra_penumbra = events[_column_label]
    if _do_eclipse(umbra_penumbra):
        if umbra_penumbra[0]:
            (detector_umb, logger_umb) = mkdet(propagator, forceenv, True, occluder)
        else:
            logger_umb = None
            detector_umb = None
        if umbra_penumbra[1]:
            (detector_pen, logger_pen) = mkdet(propagator, forceenv, False, occluder)
        else:
            logger_pen = None
            detector_pen = None
        return ((detector_umb, detector_pen), (logger_umb, logger_pen))

def _gentrans(detlogs, generator, reftime, output):
    '''Generate eclipse transitions.'''
    def transtable(loggers, state_chars, reftime, output):
        pvt = util._event_transition_table(loggers, _column_label, state_chars, reftime)
        if output=='pvt':
            return pvt
        elif hasattr(pvt,'ephemeris'):
            return pvt.ephemeris()

    if detlogs:
        generator['event detectors'][_column_label] = detlogs[0]
        generator['sun transition'] = transtable(detlogs[1], _sun_states(detlogs[1]), reftime, output)

def _statechar(detectors, spacecraft_state):
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
        return (_column_label, util._label_positive_count([umbsl, pensl], _sun_states([umbsl, pensl])))
    else:
        return ()

############################################
### Internal definitions                 ###
############################################

_column_label = 'eclipse'

def _do_eclipse(umbpen):
    return hasattr(umbpen,'__len__') and len(umbpen)==2 and any(umbpen)

def _sun_states(umbpen):
    sun_states = 'ups' # Umbra, penumbra, sun
    if umbpen[0] and umbpen[1]:
        return sun_states
    elif umbpen[0]:
        return sun_states[0] + sun_states[2]
    else:
        return sun_states[1] + sun_states[2]
