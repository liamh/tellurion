'''Add events to propagator and generate ephemeris'''

from org.orekit.propagation.events import AltitudeDetector
from tellurion.core import astro
from tellurion.ork import convert
from . import eclipse
from . import visibility

def _add_pre(events, propagator, forceenv):
    '''Add detectors and loggers for all events.'''
    # Altitude detector to stop propagating if too low; this should always be present
    propagator.addEventDetector(AltitudeDetector(float(events['altitude'].si.value), \
                                                 forceenv['sphalt']))
    # Make the detectors and loggers
    return {'eclipse': eclipse._mkdetlog(events, propagator, forceenv), \
            'visibility': visibility._mkdetlog(events, propagator, forceenv)}

def _add_post(detlogs, events, generator, reftime, output):
    '''Generate the event transition tables and add them to `generator`.'''
    eclipse._gentrans(detlogs['eclipse'], generator, reftime, output)
    visibility._gentrans(detlogs['visibility'], events, generator, reftime, output)
    return generator

def _ephemeris(generator, spacecraft_state):
    '''Create the ephemeris with a column for each event.'''
    pvt = convert._pvt(spacecraft_state)
    def aa(event):
        if event:
            if type(event) is list:
                new = dict(event)
            else:
                new = {event[0]: event[1]}
        else:
            new = {}
        pvt.aux = pvt.aux | new
    ed = generator['event detectors']
    aa(eclipse._statechar(ed.get(eclipse._column_label), spacecraft_state))
    aa(visibility._statechar(ed.get(visibility._column_label), spacecraft_state))
    return pvt
