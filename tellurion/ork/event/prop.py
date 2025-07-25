'''Add events to propagator and generate ephemeris'''

from org.orekit.propagation.events import AltitudeDetector
from tellurion.core import astro
from tellurion.ork import convert
from . import eclipse
from . import visibility

def _add(events, propagator, generator, proptime, forceenv, reftime, output):
    '''Add detectors and loggers for all events, propagate, then
    generate the event transition tables and add them to `generator`.'''

    # Altitude detector to stop propagating if too low; this should always be present
    propagator.addEventDetector(AltitudeDetector(float(events['altitude'].si.value), \
                                                 forceenv['sphalt']))

    # Make the detectors and loggers
    eccdls = eclipse._mkdetlog(events, propagator, forceenv)
    visdls = visibility._mkdetlog(events, propagator, forceenv)

    # Propagate
    generator['propfn'](generator['epoch'].shiftedBy(astro.timesec(proptime)))

    # Generate the event transition tables
    eclipse._gentrans(eccdls, generator, reftime, output)
    visibility._gentrans(visdls, events, generator, reftime, output)
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
