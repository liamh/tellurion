"""Add events to propagator and generate ephemeris"""

import funcy

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.events import AltitudeDetector

from . import eclipse, visibility


def _add_pre(events, propagator, forceenv):
    """Add detectors and loggers for all events."""
    # Altitude detector to stop propagating if too low; this should always be present
    propagator.addEventDetector(AltitudeDetector(float(events["altitude"].si.value), \
                                                 forceenv["sphalt"]))
    # Make the detectors and loggers
    return {"eclipse": eclipse._mkdetlog(events, propagator, forceenv), \
            "visibility": visibility._mkdetlog(events, propagator, forceenv)}

def _add_post(detlogs, events, generator, reftime, output):
    """Generate the event transition tables and add them to `generator`."""
    eclipse._gentrans(detlogs["eclipse"], generator, reftime, output)
    visibility._gentrans(detlogs["visibility"], events, generator, reftime, output)
    return generator

# Replace _ephemeris with _evstates, which does not convert the pvt, just generates the event states
def _evstates(generator, spacecraft_states):
    """Create the aux dictionary with event states."""
    if type(spacecraft_states) is SpacecraftState:
        spacecraft_states = [spacecraft_states]
    def aa(event):
        if event:
            if type(event) is list:
                new = dict(event)
            else:
                new = {event[0]: event[1]}
        else:
            new = {}
        return new
    ed = generator["event detectors"]
    d = funcy.merge_with(" ".join, *[aa(eclipse._statechar(ed.get(eclipse._column_label), ss)) \
                                     for ss in spacecraft_states]) | \
        funcy.merge_with(" ".join, *[aa(visibility._statechar(ed.get(visibility._column_label), ss)) \
                                     for ss in spacecraft_states])
    return d
