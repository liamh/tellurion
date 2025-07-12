"""The `'visibility'` dictionary entry in the `events` argument to
`prepare()` is a list of locations as produced by observer_location().
A transitions dictionary will be created by `prepare()` in
`['visibility']`, with the names of the observer locations serving as
indices and the start of the pass indicated by 'ny' (no visibility to
yes visibility) and end by 'yn' (yes visibility to no visibility).

See example `demod` in `demos/obsdemo.py`.
"""

import astropy.coordinates as coord
import astropy.units as u
from org.orekit.propagation.events import ElevationDetector
from org.orekit.frames import TopocentricFrame
from ..core import astro
from ..ork import geog
from . import event

############################################
### Required by the propagator (prop.py) ###
############################################

def _add(locations, propagator, gendict, proptime, forceenv, reftime, output='et'):
    '''Add detectors to the propagator giving visibility for an observer at the requested locations.

    `locations`: A list of dicts with the following keys
        `'location'`: The location of the observer (an astropy.coordinates.earth.EarthLocation)
        `'minelev'` : The minimum elevation (a u.Angle)
        `'name'`    : The name of the observing location (a string)
    '''

    def mkdet(obsloc, observer_body='earth'):
        '''Make a visibility detector for a single location and add it
        to the `propagator`.'''
        tf = TopocentricFrame(forceenv[observer_body], \
                              geog.geodpt(obsloc['location']), obsloc['name'])
        detector = ElevationDetector(tf).withConstantElevation(obsloc['minelev'].radian)
        return (detector, event._make_evdet(propagator, detector, True))

    def transtable(obsloc, logger):
        '''Add a column with the visibility transitions for a single observer location.'''
        colname = " ".join([obsloc['name'],_column_label])
        if logger:
            pvt = event._event_transition_table([logger], colname, _states, reftime)
            if pvt:
                if output=='pvt':
                    return pvt
                else:
                    return pvt.ephemeris()

    dls = [mkdet(ol) for ol in locations] # Add detectors for all locations
    gendict['propfn'](gendict['epoch'].shiftedBy(astro.timesec(proptime))) # Propagate
    names = [ol['name'] for ol in locations]
    gendict['visibilitydet'] = dict(zip(names, [dl[0] for dl in dls]))
    # Make a dict of the visibility transition ephemerides for all the locations
    gendict[_column_label] \
        = dict(zip(names, [transtable(ol, dl[1]) for (ol, dl) in zip(locations, dls)]))
    return gendict

def _statechar(detectors, spacecraft_state):
    '''A single character, one of 'n' or 'y', for each location at each time step.'''
    if detectors:
        dvis = []
        for name, det in detectors.items():
            visst = det.g(spacecraft_state)
            colname = " ".join([name, _column_label])
            dvis.append((colname, event._label_positive_count([visst], _states)))
        return dvis
    return None

############################################
### Internal definitions                 ###
############################################

_column_label = 'visibility'
_states = 'ny' # no/yes
