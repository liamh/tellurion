"""Add events to the propagator"""

import operator
import itertools
import numpy as np
import astropy.units as u
import astropy.time
from org.orekit.propagation.events import EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent
from org.orekit.propagation.events import AltitudeDetector
from ..core import posvel
from ..core import astro
from . import convert
from . import eclipse
from . import visibility

############################################
### High level, called by prop()        ####
############################################

def _add(events, propagator, gendict, proptime, forceenv, reftime, output='et'):
    '''Add detectors and loggers for all events, propagate, then
    generate the event transition tables and add them to `gendict`.'''

    # Altitude detector to stop propagating if too low; this should always be present
    propagator.addEventDetector(AltitudeDetector(float(events['altitude'].si.value), \
                                                 forceenv['sphalt']))

    # Make the detectors and loggers
    eccdls = eclipse._mkdetlog(events, propagator, forceenv)
    visdls = visibility._mkdetlog(events, propagator, forceenv)

    # Propagate
    gendict['propfn'](gendict['epoch'].shiftedBy(astro.timesec(proptime)))

    # Generate the event transition tables
    eclipse._gentrans(eccdls, gendict, reftime, output)
    visibility._gentrans(visdls, events, gendict, reftime, output)
    return gendict

def _ephemeris(gendict, spacecraft_state):
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
    ed = gendict['event detectors']
    aa(eclipse._statechar(ed.get(eclipse._column_label), spacecraft_state))
    aa(visibility._statechar(ed.get(visibility._column_label), spacecraft_state))
    return pvt

############################################
### Low level, called by specific event ####
############################################

def _make_evdet(propagator, detector, continue_prop=True):
    '''Make an event detector and add it to the propagator'''
    logger = EventsLogger()
    if continue_prop: # Continue propagating even after first event is detected
        handled = detector.withHandler(ContinueOnEvent())
    else:
        handled = detector.withHandler()
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return logger

def _event_transition_table(loggers, column_label, state_labels, reftime='epoch'):
    '''Create an ephemeris PVT with the event added to the aux dict'''
    def transition_pairs(string):
        # _transition_pairs(state_labels)
        # [['up', 'pu'], ['ps', 'sp']]
        sp = [a+b for a, b in zip(string, string[1:])]
        return [[pr, pr[::-1]] for pr in sp]
    if loggers:
        trprs = transition_pairs(state_labels)
        pvtlist = []
        for (lg, idl) in zip(loggers, trprs):
            if lg is not None:
                pvtl = _pvt_from_logger(lg, column_label, idl)
                if pvtl is not None:
                    pvtlist.append(pvtl)
        if pvtlist:
            pvtcat = posvel.pvt(pvtlist)
            return pvtcat.timeorder()

def _pvt_from_logger(logger, column_label, inc_dec_labels):
    '''A two-character transition label made from two one-character state labels'''
    loggedevents = logger.getLoggedEvents()
    def pvet(ev):
        '''A 3-tuple of posvel, event transition (2-character string with prior and posterior event state), and time.'''
        pvt = convert._pvt(ev.getState().getPVCoordinates())
        if ev.isIncreasing():
            trlabel = inc_dec_labels[0]
        else:
            trlabel = inc_dec_labels[1]
        pvt.aux = {column_label:trlabel}
        return pvt
    if loggedevents:
        return posvel.pvt([pvet(ev) for ev in loggedevents])

def _label_positive_count(detectors, labels):
    '''Determine the appropriate event label by the number of positive counts among the detectors.'''
    return(labels[sum(1 for x in detectors if x is not None and x > 0.0)])
