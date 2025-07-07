"""Internal-only functions for three-state events"""

import operator
import itertools
import numpy as np
import astropy.units as u
import astropy.time
from org.orekit.propagation.events import EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent
from ..core import posvel
from ..core import astro
from . import convert

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
