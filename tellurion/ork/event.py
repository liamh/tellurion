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
    '''Create an ephemeris table with a column of transitions'''
    trans = [_event_transition_label(lg, ind, state_labels) \
             for (lg, ind) in zip(loggers, list(range(len(loggers)))) if lg is not None]
    merged = sorted(list(itertools.chain.from_iterable(trans)), key=operator.itemgetter(2))
    if len(merged) > 0:
        pvs = u.Quantity([m[0] for m in merged])
        times = astropy.time.Time([m[2] for m in merged])
        # TODO: port to .ephemeris()
        txyz = posvel.posxyz(posvel.tsephem(pvs, times))
        txyz[column_label] = [m[1] for m in merged]
        elapsed = [dt.quantity_str for dt in np.diff(txyz['time'])]
        elapsed.insert(0,'')
        txyz.add_column(elapsed, index=1, name='elapsed')
        if reftime is not None:
            astro.fromtime(txyz, reftime=reftime, copy=False)
        return txyz
    return None

def _spairs(string, reverse=False):
    '''Make a string of successive pairs of characters'''
    if reverse:
        str=string[::-1]
        return [a+b for a, b in zip(str, str[1:])][::-1]
    else:
        str=string
        return [a+b for a, b in zip(str, str[1:])]

def _event_transition_label(logger, ind, labels):
    '''A two-character transition label made from two one-character state labels'''
    loggedevents = logger.getLoggedEvents()
    def pvet(ev):
        '''A 3-tuple of posvel, event transition (2-character string with prior and posterior event state), and time.'''
        pvt = convert._pvt(ev.getState().getPVCoordinates())
        st = _spairs(labels, not ev.isIncreasing())[ind]
        return (pvt.pv, st, pvt.time)
    return [pvet(ev) for ev in loggedevents]

def _label_positive_count(detectors, labels):
    '''Determine the appropriate event label by the number of positive counts among the detectors.'''
    return(labels[sum(1 for x in detectors if x is not None and x > 0.0)])
