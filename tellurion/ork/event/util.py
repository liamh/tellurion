"""Common utility functions to be called by various specific events"""

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.propagation.events import EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent

from tellurion.ork import convert


def _make_evdet(propagator, detector, continue_prop=True):
    """Make an event detector and add it to the propagator"""
    logger = EventsLogger()
    if continue_prop: # Continue propagating even after first event is detected
        handled = detector.withHandler(ContinueOnEvent())
    else:
        handled = detector.withHandler()
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return logger

def _event_transition_table(loggers, column_label, state_labels, reftime="epoch"):
    """Create an ephemeris PVT with the event added to the aux dict"""

    def transition_pairs(string):
        """Create a list of successive state transitions of a list of
        pairs [up, down] of two-character state transitions from the
        single-character state designation,

         transition_pairs('ups')
         [['up', 'pu'], ['ps', 'sp']]
        """
        sp = [a+b for a, b in zip(string, string[1:])]
        return [[pr, pr[::-1]] for pr in sp]

    def _pvt_from_logger(logger, column_label, inc_dec_labels):
        """Create a PVT including columns having two-character
        transition label made from two one-character state labels
        """
        loggedevents = logger.getLoggedEvents()
        if loggedevents:
            pvt = convert._pvt([ev.getState() for ev in loggedevents])
            pvt.aux = {column_label: " ".join([inc_dec_labels[0] if ev.isIncreasing() else inc_dec_labels[1]
                                   for ev in loggedevents])}
            return pvt

    if loggers:
        trprs = transition_pairs(state_labels)
        pvtlist = []
        for (lg, idl) in zip(loggers, trprs):
            if lg is not None:
                pvtl = _pvt_from_logger(lg, column_label, idl)
                if pvtl is not None:
                    pvtlist.append(pvtl)
        if pvtlist:
            return pvtlist[0].merge(pvtlist[1:])

def _label_positive_count(detectors, labels):
    """Determine the appropriate event label by the number of positive counts among the detectors."""
    return(labels[sum(1 for x in detectors if x is not None and x > 0.0)])
