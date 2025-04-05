import operator
import collections
import collections.abc
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
import orekit_jpype as orekit
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit, KeplerianOrbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector, EclipseDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent
import org.orekit.forces.gravity as okgrav

from ..core import astro
from ..core import posvel
from ..core import element
from . import force
from . import element as oelement
from . import convert

defev = {'altitude': 125.0*u.km, 'eclipse': [], 'visibility': []}

def generate(initstate, proptime, forceenv=force.deffe, events=defev):
    """Make a generator for an ephemeris, optionally include eclipse
    information. The result of this function is passed as the first
    argument to `propagate()`.

    Parameters
    ----------
    initstate: a tell.ispvter() representing the initial state

    proptime:  u.Quantity, float
      The maximum time to propagate, numbers are in seconds

    forceenv:  dict
      Forces to use; output of force.setgravity()

    stopalt:   float
      Stop propagation if altitude above spherical earth (m) drops below this threshold

    Returns
    -------
    An Orekit object that is passed to `propagate()` as the first
    argument, if no event other than altitude is included in `events`.

    If eclipse detection is added to `events`, then a dictionary with
    the generator `['ephgen']`, and an ephemeris table `['sun
    transitions']` with xyz positions is returned. The latter has a
    two-character string 'suntrans', and elapsed time from the
    previous transition 'elapsed'. There are also two detectors used
    by `propagate()` in the dictionary.

    """

    ork0 = CartesianOrbit(convert._tspvc(*initstate), \
                          forceenv['celestialframe'], forceenv['earthmu'].si.value)

    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, ork0, ork0.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
        minstep,
        maxstep,
        tolerances[0],
        tolerances[1])
    integrator.setInitialStepSize(initStep)

    initialState = SpacecraftState(ork0, forceenv['mass'])
    okprop = NumericalPropagator(integrator)
    okprop.setOrbitType(OrbitType.CARTESIAN)
    okprop.setInitialState(initialState)
    generator = okprop.getEphemerisGenerator()

    # Forces
    okprop.addForceModel(okgrav.HolmesFeatherstoneAttractionModel(forceenv['earthframe'], forceenv['gravity']))
    if 'dragforce' in forceenv:
        okprop.addForceModel(forceenv['dragforce'])

    # Events
    okprop.addEventDetector(AltitudeDetector(float(events['altitude'].si.value), forceenv['sphalt']))
    ret = {}

    if events['eclipse']: # Only includes earth as occluding body
        (ret['umbradet'], logger_umb) = _make_eclipsedet(okprop, forceenv, True)
        (ret['penumbradet'], logger_pen) = _make_eclipsedet(okprop, forceenv, False)

    # Propagate
    propagated = okprop.propagate(ork0.getDate(), ork0.getDate().shiftedBy(astro.timesec(proptime)))

    if events['eclipse']:
        # Return a dictionary with the generator ['ephgen'], and an
        # ephemeris table ['sun transitions'] with xyz positions,
        # two-character string 'suntrans', and elapsed time from the
        # previous transition 'elapsed'
        ret['ephgen'] = generator.getGeneratedEphemeris();
        umbra = _eclipse_transitions(logger_umb, True)
        penumbra = _eclipse_transitions(logger_pen, False)
        suntr = sorted(umbra + penumbra, key=operator.itemgetter(2))
        pvs = u.Quantity([m[0] for m in suntr])
        times = astropy.time.Time([m[2] for m in suntr])
        txyz = posvel.posxyz(posvel.tsephem(pvs, times))
        txyz['suntrans'] = [m[1] for m in suntr]
        elapsed = [dt.quantity_str for dt in np.diff(txyz['time'])]
        elapsed.insert(0,'')
        txyz['elapsed'] = elapsed
        ret['sun transition'] = txyz
    else:
        ret = generator.getGeneratedEphemeris();
    return ret

def propagate(generator, reltimes, include_init=True, spacecraftstate=False):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. The relative times must
    satisfy posvel.isreltime(reltimes), and if the size
    prop5m1h.shape[0] > 0, an ephemeris table is returned. If reltimes
    is a single time, then a PVT is returned. If `include`_init is
    true, then include the initial PVT in the ephemeris table.

    If an eclipse detector has been added, the `'sunlight'` value will
    be one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial
    eclipse, or 's' (full sun).
    '''
    if posvel.isreltime(reltimes):
        rtshape = reltimes.shape
        rtscalar = rtshape == ()
    else:
        raise ValueError("Reltimes must be a relative time or times: a u.Quantity with physical type 'time'")

    if type(generator) is dict:
        gen = generator['ephgen']
        umbd = generator['umbradet']
        pend = generator['penumbradet']
    else:
        gen = generator

    if rtscalar:
        # This includes the value of the event function "pvut" = position, velocity, umbra and time
        ss = gen.propagate(gen.getMinDate().shiftedBy(float(reltimes.to(u.s).value)))
        pvt = convert._pvt(ss)
        if 'umbd' in locals():
            # https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EclipseDetector.html
            # g: Compute the value of the switching function. This
            # function becomes negative when entering the region of
            # shadow and positive when exiting.
            umbsl = umbd.g(ss)
            pensl = pend.g(ss)

            if umbsl < 0.0 and pensl < 0.0:
                sunstate = 'u'
            elif umbsl*pensl < 0.0:
                sunstate = 'p'
            else:
                sunstate = 's'
            return pvt + (sunstate,)
        elif spacecraftstate:
            return ss
        else:
            return pvt
    else:
        if include_init:
            reltimes = np.insert(reltimes, 0, 0.0)
        if spacecraftstate:
            return [propagate(generator, rt, False, True) for rt in reltimes]
        data = [propagate(generator, rt, False, False) for rt in reltimes]
        pvs = [d[0] for d in data]
        times = [d[1] for d in data]
        pvsq = u.Quantity(np.asarray(pvs), pvs[0].unit)
        if 'umbd' in locals():
            ephem = posvel.posxyz(posvel.tsephem(pvsq, times))
            ephem['sunlight'] = [d[2] for d in data]
        else:
            ephem = posvel.tsephem(pvsq, times)
        return ephem

def _make_eclipsedet(propagator, forceenv, umbra):
    '''Make an eclipse detector for either umbra (`umbra=True`) or penumbra (`umbra=False`) and add it to the `propagator`.'''
    logger = EventsLogger()
    if umbra:
        # Not necessary to have .withUmbra(), it is already set that way
        eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, forceenv['earth']).withUmbra()
        handled = eclipsedet.withHandler(ContinueOnEvent())
    else:
        # Necessary to have withPenumbra(), as it is not changed in the instance
        eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'].si.value, forceenv['earth']).withPenumbra()
        handled = eclipsedet.withHandler(ContinueOnEvent())
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return (eclipsedet, logger)

def _eclipse_transitions(logger, umbra):
    '''Find the transitions in and out of eclipse'''
    loggedevents = logger.getLoggedEvents()
    def suntrans(ev):
        if umbra:
            if ev.isIncreasing():
                return 'up' # Transition from umbra to penumbra
            else:
                return 'pu' # Transition from penumbra to umbra
        else:
            if ev.isIncreasing():
                return 'ps' # Transition from penumbra to full sunlight
            else:
                return 'sp' # Transition from full sunlight to penumbra
    def pvet(ev):
        '''A 3-tuple of posvel, sun transition (2-character string with prior and posterior sun state), and time.'''
        pvt = convert._pvt(ev.getState().getPVCoordinates())
        st = suntrans(ev)
        return (pvt[0], st, pvt[1])
    return [pvet(ev) for ev in loggedevents]

def timerange(object):
    '''The time difference between the earliest (usually the initial
    time) and the latest; not always what is requested as atmospheric
    drag can shorten the timespan
    '''
    if hasattr(object,'getMaxDate') and hasattr(object,'getMinDate'):
        return (convert._okad(object.getMaxDate())-convert._okad(object.getMinDate())).to(u.s)
    else:
        raise ValueError('Cannot get timerange for this object')
