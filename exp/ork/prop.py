import operator
import collections
import collections.abc
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
import orekit
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector, EclipseDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent
import org.orekit.forces.gravity as okgrav

import astro
import posvel
import element
import ork.force
import ork.posvel

defev = {'altitude': 125.0*u.km, 'eclipse': [], 'visibility': []}

def generate(initstate, proptime, forceenv=ork.force.deffe, events=defev):
    """Make a generator for an ephemeris

    Parameters
    ----------
    initstate: u.Quantity or tuple of (u.Quantity, astropy.time.Time)
      The initial state, either a state (posvel or element set) or a
      tuple of (state, epoch). If epoch time of the initial state is not
      specified, the current time is used.

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

    A Dict of generator (`'ephgen'`) to pass to `propagate()`,
    EclipseDetector (`'eclipsedet'`), and a timetable of sun
    transitions `'sun transition'`, which are three-digit numbers
    beginning with 1, the second digit is the prior sun state (0 =
    total eclipse, 1 = partial eclipse, 2 = full sun), and the third
    digit is the posterior sun state.
    """

    if posvel.ispvt(initstate):
        pvt0 = initstate
        ork0 = ork.posvel.tspvc(*pvt0).cartesianorbit(forceenv)
    elif posvel.ispv(initstate):
        pvt0 = posvel.pvt((initstate, posvel.nowutc()))
        ork0 = ork.posvel.tspvc(*pvt0).cartesianorbit(forceenv)
    elif element.iskepels(initstate):
        if type(initstate) is tuple:
            ork0 = ork.element.keplerianorbit(*initstate)
        else:
            ork0 = ork.element.keplerianorbit(initstate, posvel.nowutc())
    elif hasattr(initstate, 'cartesianorbit'):
        ork0 = initstate.cartesianorbit()
    else:
        raise ValueError('Cannot propagate initstate')

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
	    orekit.JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be cast in Python
	    orekit.JArray_double.cast_(tolerances[1]))
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
        (ret['umbradet'], logger_umb) = make_eclipsedet(okprop, forceenv, True)
        (ret['penumbradet'], logger_pen) = make_eclipsedet(okprop, forceenv, False)

    # Propagate
    propagated = okprop.propagate(ork0.date, ork0.date.shiftedBy(astro.timesec(proptime)))
    if events['eclipse']:
        ret['ephgen'] = generator.getGeneratedEphemeris();
        umbra = eclipse_transitions(logger_umb, 'umbra')
        penumbra = eclipse_transitions(logger_pen, 'penumbra')
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

# fnval > 0 ork.prop.propagate(demoa.prop.eclipse.genev, TimeDelta('1hr 25min').datetime.seconds*u.s)

def make_eclipsedet(propagator, forceenv, umbra):
    '''Make an eclipse detector for either umbra (`umbra=True`) or
    penumbra (`umbra=False`) and add it to the `propagator`.'''
    eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'], forceenv['earth'])
    logger = EventsLogger()
    if umbra:
        handled = eclipsedet.withUmbra().withHandler(ContinueOnEvent())
    else:
        handled = eclipsedet.withPenumbra().withHandler(ContinueOnEvent())
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return (eclipsedet, logger)

def eclipse_transitions(logger, which):
    loggedevents = logger.getLoggedEvents()

    def suntrans(ev):
        if which == 'penumbra':
            if ev.increasing:
                return 'ps' # Transition from penumbra to full sunlight
            else:
                return 'sp' # Transition from full sunlight to penumbra
        else:  # umbra
            if ev.increasing:
                return 'up' # Transition from umbra to penumbra
            else:
                return 'pu' # Transition from penumbra to umbra

    def pvet(ev):
        '''A 3-tuple of posvel, sun transition (2-character string with prior and posterior sun state), and time.'''
        pvt = ev.state.pvt()
        st = suntrans(ev)
        return (pvt[0], st, pvt[1])

    return [pvet(ev) for ev in loggedevents]

# The following verision returns a tuple of (Quantity, Time), each with the same shape (number of rows)
def propagate2(generator, reltimes, include_init=True):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. If an list of relative
    times is given, an ephemeris table is returned; if reltimes is a
    single time, then a PVT is returned. If `include`_init is true,
    then include the initial PVT in the ephemeris table.

    If an eclipse detector has been added, the `'sunlight'` value will
    be one of 0.0 (umbra, or total eclipse), 1.0 (penumbra, or partial
    eclipse, or 2.0 (full sun).
    '''
    if type(reltimes) is u.Quantity and u.get_physical_type(reltimes) == 'time':
        rtshape = reltimes.shape
        rtscalar = rtshape == ()
    else:
        raise("Reltimes must be a u.Quantity with physical type 'time'")

    if type(generator) is BoundedPropagator:
        gen = generator
    else:
        gen = generator['ephgen']
        umbd = generator['umbradet'].withUmbra()  # Not necessary to have .withUmbra(), it is already set that way
        pend = generator['penumbradet'].withPenumbra() # Necessary to have withPenumbra(), as it is not changed in the instance

    if rtscalar:
        # This includes the value of the event function "pvut" = position, velocity, umbra and time
        ss = gen.propagate(gen.getMinDate().shiftedBy(float(reltimes.to(u.s).value)))
        pvt = ss.orbit.pvt()
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
        else:
            return pvt
    else:
        if include_init:
            reltimes = np.insert(reltimes, 0, 0.0)
        data = [propagate2(generator, rt) for rt in reltimes]
        pvs = [d[0] for d in data]
        pvsq = u.Quantity(np.asarray(pvs), pvs[0].unit)
        times = ork.posvel.okad(gen.getMinDate()) + astropy.time.TimeDelta(reltimes)
        ephem = posvel.posxyz(posvel.tsephem(pvsq, times))
        ephem['sunlight'] = [d[2] for d in data]
        return ephem

########################################
### OLD VERSION pre-umbra handling
########################################
def propagate(generator, reltimes, include_init=True):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. If an list of relative
    times is given, an ephemeris table is returned; if reltimes is a
    single time, then a PVT is returned. If `include`_init is true,
    then include the initial PVT in the ephemeris table.
    '''
    if type(generator) is BoundedPropagator:
        gen = generator
    else:
        gen = generator['ephgen']
        ecld = generator['umbradet']

    rts = astro.timesec(reltimes)
    if isinstance(rts, collections.abc.Iterable):
        states = [gen.propagate(gen.getMinDate().shiftedBy(rt)) for rt in rts]
        pvts = [st.orbit.pvt() for st in states]
        times = collections.deque([pvt[1] for pvt in pvts])
        dat = collections.deque([pvt[0] for pvt in pvts])
        if include_init:
            pvt0 = generator.initialState.pVCoordinates.pvt()
            dat.appendleft(pvt0[0])
            times.appendleft(pvt0[1])
        return posvel.tsephem(dat, times)
    else:
        # This includes the value of the event function "pvut" = position, velocity, umbra and time
        ss = gen.propagate(gen.getMinDate().shiftedBy(rts))
        pvt = ss.orbit.pvt()
        if 'ecld' in locals():
            pvtdict = astro.splitsq(pvt[0])
            pvtdict['umbra'] = ecld.g(ss)*u.dimensionless_unscaled
            return (astro.makesq(pvtdict), pvt[1])
        else:
            return pvt





# The following return AstroPy objects
SpacecraftState.pvt = lambda self, unitlookup=astro.prefunits: ork.posvel._pvtork(self.pVCoordinates, unitlookup)
BoundedPropagator.pvt = lambda self: self.initialState.pvt()
SpacecraftState.kepler = lambda self: self.orbit.kepler()
BoundedPropagator.kepler  = lambda self: self.initialState.kepler()
# The following return Orekit objects
SpacecraftState.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
BoundedPropagator.cartesianorbit  = lambda self: self.initialState.cartesianorbit()
SpacecraftState.keplerianorbit = lambda self: self.orbit.keplerianorbit()
BoundedPropagator.keplerianorbit  = lambda self: self.initialState.keplerianorbit()

# The time difference between the earliest (usually the initial time)
# and the latest; not always what is requested as atmospheric drag can
# shorten the timespan
BoundedPropagator.timerange = lambda self: (ork.posvel.okad(self.maxDate)-ork.posvel.okad(self.minDate)).to(u.s)
