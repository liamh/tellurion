import collections
import collections.abc
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
from astropy.timeseries import TimeSeries
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
sunstate = ['umbra', 'penumbra>umbra', 'umbra>penumbra', 'penumbra', 'fullsun>penumbra', 'penumbra>fullsun', 'fullsun']

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
    An Orekit object that is passed to `propagate()` as the first argument, if no event other than altitude is included in `events`.
    A Dict of generator (`'ephgen'`) to pass to `propagate()`, EclipseDetector (`'eclipsedet'`), and a timetable of umbra transitions (`'umbra') if `'umbra'` is included in `events`.

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
        umbra = eclipse_transitions(logger_umb, ('umbra>penumbra', 'penumbra>umbra'))
        penumbra = eclipse_transitions(logger_pen, ('penumbra>fullsun', 'fullsun>penumbra'))
        ret['sun transition'] = (umbra, penumbra)
    else:
        ret = generator.getGeneratedEphemeris();
    return ret

# demoa.prop.eclipse.suntrans[0][1][0:5] # umbra transitions
 # datetime.datetime(2025, 1, 1, 0, 35, 36, 376969)
 # datetime.datetime(2025, 1, 1, 1, 30, 8, 972149)
 # datetime.datetime(2025, 1, 1, 2, 6, 1, 462955)
 # datetime.datetime(2025, 1, 1, 3, 0, 34, 113826)
 # datetime.datetime(2025, 1, 1, 3, 36, 26, 547286)
# demoa.prop.eclipse.suntrans[0][0]['sunlight'][0:5]
#  <Quantity [1., 2., 1., 2., 1.]>
# demoa.prop.eclipse.suntrans[1][1][0:5] # penumbra transitions
 # datetime.datetime(2025, 1, 1, 0, 35, 45, 142205)
 # datetime.datetime(2025, 1, 1, 1, 30, 0, 379634)
 # datetime.datetime(2025, 1, 1, 2, 6, 10, 229761)
 # datetime.datetime(2025, 1, 1, 3, 0, 25, 519718)
 # datetime.datetime(2025, 1, 1, 3, 36, 35, 315610)
# demoa.prop.eclipse.suntrans[1][0]['sunlight'][0:5]
#  <Quantity [4., 5., 4., 5., 4.]>

# Order these
 # umbra 0 at 0h0m
 # 1 @ datetime.datetime(2025, 1, 1, 0, 35, 36, 376969)
 # should be penumbra 3 at 0h35m40s
 # 4 @ datetime.datetime(2025, 1, 1, 0, 35, 45, 142205)
 # should be fullsun 6 at 1h00m
 # 5 @ datetime.datetime(2025, 1, 1, 1, 30, 0, 379634)
 # should be penumbra 3 at 1h30m5s
 # 2 @ datetime.datetime(2025, 1, 1, 1, 30, 8, 972149)
 # should be umbra 0 at 2h0m
 # 1 @ datetime.datetime(2025, 1, 1, 2, 6, 1, 462955)
 # 4 @ datetime.datetime(2025, 1, 1, 2, 6, 10, 229761)
 # 5 @ datetime.datetime(2025, 1, 1, 3, 0, 25, 519718)
 # 2 @ datetime.datetime(2025, 1, 1, 3, 0, 34, 113826)
 # 1 @ datetime.datetime(2025, 1, 1, 3, 36, 26, 547286)
 # 4 @ datetime.datetime(2025, 1, 1, 3, 36, 35, 315610)


# fnval > 0 ork.prop.propagate(demoa.prop.eclipse.genev, TimeDelta('1hr 25min').datetime.seconds*u.s)

def make_eclipsedet(propagator, forceenv, umbra):
    '''Make an eclipse detector for either umbra (`umbra=True`) or
    penumbra (`umbra=False`) and addit to the `propagator`.'''
    eclipsedet = EclipseDetector(forceenv['sun'], forceenv['sunrad'], forceenv['earth'])
    logger = EventsLogger()
    if umbra:
        handled = eclipsedet.withUmbra().withHandler(ContinueOnEvent())
    else:
        handled = eclipsedet.withPenumbra().withHandler(ContinueOnEvent())
    loggeddet = logger.monitorDetector(handled)
    propagator.addEventDetector(loggeddet)
    return (eclipsedet, logger)

def eclipse_transitions(logger, incdec_tuple): # ('umbra>penumbra', 'penumbra>umbra')
    loggedevents = logger.getLoggedEvents()
    pvts = [augmentpvt(ev.state.pvt(), ev.increasing, incdec_tuple, 'sunlight')
            for ev in loggedevents]
    pvs = [pvt[0] for pvt in pvts]
    tms = [pvt[1] for pvt in pvts]
    return (u.Quantity(np.asarray(pvs), pvs[0].unit), \
            astropy.time.Time([tm.value for tm in tms]))

# The following verision returns a tuple of (Quantity, Time), each with the same shape (number of rows)
def propagate2(generator, reltimes, include_init=True):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. If an list of relative
    times is given, an ephemeris table is returned; if reltimes is a
    single time, then a PVT is returned. If `include`_init is true,
    then include the initial PVT in the ephemeris table.
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
        umbd = generator['umbradet']
        pend = generator['penumbradet']

    if rtscalar:
        # This includes the value of the event function "pvut" = position, velocity, umbra and time
        ss = gen.propagate(gen.getMinDate().shiftedBy(float(reltimes.to(u.s).value)))
        pvt = ss.orbit.pvt()
        if 'umbd' in locals():
            # return augmentpvt(pvt, umbd.g(ss), ('penumbra', 'umbra'), 'sunlight')
            return augmentpvt(pvt, [umbd.g(ss), pend.g(ss)], None, 'sunlight')
        else:
            return pvt
    else:
        if include_init:
            reltimes = np.insert(reltimes, 0, 0.0)
        data = [propagate2(generator, rt)[0] for rt in reltimes]
        return (u.Quantity(np.asarray(data), data[0].unit),
                ork.posvel.okad(gen.getMinDate()) + astropy.time.TimeDelta(reltimes))

def augmentpvt(pvt, indicator, negposnames, fieldname):
    '''Augment the posvel with a field `fieldname` using the index in
    sunstate that matches the name in `negposnames` (a tuple of size
    2) dependening on whether the `indicator` is negative or
    positive.
    '''
    pvtdict = astro.splitsq(pvt[0])
    if type(indicator) is list:
        if indicator[0] > 0.0 and indicator[1] > 0.0:
            fv = sunstate.index('fullsun')*u.dimensionless_unscaled
        elif indicator[0]*indicator[1] < 0.0:
            fv = sunstate.index('penumbra')*u.dimensionless_unscaled
        else:
            fv = sunstate.index('umbra')*u.dimensionless_unscaled
    else: # Transition
        if indicator <= 0.0:
            fv = sunstate.index(negposnames[0])*u.dimensionless_unscaled
        else:
            fv = sunstate.index(negposnames[1])*u.dimensionless_unscaled
            pvtdict[fieldname] = fv
    return (astro.makesq(pvtdict), pvt[1])

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
