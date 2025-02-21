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
    if events['eclipse']:
        ret['eclipsedet'] = EclipseDetector(forceenv['sun'], forceenv['sunrad'], forceenv['earth'])
        handled = ret['eclipsedet'].withUmbra().withHandler(ContinueOnEvent())
        logger = EventsLogger()
        loggeddet = logger.monitorDetector(handled)
        okprop.addEventDetector(loggeddet)

    # Propagate
    propagated = okprop.propagate(ork0.date, ork0.date.shiftedBy(astro.timesec(proptime)))
    if events['eclipse']:
        ret['ephgen'] = generator.getGeneratedEphemeris();
        loggedevents = logger.getLoggedEvents()
        ret['umbra'] = posvel.tsephem([ev.state.pvt()[0] for ev in loggedevents],
                               [ev.state.pvt()[1] for ev in loggedevents])
        ret['umbra']['entering'] = [ev.increasing for ev in loggedevents]
        ret['events'] = loggedevents
    else:
        ret = generator.getGeneratedEphemeris();
    return ret

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
        ecld = generator['eclipsedet']

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

# The following verision returns a tuple of (Quantity, Time), each with the same shape (number of rows)
def propagate2(generator, reltimes, include_init=True):
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
        ecld = generator['eclipsedet']

    rts = astro.timesec(reltimes)
    if isinstance(rts, collections.abc.Iterable):
        if include_init:
            reltimes = np.insert(reltimes, 0, 0.0)
        data = [propagate2(generator, rt)[0] for rt in reltimes]
        return (u.Quantity(np.asarray(data), data[0].unit),
                ork.posvel.okad(gen.getMinDate()) + astropy.time.TimeDelta(reltimes))
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
