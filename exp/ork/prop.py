import collections
import collections.abc
import astropy.time
from astropy.timeseries import TimeSeries
import orekit
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
import org.orekit.forces.gravity as okgrav

import astro
import posvel
import element
import ork.force
import ork.posvel

def generate(initstate, proptime, forceenv=ork.force.deffe, stopalt=125.0e3):
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
    An Orekit object that is passed to `propagate()` as the first argument

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
    okprop.addEventDetector(AltitudeDetector(stopalt, forceenv['sphalt']))

    # Propagate
    propagated = okprop.propagate(ork0.date, ork0.date.shiftedBy(astro.timesec(proptime)))
    ephgen = generator.getGeneratedEphemeris();

    return(ephgen)
# end mkephem

def propagate(generator, reltimes, include_init=True):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. If an list of relative
    times is given, an ephemeris table is returned; if reltimes is a
    single time, then a PVT is returned. If `include`_init is true,
    then include the initial PVT in the ephemeris table.
    '''
    rts = astro.timesec(reltimes)
    if isinstance(rts, collections.abc.Iterable):
        states = [generator.propagate(generator.getMinDate().shiftedBy(rt)).orbit for rt in rts]
        pvts = [st.pvt() for st in states]
        times = collections.deque([pvt[1] for pvt in pvts])
        dat = collections.deque([pvt[0] for pvt in pvts])
        if include_init:
            pvt0 = generator.initialState.pVCoordinates.pvt()
            dat.appendleft(pvt0[0])
            times.appendleft(pvt0[1])
        datdict = {posvel._eph_pos: [d[posvel._eph_pos] for d in dat], \
                   posvel._eph_vel: [d[posvel._eph_vel] for d in dat]}
        ts = TimeSeries(time=times, data=datdict)
        ts[posvel._eph_pos].info.format = posvel._pos_format
        ts[posvel._eph_vel].info.format = posvel._vel_format
        return ts
    else:
        return generator.propagate(generator.getMinDate().shiftedBy(rts)).orbit.pvt()

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
