import collections
import collections.abc
import astropy.time
from astropy.timeseries import TimeSeries
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
import org.orekit.forces.gravity as okgrav

import astro
import orekit
import pvork
import posvel
import dttm

def generate(initpvt, proptime, force=pvork.deffe, stopalt=125.0e3):
    """
        Make a generator for an ephemeris; the output is the first argument of propagate()

        mkephem creates a BoundedPropagator but computes no actual states (see prop())
        initpvt:  the initial state
        proptime: the maximum time (s) to propagate
        force:    forces to use
        stopalt:  lowest altitude above spherical earth (m) to propagate
    """

    if posvel.ispv(initpvt):
        pvt0 = posvel.makepvt((initpvt, dttm.nowutc()))
    else:
        pvt0 = initpvt
    ork0 = pvork.orkpvt(*pvt0).cartesian()

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

    initialState = SpacecraftState(ork0, force['mass'])
    okprop = NumericalPropagator(integrator)
    okprop.setOrbitType(OrbitType.CARTESIAN)
    okprop.setInitialState(initialState)
    generator = okprop.getEphemerisGenerator()

    # Forces
    okprop.addForceModel(okgrav.HolmesFeatherstoneAttractionModel(force['earthframe'], force['gravity']))
    if 'dragforce' in force:
        okprop.addForceModel(force['dragforce'])

    # Events
    okprop.addEventDetector(AltitudeDetector(stopalt, force['sphalt']))

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
        times = collections.deque([pvork.pvtork(st)[1] for st in states])
        dat = collections.deque([pvork.pvtork(st)[0] for st in states])
        if include_init:
            pvt0 = pvork.pvtork(generator.initialState.pVCoordinates)
            dat.appendleft(pvt0[0])
            times.appendleft(pvt0[1])
        datdict = {posvel._eph_pos: [d[posvel._eph_pos] for d in dat], \
                   posvel._eph_vel: [d[posvel._eph_vel] for d in dat]}
        return TimeSeries(time=times, data=datdict)
    else:
        return pvork.pvtork(generator.propagate(generator.getMinDate().shiftedBy(rts)).orbit)
