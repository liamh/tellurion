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

def mkephem(orbitinit, proptime, force, stopalt=125.0e3):
    """
        Time series for ephemeris

        Propagate from epoch for a specified time
        mkephem creates a BoundedPropagator but computes no actual states (see prop())
        orbitinit:    the initial state
        proptime: the maximum time (s) to propagate
        force:    forces to use
        stopalt:  lowest altitude above spherical earth (m) to propagate
    """

    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, orbitinit, orbitinit.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
        minstep,
        maxstep,
	    orekit.JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be cast in Python
	    orekit.JArray_double.cast_(tolerances[1]))
    integrator.setInitialStepSize(initStep)

    initialState = SpacecraftState(orbitinit, force['mass'])
    okprop = NumericalPropagator(integrator)
    okprop.setOrbitType(OrbitType.CARTESIAN)
    okprop.setInitialState(initialState)
    generator = okprop.getEphemerisGenerator()

    # Forces
    okprop.addForceModel(okgrav.HolmesFeatherstoneAttractionModel(force['earthframe'], force['gravity']))
    if 'dragforce' in force:
        okprop.addForceModel(force['dragforce'])

    # Events
    # okprop.addEventDetector(AltitudeDetector(stopalt, force['sphalt']))

    # Propagate
    propagated = okprop.propagate(orbitinit.date, orbitinit.date.shiftedBy(astro.timesec(proptime)))
    ephgen = generator.getGeneratedEphemeris();

    return(ephgen)
# end mkephem

def ephem(initpvt, reltimes, forceenv=pvork.deffe):
    if type(initpvt) is tuple:
        initcart = pvork.orkpvt(*initpvt).cartesian()
    elif posvel.isephrow(initpvt):
        return ephem((posvel.makepv(initpvt[posvel._eph_pos], initpvt[posvel._eph_vel]), \
                      initpvt['time']), \
                     reltimes, forceenv)
    else:
        initcart = pvork.orkpvt(initpvt, dttm.nowutc()).cartesian()
    bp = mkephem(initcart, max(reltimes), forceenv)
    states = [bp.propagate(bp.getMinDate().shiftedBy(astro.timesec(rt))).orbit \
              for rt in reltimes]
    times = [pvork.pvtork(st)[1] for st in states]
    dat = [pvork.pvtork(st)[0] for st in states]
    datdict = {posvel._eph_pos: [d['p'] for d in dat], posvel._eph_vel: [d['v'] for d in dat]}
    return TimeSeries(time=times, data=datdict)
