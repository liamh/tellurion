from orbit import *
import warnings

# Propagation and ephemeris
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector

# Propagate from epoch for a specified time
def mkephem(orbit,proptime,spacecraft=None,stopalt=125.0e3):
    [orb, tree] = thingofclass(orbit, Orbit)

    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, orb, orb.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
    	       minstep,
    	       maxstep,
	       JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
	       JArray_double.cast_(tolerances[1]))
    integrator.setInitialStepSize(initStep)

    # Initialize the spacecraft state
    # satellite_mass = 100.0  # The models need a spacecraft mass, unit kg.
    if spacecraft is None:
        initialState = SpacecraftState(orb, 100.0)
    else:
        initialState = SpacecraftState(orb, spacecraft['mass'])
    prop = NumericalPropagator(integrator)
    prop.setOrbitType(okc['cartesian'])
    prop.setInitialState(initialState)
    generator = prop.getEphemerisGenerator()

    # Forces
    prop.addForceModel(HolmesFeatherstoneAttractionModel(forcedflt['earthframe'], forcedflt['gravity']))
    if spacecraft is not None:
        prop.addForceModel(spacecraft['dragforce'])

    # Events
    prop.addEventDetector(AltitudeDetector(stopalt, sphalt))

    # Propagate
    propagated = prop.propagate(orb.date, orb.date.shiftedBy(proptime))
    ephemeris = generator.getGeneratedEphemeris();

    if type(orbit) is Ecis:
        orbit.update(prop=newtree('ephemeris',ephemeris))
        orbit.prop.maxtime = proptime
        newname = f"prop{int(proptime)}s"
        tree[newname] = orbit.pop('prop')
        return(tree[newname])
    else:
        return(ephemeris)

# Propagate to the relative time requested
# orbit: tree with orbit in it
# reltime: time (seconds) past the earliest time of the propagation
def prop(orbit, reltime, spacecraft=None):
    [bp, tree] = thingofclass(orbit, BoundedPropagator)
    if bp is None:
        bp = mkephem(orbit, reltime, spacecraft)
        [bp, tree] = thingofclass(orbit, BoundedPropagator)
    state = bp.propagate(bp.getMinDate().shiftedBy(reltime)).orbit
    if type(orbit) is Ecis:
        name = f"state{int(reltime)}s"
        tree[name] = state
    return(state)
