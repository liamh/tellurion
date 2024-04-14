from orbit import *
from cartprodparam import *
import warnings

# Propagation and ephemeris
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector

# Propagate from epoch for a specified time
def mkephem(orbit, proptime, force, stopalt=125.0e3):
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

    initialState = SpacecraftState(orb, force['mass'])
    prop = NumericalPropagator(integrator)
    prop.setOrbitType(okc['cartesian'])
    prop.setInitialState(initialState)
    generator = prop.getEphemerisGenerator()

    # Forces
    prop.addForceModel(HolmesFeatherstoneAttractionModel(force['earthframe'], force['gravity']))
    if 'dragforce' in force:
        prop.addForceModel(force['dragforce'])

    # Events
    prop.addEventDetector(AltitudeDetector(stopalt, force['sphalt']))

    # Propagate
    pt = timesec(proptime) # Won't handle lists yet
    propagated = prop.propagate(orb.date, orb.date.shiftedBy(pt))
    ephemeris = generator.getGeneratedEphemeris();

    if type(orbit) is Ecis:
        orbit.update(prop=newtree('ephemeris',ephemeris))
        orbit.prop.maxtime = pt
        orbit.prop.default = force
        newname = f"prop{int(pt)}s"
        tree[newname] = orbit.pop('prop')
        return(tree[newname])
    else:
        return(ephemeris)

# Propagate to the relative time requested
# orbit: tree with orbit in it
# reltime: time (seconds) past the earliest time of the propagation
def prop(orbit, reltime, stopalt=125.0e3):
    [bp, tree] = thingofclass(orbit, BoundedPropagator)
    if tree is None:
        default = None
    else:
        default = tree.default
    maxtime = reltime # max(reltime)
    if bp is None:
        bp = mkephem(orbit, maxtime, default, stopalt)
        [bp, tree] = thingofclass(orbit, BoundedPropagator)
    #state = [bp.propagate(bp.getMinDate().shiftedBy(float(rt))).orbit for rt in reltime]
    state = bp.propagate(bp.getMinDate().shiftedBy(float(reltime))).orbit
#    if type(orbit) is Ecis:
#        name = "table" # f"state{int(reltime)}s"
#        tree[name] = state
    return(state)

############### [2024-04-07 Sun 22:41] This only works for scalar time, and only once
#from propagate import *
# ex1day = prop(ex1,86400.0)
# prop(ex1,43200.0) # show state at half day
# prop(ex1,86400.0) # show state at full day

# use timeseries to produce table of values



# from astropy.timeseries import TimeSeries
# ts1 = TimeSeries(time_start='2016-03-22T12:30:31', time_delta=3 * u.s, n_samples=5)


# >>>>>>>>>> Apply to list too
# Convert a Quantity to seconds as a Python float
def timesec(t):
    if type(t) is Quantity and get_physical_type(t) == 'time':
        pt = float(t.si.value) # convert to seconds and get the value_unit
    else:
        pt = float(t) # assume seconds
    return(pt)
