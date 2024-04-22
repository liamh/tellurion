##### Propagate orbits
## Main function: ephemeris() will generate an ephemeris of uniform step times from an orbit
##  ephemeris(orbit, tstep, nsteps):
## Returns a TimeSeries
##
## Example propagate ex1 for 10 minutes from epoch to epoch + 2 hours inclusive
##  ephemeris(ex1,10*u.min,12)
## Returns the ephemeris TimeSeries and saves in ex1.prop7200s.ephem600s12
## Get a single PVT from the sixth row
##  ex1.prop7200s.ephem600s12[5].posveltime()
## Get the whole time series as a list of PVTs
##  ex1.prop7200s.ephem600s12.posveltime()

from orbit import *
from cartprodparam import *
import warnings
import collections.abc

# Propagation and ephemeris
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
# Time series for ephemeris
from astropy.timeseries import TimeSeries

# Propagate from epoch for a specified time
# mkephem creates a BoundedPropagator but computes no actual states (see prop())
#  orbit:    the initial state
#  proptime: the maximum time (s) to propagate
#  force:    forces to use
#  stopalt:  lowest altitude above spherical earth (m) to propagate
def mkephem(orbit, proptime, force, stopalt=125.0e3):
    [orb, tree] = ecis.thingofclass(orbit, Orbit)

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
    pt = timesec(proptime)
    propagated = prop.propagate(orb.date, orb.date.shiftedBy(pt))
    ephgen = generator.getGeneratedEphemeris();

    if type(orbit) is ecis.Ecis:
        orbit.update(prop=ecis.newtree('ephgen',ephgen))
        orbit.prop.maxtime = pt
        orbit.prop.forceenv = force
        newname = f"prop{int(pt)}s"
        tree[newname] = orbit.pop('prop')
        return(tree[newname])
    else:
        return(ephgen)

# Propagate the BoundedPropagator
# orbit: tree with orbit in it
# reltimes: times (seconds) past the initial time; may be a number (seconds), Quantity with type 'time' or list of these
def prop(orbit, reltimes, stopalt=125.0e3):
    [bp, tree] = ecis.thingofclass(orbit, BoundedPropagator)
    if tree is None:
        forceenv = None
    else:
        forceenv = tree.forceenv
    if isinstance(reltimes, collections.abc.Iterable):
        maxtime = max(reltimes)
    else:
        maxtime = reltimes
    if bp is None:
        bp = mkephem(orbit, maxtime, forceenv, stopalt)
        [bp, tree] = ecis.thingofclass(orbit, BoundedPropagator)
    if isinstance(reltimes, collections.abc.Iterable):
        state = [bp.propagate(bp.getMinDate().shiftedBy(timesec(rt))).orbit for rt in reltimes]
    else:
        state = bp.propagate(bp.getMinDate().shiftedBy(timesec(reltimes))).orbit
    return(state)

# Compute an ephemeris table from an orbital state, assuming a
# BoundedPropagator exists in the tree.
def ephemeris(orbit, tstep, nsteps):
    tss = int(timesec(tstep))
    times = [timesec(dt) for dt in rangi(0,tss*nsteps,tss)]
    pvts = [porb.posveltime().pvt for porb in prop(ex1,times)]
    ts = TimeSeries(time=[pvt.time for pvt in pvts],
               data={'position': [pvt['p'] for pvt in pvts],
                     'velocity': [pvt['v'] for pvt in pvts]})
    ts['position'].info.format = '7.0f'
    ts['velocity'].info.format = '5.3f'
    [bp, tree] = ecis.thingofclass(orbit, BoundedPropagator)
    tree[f"ephem{tss}s{nsteps}"] = ts
    return ts

TimeSeries.posveltime = lambda ts: [ts[row].posveltime() for row in range(0,len(ts))]
astropy.table.row.Row.posveltime = \
    lambda row: PVT(row['position'], row['velocity'],row['time'],
                    units=[row.columns['position'].unit, row.columns['velocity'].unit])

# Convert a Quantity to seconds as a Python float
def timesec(t):
    if type(t) is Quantity and get_physical_type(t) == 'time':
        pt = float(t.si.value) # convert to seconds and get the value_unit
    else:
        pt = float(t) # assume seconds
    return(pt)
