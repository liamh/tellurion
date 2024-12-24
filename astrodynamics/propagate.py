"""Propagate orbits; main function: prop()"""

import astropy.units as u
import astropy.timeseries as apts
import astropy.table as aptbl
import orekit
import numpy as np
import warnings
import collections.abc

from . import cartprodparam
from . import ecis
from . import orbit

# Propagation and ephemeris
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
import org.orekit.forces.gravity as okgrav


def prop(orbitinit, reltimes, stopalt=125.0e3, maximum_tof=0.0):
    """Propagate an orbit
    orbitinit: tree with orbit in it
    reltimes: times (seconds) past the initial time; may be a number (seconds), Quantity with type 'time' or list of these

      The set of propagated points can be added to at any time, but be sure to set maximum_tof
      to the highest possible value when prop() is first called; otherwise, previously computed
      values might be lost. See example.propdemo().
      The class variable .pvt will have a list of pvts
      The class variable .orbit will have a list of orbits
      The class variable .ephem will have an AstroPy time series ephemeris table
      Pure numpy (no AstroPy or Orekit) is obtained from the ephemeris table with .makenp()
    """
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
        [orb, tree] = ecis.thingofclass(orbitinit, Orbit)

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
    	       orekit.JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be cast in Python
    	       orekit.JArray_double.cast_(tolerances[1]))
        integrator.setInitialStepSize(initStep)

        initialState = SpacecraftState(orb, force['mass'])
        prop = NumericalPropagator(integrator)
        prop.setOrbitType(orbit._okc['cartesian'])
        prop.setInitialState(initialState)
        generator = prop.getEphemerisGenerator()

        # Forces
        prop.addForceModel(okgrav.HolmesFeatherstoneAttractionModel(force['earthframe'], force['gravity']))
        if 'dragforce' in force:
            prop.addForceModel(force['dragforce'])

        # Events
        prop.addEventDetector(AltitudeDetector(stopalt, force['sphalt']))

        # Propagate
        pt = timesec(proptime)
        propagated = prop.propagate(orb.date, orb.date.shiftedBy(pt))
        ephgen = generator.getGeneratedEphemeris();

        if type(orbitinit) is ecis.Ecis:
            orbitinit.update(prop=ecis.newtree('ephgen',ephgen))
            orbitinit.prop.maxtime = pt
            orbitinit.prop.forceenv = force
            newname = "prop" #f"prop{int(pt)}s"
            tree[newname] = orbitinit.pop('prop')
            return(tree[newname])
        else:
            return(ephgen)
    # end mkephem

    reltimes = timearray(reltimes)
    maxtof = timearray(maximum_tof)
    if isinstance(reltimes, collections.abc.Iterable):
        maxtime = max(max(reltimes), maxtof)
    else:
        maxtime = max(reltimes, maxtof)
    [bp, tree] = ecis.thingofclass(orbitinit, BoundedPropagator)
    if tree is None:
        forceenv = None
    else:
        forceenv = tree.forceenv
    if bp is None or ("maxtime" in tree and type(tree.maxtime) is u.Quantity and maxtime > tree.maxtime):
        bp = mkephem(orbitinit, maxtime, forceenv, stopalt)
        [bp, tree] = ecis.thingofclass(orbitinit, BoundedPropagator)
    if isinstance(reltimes, collections.abc.Iterable):
        state = [bp.propagate(bp.getMinDate().shiftedBy(timesec(rt))).orbit for rt in reltimes]
        pvt=[orbit.PVT(st) for st in state]
    else:
        state = bp.propagate(bp.getMinDate().shiftedBy(timesec(reltimes))).orbit
        pvt=orbit.PVT(state)
    pvts = cartprodparam.ensurelist(tree.get('pvt')) + cartprodparam.ensurelist(pvt)
    pvts.sort(key = lambda s: s.pvtq.time)
    tree.update(pvt = pvts)
    states = cartprodparam.ensurelist(tree.get('orbit')) + cartprodparam.ensurelist(state)
    states.sort(key = lambda s: orbit.PVT(s).pvtq.time)
    tree.update(orbit=states)
    tree.update(ephem=ephts(tree))
    return(pvt)

# Make a Quantity with array value and time unit
def timearray(times):
    if type(times) is u.Quantity:
        # Weirdly, a scalar Quantity is iterable but you can't call max on it because it's not iterable
        if times.isscalar:
            return([times.value]*times.unit) # make it an array of one (singleton)
        else:
            return (times)
    elif not isinstance(times, collections.abc.Iterable):
        return([times]*u.second) # make it an array of one (singleton)
    else:
        return(times*u.second)


# Make a time series (ephemeris table) from the propagated ephemeris
def ephts(prop):
    ts = apts.TimeSeries(time=[pvt.pvtq.time for pvt in prop.pvt],
                    data={'position': [pvt.pvtq['p'] for pvt in prop.pvt],
                          'velocity': [pvt.pvtq['v'] for pvt in prop.pvt]})
    ts['position'].info.format = '9.3f'
    ts['velocity'].info.format = '9.6f'
    return ts

# Make a list or single PVTs from an ephemeris table
# possibly add later optional mintime, maxtime arguments
def pvts(ephem):
    if type(ephem) == aptbl.row.Row:
        return orbit.PVT(ephem)
    else:
        return [orbit.PVT(row) for row in list(ephem.iterrows())]

# Return a list of numpy arrays and datetimes from the
# ephemeris table or a row of it.
def makenp(ts):
    return (np.concatenate((ts['position'].value, ts['velocity'].value), axis=1), ts['time'].datetime64)
# .makenp() convert to numpy; units are same as ephemeris table but not specified in the result
apts.TimeSeries.makenp = makenp
aptbl.row.Row.makenp = makenp

# Convert a Quantity to seconds as a Python float
def timesec(t):
    if type(t) is u.Quantity and u.get_physical_type(t) == 'time':
        pt = float(t.si.value) # convert to seconds and get the value_unit
    else:
        pt = float(t) # assume seconds
    return(pt)

## Time series of orbital elements
def tselements(ephem, elements):
    return apts.TimeSeries(time=[pvt.pvtq.time for pvt in ephem.pvt],
                      data=[dict(zip(elements,
                                     orbit.elementval(orb, elements,
                                                ephem.forceenv["earthrad"])))
                            for orb in ephem.orbit])
