##### Propagate orbits
## Main function: prop()

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
        newname = "prop" #f"prop{int(pt)}s"
        tree[newname] = orbit.pop('prop')
        return(tree[newname])
    else:
        return(ephgen)

# Propagate the BoundedPropagator
# orbit: tree with orbit in it
# reltimes: times (seconds) past the initial time; may be a number (seconds), Quantity with type 'time' or list of these
#
# The set of propagated points can be added to at any time, but be sure to set maximum_tof
# to the highest possible value when prop() is first called; otherwise, previously computed
# values might be lost.
# The class variable .pvt will have a list of pvts
# The class variable .orbit will have a list of orbits
# The class variable .ephem will have an AstroPy time series ephemeris table
# Pure numpy (no AstroPy or Orekit) is obtained from the ephemeris table with .makenp()
#
# Example
# In [2]: prop(ex1,300.0,maximum_tof=86400.0)
# Out[2]: <PVT position: [4581.814910467976, 4512.262287102674, 1616.8263139353987] (km) velocity:[-4.8914488794011595, 3.141592098656343, 5.166423005465441] (km/s) epoch 2022-06-01 12:05:00 (UTC)>

# In [3]: prop(ex1,[1200.0,1500.0,1800.0])
# Out[3]:
# [<PVT position: [-1359.9602212056286, 4580.209288805089, 4646.557709357471] (km) velocity:[-7.074061594952203, -2.989070151590076, 0.948420112505302] (km/s) epoch 2022-06-01 12:20:00 (UTC)>,
#  <PVT position: [-3358.332160090897, 3427.3199510360746, 4647.312224539923] (km) velocity:[-6.115173953646672, -4.617481772761261, -0.9412695398995137] (km/s) epoch 2022-06-01 12:25:00 (UTC)>,
#  <PVT position: [-4956.791809682075, 1865.869024125194, 4094.2858798679345] (km) velocity:[-4.4362841450361055, -5.686777789101822, -2.7067519585212185] (km/s) epoch 2022-06-01 12:30:00 (UTC)>]

# In [4]: prop(ex1,600.0)
# Out[4]: <PVT position: [2865.4998895513913, 5161.045684511684, 3036.84672810346] (km) velocity:[-6.432386509240479, 1.1404155012024495, 4.2038220495311265] (km/s) epoch 2022-06-01 12:10:00 (UTC)>

# In [5]: ex1.prop.pvt
# Out[5]:
# [<PVT position: [4581.814910467976, 4512.262287102674, 1616.8263139353987] (km) velocity:[-4.8914488794011595, 3.141592098656343, 5.166423005465441] (km/s) epoch 2022-06-01 12:05:00 (UTC)>,
#  <PVT position: [2865.4998895513913, 5161.045684511684, 3036.84672810346] (km) velocity:[-6.432386509240479, 1.1404155012024495, 4.2038220495311265] (km/s) epoch 2022-06-01 12:10:00 (UTC)>,
#  <PVT position: [-1359.9602212056286, 4580.209288805089, 4646.557709357471] (km) velocity:[-7.074061594952203, -2.989070151590076, 0.948420112505302] (km/s) epoch 2022-06-01 12:20:00 (UTC)>,
#  <PVT position: [-3358.332160090897, 3427.3199510360746, 4647.312224539923] (km) velocity:[-6.115173953646672, -4.617481772761261, -0.9412695398995137] (km/s) epoch 2022-06-01 12:25:00 (UTC)>,
#  <PVT position: [-4956.791809682075, 1865.869024125194, 4094.2858798679345] (km) velocity:[-4.4362841450361055, -5.686777789101822, -2.7067519585212185] (km/s) epoch 2022-06-01 12:30:00 (UTC)>]

# In [6]: ex1.prop.orbit
# Out[6]:
# [<Orbit: Cartesian parameters: {P(4581814.910467976, 4512262.287102674, 1616826.3139353986), V(-4891.4488794011595, 3141.592098656343, 5166.423005465442)}>,
#  <Orbit: Cartesian parameters: {P(2865499.8895513914, 5161045.684511684, 3036846.72810346), V(-6432.386509240479, 1140.4155012024494, 4203.822049531126)}>,
#  <Orbit: Cartesian parameters: {P(-1359960.2212056285, 4580209.288805089, 4646557.709357471), V(-7074.061594952203, -2989.070151590076, 948.420112505302)}>,
#  <Orbit: Cartesian parameters: {P(-3358332.160090897, 3427319.9510360747, 4647312.224539923), V(-6115.173953646672, -4617.481772761261, -941.2695398995137)}>,
#  <Orbit: Cartesian parameters: {P(-4956791.809682075, 1865869.024125194, 4094285.879867934), V(-4436.284145036106, -5686.777789101822, -2706.7519585212185)}>]

# In [8]: ex1.prop.ephem
# Out[8]:
# <TimeSeries length=5>
#         time               position               velocity
#                               km                   km / s
#         Time              float64[3]             float64[3]
# ------------------- ---------------------- ----------------------
# 2022-06-01 12:05:00  4581.815 ..  1616.826 -4.891449 ..  5.166423
# 2022-06-01 12:10:00  2865.500 ..  3036.847 -6.432387 ..  4.203822
# 2022-06-01 12:20:00 -1359.960 ..  4646.558 -7.074062 ..  0.948420
# 2022-06-01 12:25:00 -3358.332 ..  4647.312 -6.115174 .. -0.941270
# 2022-06-01 12:30:00 -4956.792 ..  4094.286 -4.436284 .. -2.706752

# In [8]: ex1.prop.ephem.makenp()
# Out[8]:
# [array([[ 4581.81491047,  4512.2622871 ,  1616.82631394],
#         [ 2865.49988955,  5161.04568451,  3036.8467281 ],
#         [-1359.96022121,  4580.20928881,  4646.55770936],
#         [-3358.33216009,  3427.31995104,  4647.31222454],
#         [-4956.79180968,  1865.86902413,  4094.28587987]]),
#  array([[-4.89144888,  3.1415921 ,  5.16642301],
#         [-6.43238651,  1.1404155 ,  4.20382205],
#         [-7.07406159, -2.98907015,  0.94842011],
#         [-6.11517395, -4.61748177, -0.94126954],
#         [-4.43628415, -5.68677779, -2.70675196]]),
#  array(['2022-06-01T12:05:00.000000000', '2022-06-01T12:10:00.000000000',
#         '2022-06-01T12:20:00.000000000', '2022-06-01T12:25:00.000000000',
#         '2022-06-01T12:30:00.000000000'], dtype='datetime64[ns]')]

# In [29]: pvts(ex1.prop.ephem) # same as ex1.prop.pvt
# Out[29]:
# [<PVT position: [4581.814910467976, 4512.262287102674, 1616.8263139353987] (km) velocity:[-4.8914488794011595, 3.141592098656343, 5.166423005465441] (km / s) epoch 2022-06-01 12:05:00 (UTC)>,
#  <PVT position: [2865.4998895513913, 5161.045684511684, 3036.84672810346] (km) velocity:[-6.432386509240479, 1.1404155012024495, 4.2038220495311265] (km / s) epoch 2022-06-01 12:10:00 (UTC)>,
#  <PVT position: [-1359.9602212056286, 4580.209288805089, 4646.557709357471] (km) velocity:[-7.074061594952203, -2.989070151590076, 0.948420112505302] (km / s) epoch 2022-06-01 12:20:00 (UTC)>,
#  <PVT position: [-3358.332160090897, 3427.3199510360746, 4647.312224539923] (km) velocity:[-6.115173953646672, -4.617481772761261, -0.9412695398995137] (km / s) epoch 2022-06-01 12:25:00 (UTC)>,
#  <PVT position: [-4956.791809682075, 1865.869024125194, 4094.2858798679345] (km) velocity:[-4.4362841450361055, -5.686777789101822, -2.7067519585212185] (km / s) epoch 2022-06-01 12:30:00 (UTC)>]

def prop(orbit, reltimes, stopalt=125.0e3, maximum_tof=0.0):
    if isinstance(reltimes, collections.abc.Iterable):
        maxtime = max(max(reltimes), maximum_tof)
    else:
        maxtime = max(reltimes, maximum_tof)
    [bp, tree] = ecis.thingofclass(orbit, BoundedPropagator)
    if tree is None:
        forceenv = None
    else:
        forceenv = tree.forceenv
    if bp is None or ("maxtime" in tree and type(tree.maxtime) is float and maxtime > tree.maxtime):
        bp = mkephem(orbit, maxtime, forceenv, stopalt)
        [bp, tree] = ecis.thingofclass(orbit, BoundedPropagator)
    if isinstance(reltimes, collections.abc.Iterable):
        state = [bp.propagate(bp.getMinDate().shiftedBy(timesec(rt))).orbit for rt in reltimes]
        pvt=[PVT(st) for st in state]
    else:
        state = bp.propagate(bp.getMinDate().shiftedBy(timesec(reltimes))).orbit
        pvt=PVT(state)
    pvts = ensurelist(tree.get('pvt')) + ensurelist(pvt)
    pvts.sort(key = lambda s: s.pvtq.time)
    tree.update(pvt = pvts)
    states = ensurelist(tree.get('orbit')) + ensurelist(state)
    states.sort(key = lambda s: PVT(s).pvtq.time)
    tree.update(orbit=states)
    tree.update(ephem=ephts(tree))
    return(pvt)

# Make a time series (ephemeris table) from the propagated ephemeris
def ephts(prop):
    ts = TimeSeries(time=[pvt.pvtq.time for pvt in prop.pvt],
                    data={'position': [pvt.pvtq['p'] for pvt in prop.pvt],
                          'velocity': [pvt.pvtq['v'] for pvt in prop.pvt]})
    ts['position'].info.format = '9.3f'
    ts['velocity'].info.format = '9.6f'
    return ts

# Make a list or single PVTs from an ephemeris table
# possibly add later optional mintime, maxtime arguments
def pvts(ephem):
    if type(ephem) == astropy.table.row.Row:
        return PVT(ephem)
    else:
        return [PVT(row) for row in list(ephem.iterrows())]

# Return a list of numpy arrays and datetimes from the
# ephemeris table or a row of it.
def makenp(ts):
    return (np.concatenate((ts['position'].value, ts['velocity'].value), axis=1), ts['time'].datetime64)
# .makenp() convert to numpy; units are same as ephemeris table but not specified in the result
TimeSeries.makenp = makenp
astropy.table.row.Row.makenp = makenp

# Convert a Quantity to seconds as a Python float
def timesec(t):
    if type(t) is Quantity and get_physical_type(t) == 'time':
        pt = float(t.si.value) # convert to seconds and get the value_unit
    else:
        pt = float(t) # assume seconds
    return(pt)
