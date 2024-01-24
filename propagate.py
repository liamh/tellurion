from setup import *

# Propagate from epoch for a specified time
def prop(orbit,proptime):
    # Set parameters
    minStep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, orbit, orbit.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
    	       minStep,
    	       maxstep,
	       JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
	       JArray_double.cast_(tolerances[1]))
    integrator.setInitialStepSize(initStep)

    # Initialize the spacecraft state
    satellite_mass = 100.0  # The models need a spacecraft mass, unit kg.
    initialState = SpacecraftState(orbit, satellite_mass)
    prop = NumericalPropagator(integrator)
    prop.setOrbitType(OrbitType.CARTESIAN)
    prop.setInitialState(initialState)
    generator = prop.getEphemerisGenerator()

    # Forces
    gravityProvider = GravityFieldFactory.getNormalizedProvider(10, 10) # 10x10
    prop.addForceModel(HolmesFeatherstoneAttractionModel(itrf, gravityProvider))

    # Propagate
    propagated = prop.propagate(orbit.date, orbit.date.shiftedBy(proptime))
    ephemeris = generator.getGeneratedEphemeris();
    return(ephemeris)

def cartorb(posv3d, velv3d, datetime):
    pvt = TimeStampedPVCoordinates(datetime, posv3d, velv3d) # Make the PVT initial state
    return(CartesianOrbit(pvt, FramesFactory.getEME2000(), Constants.WGS84_EARTH_MU))

# Create the PVT orbit
def orbitpvt(pos, vel, datetime):
    cartorb(Vector3D(pos), Vector3D(vel), datetime)

# The position-velocity-time for the state
def orbpvt(orbit):
    return(orbit.pVCoordinates)

# The geocentric distance of the orbit
def posmag(orbit):
    return(orbit.pVCoordinates.position.norm)

# Lookup the state at a particular time that is with the range bounded by the minimum and maximum times
# eph: output from prop()
# reltime: time (seconds) past the earliest time of the propagation
def ephlookup(eph, reltime):
    return(eph.propagate(eph.getMinDate().shiftedBy(reltime)).orbit)

# ex1p = [5740132.68349499, 3314067.15, 0.0]
# ex1v = [-2750.82683526322, 4764.5718414998, 5501.65367052644]
# ex1t = AbsoluteDate(2022, 6, 1, 12, 0, 00.000, utc)
# ex1pvt = orbitpvt(ex1p,ex1v,ex1t)
# ex1eph1hr = prop(ex1pvt,3600.0)
# ex1prop1200s = ephlookup(ex1eph1hr,1200.0)
# orbpvt(ex1prop1200s)
# Out[20]: <TimeStampedPVCoordinates: {2022-06-01T12:20:00.000, P(-1363975.7207431477, 4575863.159211461, 4639461.68560796), V(-7076.625054306512, -2994.509027288756, 933.643196058037), A(1.8386249331276936, -6.167921710610659, -6.272353908960307)}>
# orbkep(ex1prop1200s)
# Out[21]: <Orbit: Keplerian parameters: {a: 6662696.447232186; e: 0.005515508998697014; i: 44.958136752610166; pa: -1.8475032055235747; raan: 29.924998771236737; v: 82.3318582445891;}>
# posmag(ex1prop1200s)
# Out[22]: 6657594.021178353
