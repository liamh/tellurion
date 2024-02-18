from orbit import *

# Propagation and ephemeris
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from orekit import JArray_double
from org.orekit.forces.drag import AbstractDragForceModel, DragForce
from org.orekit.models.earth.atmosphere import Atmosphere, HarrisPriester, DTM2000, NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.forces.drag import IsotropicDrag

# Orbital environment constants
envct['earth'] = OneAxisEllipsoid(envct['earthrad'], envct['earthflat'],  envct['earthframe'])
envct['sphearth'] = OneAxisEllipsoid(envct['earthrad'], 0.0,  envct['earthframe'])
envct['swdata'] = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")
envct['stopalt'] = 125.0e3  # Altitude at which propagation should stop

# Atmospheric density models
envct['hp'] = HarrisPriester(envct['sun'], envct['earth']) # Harris-Priester atmospheric density model
envct['dtm'] = DTM2000(envct['swdata'], envct['sun'], envct['earth']) # DTM2000 atmospheric density model
envct['msis'] = NRLMSISE00(envct['swdata'], envct['sun'], envct['earth'])

# Create an example spacecraft with B = C_D A/m = 0.01 m^2/kg
scB010 = {'mass': 100.0,  # The models need a spacecraft mass, unit kg.
        'dragarea': 1.0, # Cross-sectional area perpendicular to atmosphere direction, m^2
        'dragcoef': 1.0 # Coefficient of drag
        }
scB010['drag'] = IsotropicDrag(scB010['dragarea'], scB010['dragcoef'])
scB010['atmdens'] = envct['hp']
scB010['dragforce'] = DragForce(scB010['atmdens'], scB010['drag']);

def atmdens(location, time, model = 'hp'):
    return(envct[model].getDensity(time, location, envct['celestialframe']))

def setgravity(degree, order):
    envct['gravity-degree-order'] = [degree, order]
    envct['gravity'] = GravityFieldFactory.getNormalizedProvider(degree, order)

# Set the default gravity 0x0
setgravity(0, 0)

# Propagate from epoch for a specified time
def prop(orbit,proptime,spacecraft=scB010):
    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, orbit, orbit.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
    	       minstep,
    	       maxstep,
	       JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
	       JArray_double.cast_(tolerances[1]))
    integrator.setInitialStepSize(initStep)

    # Initialize the spacecraft state
    # satellite_mass = 100.0  # The models need a spacecraft mass, unit kg.
    initialState = SpacecraftState(orbit, spacecraft['mass'])
    prop = NumericalPropagator(integrator)
    prop.setOrbitType(okc['cartesian'])
    prop.setInitialState(initialState)
    generator = prop.getEphemerisGenerator()

    # Forces
    prop.addForceModel(HolmesFeatherstoneAttractionModel(envct['earthframe'], envct['gravity']))
    prop.addForceModel(spacecraft['dragforce'])

    # Events
    altdet = AltitudeDetector(envct['stopalt'], envct['sphearth'])
    prop.addEventDetector(altdet)

    # Propagate
    propagated = prop.propagate(orbit.date, orbit.date.shiftedBy(proptime))
    ephemeris = generator.getGeneratedEphemeris();
    return(ephemeris)

# Lookup the state at a particular time that is with the range bounded by the minimum and maximum times
# eph: output from prop()
# reltime: time (seconds) past the earliest time of the propagation
def ephlookup(eph, reltime):
    return(eph.propagate(eph.getMinDate().shiftedBy(reltime)).orbit)

# Example, see example in orbit.pv
# ex1eph1hr = prop(ex1orb,3600.0)
# ex1prop1200s = ephlookup(ex1eph1hr,1200.0)
# pvt(ex1prop1200s)
# Out[20]: <TimeStampedPVCoordinates: {2022-06-01T12:20:00.000, P(-1363975.7207431477, 4575863.159211461, 4639461.68560796), V(-7076.625054306512, -2994.509027288756, 933.643196058037), A(1.8386249331276936, -6.167921710610659, -6.272353908960307)}>
# pvtpy(pvt(ex1prop1200s))
# [-1359953.707299648, 4580203.17356122, 4646548.529260602],
#  [-7074.051859352501, -2989.0818948108367, 948.4044990087592],
#  datetime.datetime(2022, 6, 1, 12, 20)]
# convert(ex1prop1200s,'kep')
# Out[21]: <Orbit: Keplerian parameters: {a: 6662696.447232186; e: 0.005515508998697014; i: 44.958136752610166; pa: -1.8475032055235747; raan: 29.924998771236737; v: 82.3318582445891;}>
# posmag(ex1prop1200s)
# Out[22]: 6657594.021178353
