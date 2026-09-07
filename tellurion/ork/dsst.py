from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.propagation.semianalytical.dsst import DSSTPropagator
from org.orekit.propagation.semianalytical.dsst.forces import (
    DSSTAtmosphericDrag,
    DSSTSolarRadiationPressure,
    DSSTTesseral,
    DSSTZonal,
)
from org.orekit.utils import Constants

from tellurion.ork.prop import _additional, _convert_to_orbit, _make_generator


def dsstprep(initstate, proptime, events, forceenv, reftime, output):
    """Prepare the DSST propagator"""
    ork0 = _convert_to_orbit(initstate, forceenv)

    # DSST integration parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0e-3
    tolerances = NumericalPropagator.tolerances(positionTolerance, ork0, ork0.getType())

    # Initialize the integrator for DSST's numerical part
    integrator = DormandPrince853Integrator(
        minstep, maxstep, tolerances[0], tolerances[1]
    )
    integrator.setInitialStepSize(initStep)

    # Create initial spacecraft state
    initialState = SpacecraftState(ork0, forceenv["mass"])

    # Create DSST propagator
    propagator = DSSTPropagator(integrator)
    propagator.setInitialState(initialState)

    # Add force models
    if forceenv.get("gravity-degree-order"):
        _, order = forceenv["gravity-degree-order"]
        zonal = DSSTZonal(forceenv["gravity-unnorm"])
        propagator.addForceModel(zonal)

        if order > 0:
            tesseral = DSSTTesseral(
                forceenv["earthframe"],
                Constants.WGS84_EARTH_ANGULAR_VELOCITY,
                forceenv["gravity-unnorm"],
            )
            propagator.addForceModel(tesseral)

    if "dragforce" in forceenv:
        atm_drag = DSSTAtmosphericDrag(
            forceenv["dragmodel"], forceenv["cdrag"], forceenv["area"], forceenv["mass"]
        )
        propagator.addForceModel(atm_drag)

    if "srpmodel" in forceenv:
        srp = DSSTSolarRadiationPressure(
            forceenv["srpmodel"], forceenv["crp"], forceenv["area"], forceenv["mass"]
        )
        propagator.addForceModel(srp)

    # Get ephemeris generator BEFORE creating the main generator
    ephgen = propagator.getEphemerisGenerator()

    # Create generator with lambda using generator["epoch"] (will be
    # populated by _make_generator)
    generator = _make_generator(
        ork0, lambda propto: propagator.propagate(generator["epoch"], propto)
    )

    # Perform initial propagation and additional calculations
    _additional(events, propagator, generator, proptime, forceenv, reftime, output)

    # NOW get the generated ephemeris (after _additional has run)
    gge = ephgen.getGeneratedEphemeris()

    # Replace propfn to use interpolation for subsequent queries
    generator["propfn"] = lambda propto: gge.propagate(propto)
    generator["mindate"] = gge.getMinDate()
    generator["maxdate"] = gge.getMaxDate()

    return generator
