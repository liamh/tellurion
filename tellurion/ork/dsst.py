from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import PropagationType, SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.propagation.semianalytical.dsst import DSSTPropagator
from org.orekit.propagation.semianalytical.dsst.forces import (
    DSSTAtmosphericDrag,
    DSSTSolarRadiationPressure,
    DSSTTesseral,
    DSSTZonal,
)
from org.orekit.utils import Constants

from tellurion.ork import convert
from tellurion.ork.prop import _additional, _convert_to_orbit, _make_generator


def dsstprep(initstate, proptime, events, forceenv, reftime, output, state_type="mean"):
    """Prepare the DSST propagator"""
    ork0 = _convert_to_orbit(initstate, forceenv)

    # DSST integration parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0e-3
    tolerances = NumericalPropagator.tolerances(positionTolerance, ork0, ork0.getType())

    # Create initial spacecraft state
    initialState = SpacecraftState(ork0, forceenv["mass"])

    def make_propagator(propagation_type):
        integrator = DormandPrince853Integrator(
            minstep, maxstep, tolerances[0], tolerances[1]
        )
        integrator.setInitialStepSize(initStep)
        propagator = DSSTPropagator(integrator, propagation_type)
        propagator.setInitialState(initialState)

        degree, order = forceenv.get("gravity-degree-order", (0, 0))
        if degree >= 2:
            propagator.addForceModel(
                DSSTZonal(forceenv["earthframe"], forceenv["gravity-unnorm"])
            )
            if order > 0:
                propagator.addForceModel(
                    DSSTTesseral(
                        forceenv["earthframe"],
                        Constants.WGS84_EARTH_ANGULAR_VELOCITY,
                        forceenv["gravity-unnorm"],
                    )
                )

        if "dragforce" in forceenv:
            propagator.addForceModel(
                DSSTAtmosphericDrag(forceenv["dragforce"], forceenv["earthmu"].si.value)
            )

        if "srpmodel" in forceenv:
            srp = DSSTSolarRadiationPressure(
                forceenv["srpmodel"], forceenv["crp"], forceenv["area"], forceenv["mass"]
            )
            propagator.addForceModel(srp)
        return propagator

    # Event detectors must see osculating states, independently of output.
    propagator = make_propagator(PropagationType.OSCULATING)

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

    if state_type == "mean":
        mean_propagator = make_propagator(PropagationType.MEAN)
        mean_ephgen = mean_propagator.getEphemerisGenerator()
        final_date = generator["final"]["state"].getDate()
        mean_propagator.propagate(generator["epoch"], final_date)
        mean_eph = mean_ephgen.getGeneratedEphemeris()
        generator["osculating propfn"] = generator["propfn"]
        generator["propfn"] = lambda date: mean_eph.propagate(date)
        generator["pvt0"] = convert._pvt(generator["propfn"](generator["epoch"]))
        generator["final"]["state"] = generator["propfn"](final_date)
        generator["final"]["pvt"] = convert._pvt(generator["final"]["state"])

    return generator
