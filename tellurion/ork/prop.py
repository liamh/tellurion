"""
Propagate orbit states in two steps: `prepare()` and `propagate()`.

Overview of Propagators
-----------------------
1. **Kepler Analytic**: Specified with ``kepleranalytic()`` when initial
   state is a ``PositionVelocityT``.
2. **Brouwer-Lyddane (J2 perturbations)**: Specified with ``setgravity()``
   when initial state is a ``PositionVelocityT``.
3. **DSST (semi-analytical)**: Specified with ``setgravity()`` and optional
   ``dragforce()`` / ``srpmodel()`` when initial state is a
   ``PositionVelocityT``. Ideal for long-term propagation.
4. **Numerical Integration**: Specified by ``setgravity()`` with optional
   atmospheric drag using ``dragforce()`` when initial state is a
   ``PositionVelocityT``.
5. **SGP4 Mean Element Propagation**: Automatically used if the initial
   state is a ``MeanElementSetT``.

Examples
--------
.. code-block:: python

    import astropy.units as u
    import numpy as np

    isssent = tell.spacetrack_latest(stclient, [25544, 41335])
    sentprep = tork.SGP4prep(
        isssent["SENTINEL 3A"],
        1 * u.day,
        {
            "altitude": 125.0 * u.km,
            "eclipse": [True, True],
            "visibility": [],
        },
    )
    tork.propagate(sentprep, np.linspace(5.0 * u.minute, 60.0 * u.minute, 12), True)

Notes
-----
**Supported Events**

Events are configured in the ``events`` dictionary using the following keys:

- ``'altitude'``: Terminate propagation when altitude drops below threshold.
- ``'eclipse'``: Track eclipsing events.
- ``'visibility'``: Calculate visibility from earth locations.

Setting ``'stm'`` to ``True`` in ``events`` requests the final state
transition matrix, which will be stored under ``gen['final']['stm']``.

**Custom Orekit Events**

Adding new event detectors requires defining three support functions:
``_mkdetlog()``, ``_gentrans()``, and ``_statechar()``.
"""

import astropy.units as u
import org.orekit.forces.gravity as okgrav
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical import (
    BrouwerLyddanePropagator,
    KeplerianPropagator,
)
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.utils import (
    AbsolutePVCoordinates,
    TimeStampedPVCoordinates,
)

import tellurion.astro.time as atime
from tellurion.core import MeanElementSetT, PositionVelocityT, element
from tellurion.ork.element import (
    _circularorbit,
    _equinoctialorbit,
    _keporb_from_components,
)

from . import convert, force, jacobian
from .event import prop as event

defev = {"altitude": 125.0 * u.km, "eclipse": [], "visibility": [], "stm": False}


def prepare(
    initstate,
    proptime,
    events=defev,
    forceenv=force.deffe,
    reftime="epoch",
    output="et",
    propagator="auto",
):
    """
    Prepare a propagation

    Parameters
    ----------
    initstate  : initial state, must be in a form acceptable to the propagator
                 chosen
    proptime   : time of propagation
    events     : dict
    forceenv   : dict
    reftime    : reference time used to fill the relevant column of the ephemeris table
    output     : string
        'et' - ephemeris table
        'ss' - Orekit SpacecraftState
    propagator : str
        'auto' - automatically select based on initstate and forceenv
        'keplerian' - two-body analytic
        'brouwer-lyddane' - J2 perturbations analytic
        'dsst' - semi-analytical with multiple perturbations
        'numerical' - numerical integration
        'sgp4' - SGP4 mean elements (requires MeanElementSetT)
    """
    if propagator == "auto":
        propagator = _select_propagator(initstate, forceenv)

    propagator_funcs = {
        "sgp4": SGP4prep,
        "keplerian": kaprep,
        "brouwer-lyddane": blprep,
        "numerical": niprep,
    }

    # Lazy import for dsst to avoid circular dependency
    if propagator == "dsst":
        from .dsst import dsstprep

        propagator_funcs["dsst"] = dsstprep

    propagator_func = propagator_funcs.get(propagator)

    if not propagator_func:
        raise ValueError(f"Unknown propagator: {propagator}")

    # ... rest of validation code
    return propagator_func(initstate, proptime, events, forceenv, reftime, output)


def _select_propagator(initstate, forceenv):
    """Auto-select propagator based on state type and force environment"""
    if isinstance(initstate, MeanElementSetT):
        return "sgp4"
    elif forceenv.get("gravity-degree-order"):
        return "numerical"
    else:
        return "keplerian"


def SGP4prep(meanels, proptime, events, forceenv, reftime, output):
    """Propagate mean elements using SGP4"""
    if hasattr(meanels, "model") and meanels.model == "SGP4":
        tle = TLE(*meanels.tle)
        propagator = TLEPropagator.selectExtrapolator(tle)
        generator = _make_generator(
            tle,
            lambda propto: propagator.getPVCoordinates(
                propto, forceenv["celestialframe"]
            ),
        )
        generator["pvt0"] = convert._pvt(generator["propfn"](generator["epoch"]))
        _additional(events, propagator, generator, proptime, forceenv, reftime, output)
        return generator
    raise ValueError("Can only propagate SGP4 mean elements with SGP4")


# -----------------------------------------------------------------------
# Individual propagators
# -----------------------------------------------------------------------


# All propagators except SGP4 can take ElementSetT or
# PositionVelocityT initial state
def _convert_to_orbit(state, forceenv):
    """Convert the ElementSetT or PositionVelocityT state to an Orekit Orbit"""
    if element.iskepels(state):
        conv = _keporb_from_components(state.elements, state.time)
    elif element.isequels(state):
        conv = _equinoctialorbit(state, forceenv=forceenv)
    elif element.iscircels(state):
        conv = _circularorbit(state, forceenv=forceenv)
    elif isinstance(state, PositionVelocityT):
        conv = CartesianOrbit(
            convert._tspvc(state),
            forceenv["celestialframe"],
            forceenv["earthmu"].si.value,
        )
    else:
        raise TypeError(
            "Can only convert Keplerian, Equinoctial, or Circular "
            "ElementSetT, or PositionVelocityT objects"
        )
    return conv


def kaprep(initstate, proptime, events, forceenv, reftime, output):
    """Prepare the Keplerian (two-body) analytic propagator"""
    ork0 = _convert_to_orbit(initstate, forceenv)
    propagator = KeplerianPropagator(ork0, forceenv["earthmu"].si.value)
    generator = _make_generator(ork0, lambda propto: propagator.propagate(propto))
    generator["pvt0"] = convert._pvt(generator["propfn"](generator["epoch"]))
    _additional(events, propagator, generator, proptime, forceenv, reftime, output)
    return generator


# Prepare Brouwer-Lyddane with defaults
def blprep_def(initstate, proptime, events, forceenv, reftime, output):
    """Prepare the Brouwer-Lyddane analytic propagator (J2 perturbations)"""
    ork0 = _convert_to_orbit(initstate, forceenv)
    # M2 is along-track acceleration primarily caused by atmospheric
    # drag, like the B* term in SGP4. Here we set it to zero.
    m2 = 0.0
    propagator = BrouwerLyddanePropagator(
        ork0, float(forceenv["earthmu"].si.value), forceenv["gravity-unnorm"], m2
    )
    generator = _make_generator(ork0, lambda propto: propagator.propagate(propto))
    generator["pvt0"] = convert._pvt(generator["propfn"](generator["epoch"]))
    _additional(events, propagator, generator, proptime, forceenv, reftime, output)
    return generator


def blprep(initstate, proptime, events, forceenv, reftime, output):
    """Prepare the Brouwer-Lyddane analytic propagator (J2 perturbations)"""

    ork0 = _convert_to_orbit(initstate, forceenv)
    m2 = 0.0
    epsilon = 1e-5  # Reduced tolerance for near-circular orbit convergence
    max_iterations = 500

    provider = forceenv["gravity-unnorm"]
    harmonics = provider.onDate(ork0.getDate())

    # Use computeMeanOrbit with epsilon to convert osculating to mean orbit
    mean_orbit = BrouwerLyddanePropagator.computeMeanOrbit(
        ork0, provider, harmonics, m2, epsilon, max_iterations
    )

    # Now create the propagator with the mean orbit
    propagator = BrouwerLyddanePropagator(
        mean_orbit,
        1.0,  # default mass
        provider,
        m2,
    )
    generator = _make_generator(mean_orbit, lambda propto: propagator.propagate(propto))
    generator["pvt0"] = convert._pvt(generator["propfn"](generator["epoch"]))
    _additional(events, propagator, generator, proptime, forceenv, reftime, output)
    return generator


def niprep(initstate, proptime, events, forceenv, reftime, output):
    """Make a generator for an ephemeris, optionally include eclipse
    information. The result of this function is passed as the first
    argument to `propagate()`.

    Parameters
    ----------
    initstate: a PositionVelocityT representing the initial state

    proptime:  u.Quantity, float
      The maximum time to propagate, numbers are in seconds

    forceenv:  dict
      Forces to use; output of force.setgravity()

    stopalt:   float
      Stop propagation if altitude above spherical earth (m) drops below
    this threshold

    Returns
    -------
    An Orekit object that is passed to `propagate()` as the first
    argument, if no event other than altitude is included in `events`.

    If eclipse detection is added to `events`, then a dictionary with
    the generator `['ephgen']`, and an ephemeris table `['sun
    transitions']` with xyz positions is returned. The latter has a
    two-character string 'suntrans', and elapsed time from the
    previous transition 'elapsed'. There are also two detectors used
    by `propagate()` in the dictionary.

    """

    ork0 = _convert_to_orbit(initstate, forceenv)

    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0e-3
    tolerances = NumericalPropagator.tolerances(positionTolerance, ork0, ork0.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
        minstep, maxstep, tolerances[0], tolerances[1]
    )
    integrator.setInitialStepSize(initStep)

    initialState = SpacecraftState(ork0, forceenv["mass"])
    propagator = NumericalPropagator(integrator)
    propagator.setResetAtEnd(False)
    propagator.setOrbitType(OrbitType.CARTESIAN)
    propagator.setInitialState(initialState)
    ephgen = propagator.getEphemerisGenerator()

    # Forces
    propagator.addForceModel(
        okgrav.HolmesFeatherstoneAttractionModel(
            forceenv["earthframe"], forceenv["gravity"]
        )
    )
    if "dragforce" in forceenv:
        propagator.addForceModel(forceenv["dragforce"])

    # Make generator, additional calculations, and initial propagation
    generator = _make_generator(
        ork0, lambda propto: propagator.propagate(generator["epoch"], propto)
    )
    _additional(events, propagator, generator, proptime, forceenv, reftime, output)
    gge = ephgen.getGeneratedEphemeris()
    generator["propfn"] = lambda propto: gge.propagate(
        propto
    )  # interpolate in the previous integration result
    generator["mindate"] = gge.getMinDate()
    generator["maxdate"] = gge.getMaxDate()

    return generator


def _make_generator(getdatefrom, propfn):
    generator = {}
    generator["event detectors"] = {}
    generator["epoch"] = getdatefrom.getDate()
    generator["propfn"] = propfn
    return generator


def _additional(events, propagator, generator, proptime, forceenv, reftime, output):
    """Additional calculations when propagating initially"""

    # 1) Add pre-propagation actions
    detlogs = event._add_pre(events, propagator, forceenv)  # Events
    compute_pjac = force.compute_drag_pjac(forceenv)
    if events.get("stm") or compute_pjac:
        harvester = jacobian._add_stm(propagator)  # State-transition matrix

    # 2) Propagate, saving output (SpacecraftState)
    ss = generator["propfn"](generator["epoch"].shiftedBy(atime.timesec(proptime)))
    generator["final"] = {"state": ss, "pvt": convert._pvt(ss)}

    # 3) Add post-propagation actions and save results to `generator`
    event._add_post(detlogs, events, generator, reftime, output)
    if events.get("stm"):
        generator["final"]["stm"] = jacobian.stm(harvester, generator["final"]["state"])
    if compute_pjac:
        generator["final"]["parameters jacobian"] = jacobian.pjac(
            harvester, generator["final"]["state"]
        )


def propagate(generator, reltimes, include_init=True, reftime="epoch", output="et"):
    """From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. The relative times must
    be quantities with physical type `'time'` , and if the size
    prop5m1h.shape[0] > 0, an ephemeris table is returned. If reltimes
    is a single time, then a PVT is returned. If `include`_init is
    true, then include the initial PVT in the ephemeris table.

    Output is specified as one of
      `et`: Ephemeris table (an AstroPy time series; default)
      `pvt`: a PVT
      `ss`: SpacecraftState

    If an eclipse detector has been added, the `'sunlight'` value will
    be one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial
    eclipse, or 's' (full sun).
    """

    # Compute the reference time `reft`, an astropy.time.Time, then
    # create all the absolute times
    if reftime == "epoch" and generator.get("epoch"):
        reft = convert._abstime_from_okad(generator.get("epoch"))
    elif reftime == "epoch" and generator.get("mindate"):
        reft = convert._abstime_from_okad(generator.get("mindate"))
    else:
        reft = atime.abstime(reftime)
    atimes = atime.abstime(reltimes, reft)
    if include_init:
        atimes = atime.time_concat(reft, atimes)

    # Generate a spacecraft state at each absolute time, then
    # optionally create a PositionVelocityT and add event states, and
    # convert to ephemeris table
    if atimes.isscalar:
        ss = _to_spacecraft_state(generator["propfn"](convert._abstime_to_okad(atimes)))
    else:
        ss = [
            _to_spacecraft_state(generator["propfn"](convert._abstime_to_okad(at)))
            for at in atimes
        ]
    if output == "ss":
        return ss
    pvt = convert._pvt(ss)
    pvt.aux = event._evstates(generator, ss)
    if output == "pvt":
        return pvt.timeorder()
    else:
        return pvt.timeorder().ephemeris()


def timerange(object):
    """The time difference between the earliest (usually the initial
    time) and the latest; not always what is requested as atmospheric
    drag can shorten the timespan
    """
    if object.get("maxdate") and object.get("mindate"):
        return (
            convert._abstime_from_okad(object.get("maxdate"))
            - convert._abstime_from_okad(object.get("mindate"))
        ).to(u.s)
    else:
        raise ValueError("Cannot get timerange for this object")


def _to_spacecraft_state(tspvc, forceenv=force.deffe):
    """From the TimeStampedPVCoordinates, create a SpacecraftState"""
    if isinstance(tspvc, TimeStampedPVCoordinates):
        apvc = AbsolutePVCoordinates(forceenv["celestialframe"], tspvc)
        return SpacecraftState(apvc)
    else:
        return tspvc
