"""Multibody propagation in the inertial solar-system-barycentric frame.

The spacecraft state is propagated as ``AbsolutePVCoordinates`` in ICRF
(origin at the solar-system barycenter). Every attracting body contributes
one ``SingleBodyAbsoluteAttraction``; no central attraction is added, and
because the frame is inertial no inertial-force term is needed. This is
independent of `tellurion.prepare` / `tellurion.propagate`, which remain
the Earth-orbit propagators.
"""

import collections.abc

import astropy.units as u

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.bodies import CelestialBodyFactory
from org.orekit.forces.gravity import SingleBodyAbsoluteAttraction
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.utils import AbsolutePVCoordinates, TimeStampedPVCoordinates

import tellurion.astro.time as atime
from tellurion.core.frame import FrameDescriptor
from tellurion.ork import convert
from tellurion.ork.frame import resolve_orekit_frame, transform_pvt

SSB_ICRF = FrameDescriptor("Solar System Barycenter", "ICRF")
DEFAULT_BODIES = ("Sun", "Earth", "Moon")


def _body(name):
    if not isinstance(name, str):
        return name  # assume an Orekit CelestialBody
    key = name.strip()
    for candidate in (key, key.title(), key.capitalize()):
        try:
            return CelestialBodyFactory.getBody(candidate)
        except Exception:
            continue
    raise ValueError(f"Unknown celestial body {name!r}")


def mbprepare(pvt, duration, bodies=DEFAULT_BODIES, rtol=1e-12, atol=1e-3,
              minstep=1e-3, maxstep=86400.0):
    """Prepare a multibody propagation in the barycentric inertial frame.

    Parameters
    ----------
    pvt : `~tellurion.PositionVelocityTime`
        Initial state; must be a scalar with an epoch and a ``frame``
        descriptor (e.g. Earth/EME2000 or Solar System Barycenter/ICRF).
        It is converted to barycentric ICRF internally.
    duration : `~astropy.units.Quantity`
        Time span to cover from the epoch; negative values propagate
        backward.
    bodies : sequence of str, optional
        Attracting bodies by Orekit ``CelestialBodyFactory`` name, e.g.
        ``["Sun", "Earth", "Moon", "Jupiter"]``.  Default Sun, Earth, Moon.
    rtol, atol : float, optional
        Relative and absolute tolerances of the Dormand-Prince 8(5,3)
        integrator (SI units).
    minstep, maxstep : float, optional
        Integrator step limits in seconds.

    Returns
    -------
    dict
        Generator for `mbpropagate`.
    """
    if not pvt.isscalar:
        raise ValueError("pvt must be a single state")
    if pvt.frame is None:
        raise ValueError("pvt must have a frame descriptor")
    if isinstance(bodies, str):
        bodies = [bodies]
    if not len(bodies):
        raise ValueError("At least one attracting body is required")

    icrf = resolve_orekit_frame(SSB_ICRF)
    bary = transform_pvt(pvt, SSB_ICRF)
    epoch = convert._abstime_to_okad(bary.time)
    span = u.Quantity(duration).to(u.s).value

    state = SpacecraftState(AbsolutePVCoordinates(icrf, convert._tspvc(bary)))
    integrator = DormandPrince853Integrator(
        float(minstep), float(maxstep), float(atol), float(rtol)
    )
    propagator = NumericalPropagator(integrator)
    propagator.setOrbitType(None)  # propagate absolute PV, not orbital elements
    propagator.setIgnoreCentralAttraction(True)
    names = []
    for body in bodies:
        cb = _body(body)
        names.append(str(cb.getName()))
        propagator.addForceModel(SingleBodyAbsoluteAttraction(cb))
    propagator.setInitialState(state)

    ephgen = propagator.getEphemerisGenerator()
    propagator.propagate(epoch.shiftedBy(float(span)))
    return {
        "epoch": epoch,
        "mindate": epoch if span >= 0 else epoch.shiftedBy(float(span)),
        "maxdate": epoch.shiftedBy(float(span)) if span >= 0 else epoch,
        "frame": SSB_ICRF,
        "bodies": names,
        "icrf": icrf,
        "ephemeris": ephgen.getGeneratedEphemeris(),
        "init": bary,
    }


def mbpropagate(generator, reltimes, include_init=False, frame=None):
    """Evaluate a `mbprepare` generator at the requested times.

    Parameters
    ----------
    generator : dict
        From `mbprepare`.
    reltimes : `~astropy.units.Quantity` or `~astropy.time.Time`
        Times relative to the epoch, or absolute times.
    include_init : bool, optional
        Also return the initial state.
    frame : `~tellurion.FrameDescriptor`, optional
        Frame of the result; default barycentric ICRF.  Earth-centered
        EME2000 is ``FrameDescriptor("Earth", "EME2000")``.

    Returns
    -------
    `~tellurion.PositionVelocityTime`
        A time series (or scalar) tagged with its frame.
    """
    reft = convert._abstime_from_okad(generator["epoch"])
    atimes = atime.abstime(reltimes, reft)
    if include_init:
        atimes = atime.time_concat(reft, atimes)
    icrf = generator["icrf"]
    eph = generator["ephemeris"]

    def one(t):
        date = convert._abstime_to_okad(t)
        apv = eph.propagate(date).getAbsPVA()
        return TimeStampedPVCoordinates(
            date, apv.getPosition(date, icrf), apv.getVelocity(date, icrf)
        )

    if atimes.isscalar:
        result = convert._pvt(one(atimes))
    else:
        result = convert._pvt([one(t) for t in atimes])
    result.frame = generator["frame"]
    if frame is not None:
        result = transform_pvt(result, frame)
    return result
