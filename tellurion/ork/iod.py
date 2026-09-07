from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.estimation.iod import IodLambert

from tellurion.astro import quantity_utils as qu
from tellurion.core import posvel
from tellurion.ork import convert, force

# Currently does not get long-period solutions
# This should be possible in Orekit 14: https://forum.orekit.org/t/orekit-lambert-solver/4265/16

# pt1 = tell.pvtcart(demoa.propa.cartephem[0]['position'],
# demoa.propa.cartephem[0]['time']) pt2 =
# tell.pvtcart(demoa.propa.cartephem[-1]['position'],
# demoa.propa.cartephem[-1]['time'])


def lambert(
    p1: posvel.PositionT,
    p2: posvel.PositionT,
    shortway: bool = True,
    n_rev: int = 0,
    forceenv=force.deffe,
) -> tuple[posvel.PositionVelocityT, posvel.PositionVelocityT]:
    """Solve Lambert's problem to determine the orbital transfer between two points.

    Uses Orekit's ``IodLambert`` solver to compute the trajectory connecting two
    positions at their respective epochs.

    Parameters
    ----------
    p1 : posvel.PositionT
        Initial position and epoch.
    p2 : posvel.PositionT
        Final position and epoch.
    shortway : bool, optional
        Motion direction flag. If True, solves for short way direction;
        if False, solves for long way direction. Default is True.
    n_rev : int, optional
        Number of complete revolutions around the central body. Default is 0.
    forceenv : dict or ForceEnvironment, optional
        Force environment configuration containing the central body gravitational
        parameter (``'earthmu'``) and inertial reference frame (``'celestialframe'``).
        Default is ``force.deffe``.

    Returns
    -------
    pvt : posvel.PositionVelocityT
        Complete state vector (position and velocity) at times ``t1`` and ``t2``.

    Notes
    -----
    For multi-revolution transfers (``n_rev > 0``), Orekit's ``IodLambert`` solver
    currently yields only the short-period (low-energy) branch and does not return
    the corresponding long-period (high-energy) solution.

    See Also
    --------
    tellurion.ork.convert._tspvc : Convert PositionT/PositionVelocityT
    to Orekit coordinates.
    tellurion.ork.convert._pvt : Convert Orekit Orbit objects to PositionVelocityT.

    """
    # Convert inputs to Orekit TimeStampedPVCoordinates
    tspvc1 = convert._tspvc(p1)
    tspvc2 = convert._tspvc(p2)

    # Estimate orbit at t1
    solver = IodLambert(qu.sifloat(forceenv["earthmu"]))
    orbit1 = solver.estimate(
        forceenv["celestialframe"],
        shortway,
        n_rev,
        tspvc1.getPosition(),
        tspvc1.getDate(),
        tspvc2.getPosition(),
        tspvc2.getDate(),
    )

    # Propagate state to t2
    dt = tspvc2.getDate().durationFrom(tspvc1.getDate())
    orbit2 = orbit1.shiftedBy(dt)

    # Convert Orekit Orbit objects back to a single PositionVelocityT with two elements
    return convert._pvt(orbit1).concatenate(convert._pvt(orbit2))
