import collections.abc

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.attitudes import LofOffset
from org.orekit.frames import LOFType
from org.orekit.orbits import CartesianOrbit
from org.orekit.propagation import SpacecraftState

from tellurion.astro import units
from tellurion.core.frame import FrameDescriptor, FrameKind
from tellurion.ork import convert, prop


def ntw(relsc, refsc, unitlookup=units.prefunits):
    """Relative state of a secondary spacecraft in the NTW frame.

    The NTW local orbital frame is centered on the reference spacecraft
    and has axes:

    * **N** -- in the orbital plane, normal to the velocity vector
    * **T** -- tangential, along the velocity vector (in-track)
    * **W** -- along the orbital angular momentum (cross-track)

    Parameters
    ----------
    relsc : `~tellurion.PositionVelocityTime` or `org.orekit.propagation.SpacecraftState`
        The secondary ("deputy"/"chaser") state(s) whose relative
        coordinates are wanted.  May be a time series.
    refsc : `~tellurion.PositionVelocityTime` or `org.orekit.propagation.SpacecraftState`
        The reference ("chief"/"target") state(s) defining the frame
        origin and orientation.  Must be inertial and have the same
        epochs as ``relsc`` (the same length if a time series).
    unitlookup : callable or mapping, optional
        Selects the units of the output quantities; default
        `tellurion.astro.units.prefunits`.  Pass `tellurion.siunits`
        for SI units.

    Returns
    -------
    `~tellurion.PositionVelocityTime`
        Position and velocity of ``relsc`` relative to ``refsc``,
        expressed in the NTW axes, with ``name == "NTW"``.  A time
        series if the inputs are.

    See Also
    --------
    lvlh : Relative state in the radial / along-track / cross-track frame.
    inv_ntw : Inverse transformation.

    Notes
    -----
    Velocity is the velocity of ``relsc`` relative to ``refsc`` as seen
    in the rotating frame, i.e. it includes the frame rotation
    (transport) term.

    Examples
    --------
    >>> rel = tell.ntw(pvt_b, pvt_a)           # doctest: +SKIP
    >>> rel = tell.ntw(pvt_b, pvt_a, tell.siunits)  # doctest: +SKIP
    """
    rf = _relframe(LOFType.NTW, relsc, refsc, unitlookup)
    rf.name = "NTW"
    rf.frame = _relative_frame_descriptor(refsc, "NTW")
    return rf


def lvlh(relsc, refsc, unitlookup=units.prefunits):
    """Relative state of a secondary spacecraft in the LVLH (RSW) frame.

    The local-vertical/local-horizontal frame is centered on the
    reference spacecraft and has axes:

    * **R** -- radial, from the central body toward the reference
    * **S** -- along-track, completing the right-handed set
    * **W** -- cross-track, along the orbital angular momentum

    Parameters, Returns
    -------------------
    As for `ntw`; the returned object has ``name == "LVLH"``.

    Parameters
    ----------
    relsc, refsc : `~tellurion.PositionVelocityTime` or `org.orekit.propagation.SpacecraftState`
        Secondary and reference inertial states.
    unitlookup : callable or mapping, optional
        Output unit selection; default `tellurion.astro.units.prefunits`.

    Returns
    -------
    `~tellurion.PositionVelocityTime`
        Relative position and velocity in the RSW axes.

    See Also
    --------
    ntw, inv_lvlh

    Notes
    -----
    This is the convention of Curtis, *Orbital Mechanics for Engineering
    Students*, 4th ed., Example 7.1 (see ``demos/propdemo.py``).
    """
    rf = _relframe(LOFType.LVLH, relsc, refsc, unitlookup)
    rf.name = "LVLH"
    rf.frame = _relative_frame_descriptor(refsc, "LVLH")
    return rf



def inv_ntw(relpvt, refsc):
    """Convert NTW relative states back to inertial states.

    Inverse of `ntw`.

    Parameters
    ----------
    relpvt : `~tellurion.PositionVelocityTime`
        Relative state(s) in NTW axes, as returned by `ntw`.
    refsc : `~tellurion.PositionVelocityTime` or `org.orekit.propagation.SpacecraftState`
        The inertial reference state(s) used to define the frame, at the
        same epochs as ``relpvt``.

    Returns
    -------
    `~tellurion.PositionVelocityTime`
        Inertial state(s) of the secondary spacecraft.

    See Also
    --------
    ntw
    """
    return _inv_relframe(LOFType.NTW, relpvt, refsc)


def inv_lvlh(relpvt, refsc):
    """Convert LVLH (RSW) relative states back to inertial states.

    Inverse of `lvlh`.

    Parameters
    ----------
    relpvt : `~tellurion.PositionVelocityTime`
        Relative state(s) in LVLH axes, as returned by `lvlh`.
    refsc : `~tellurion.PositionVelocityTime` or `org.orekit.propagation.SpacecraftState`
        The inertial reference state(s), at the same epochs as ``relpvt``.

    Returns
    -------
    `~tellurion.PositionVelocityTime`
        Inertial state(s) of the secondary spacecraft.

    See Also
    --------
    lvlh
    """
    return _inv_relframe(LOFType.LVLH, relpvt, refsc)


def _relframe(frame, relsc, refsc, unitlookup=units.prefunits):
    ref = prop._to_spacecraft_state(refsc)
    rel = prop._to_spacecraft_state(relsc)
    if isinstance(ref, collections.abc.Iterable):
        rf = [
            _relframe(frame, rel, ref, unitlookup) for (rel, ref) in zip(rel, ref)
        ]
        return rf[0].merge(rf[1:])
    # See https://forum.orekit.org/t/1478/2
    relframe = LofOffset(ref.getFrame(), frame)
    converted = SpacecraftState(
        ref.getOrbit(),
        relframe.getAttitude(ref.getOrbit(), ref.getDate(), ref.getFrame()),
    )
    pvrel = converted.toTransform().transformPVCoordinates(
        rel.getPVCoordinates(ref.getFrame())
    )
    return convert._pvt(pvrel)


def _inv_relframe(frame, relsc, refsc):
    ref = prop._to_spacecraft_state(refsc)
    if isinstance(ref, collections.abc.Iterable):
        inertial = [_inv_relframe(frame, rel, ref) for rel, ref in zip(relsc, ref)]
        return inertial[0].merge(inertial[1:])

    relframe = LofOffset(ref.getFrame(), frame)
    converted = SpacecraftState(
        ref.getOrbit(),
        relframe.getAttitude(ref.getOrbit(), ref.getDate(), ref.getFrame()),
    )

    # Use tellurion conversion helper to build TimeStampedPVCoordinates
    rel_tspvc = convert._tspvc(relsc.pvt())

    # _relframe uses inertial->LOF; inverse is LOF->inertial
    sec_pv = converted.toTransform().getInverse().transformPVCoordinates(rel_tspvc)

    sec_orbit = CartesianOrbit(
        sec_pv,
        ref.getFrame(),
        ref.getDate(),
        ref.getOrbit().getMu(),
    )
    return convert._pvt(SpacecraftState(sec_orbit))


def _relative_frame_descriptor(refsc, axes):
    reference = getattr(refsc, "name", None)
    return FrameDescriptor(
        origin=reference or "reference spacecraft",
        axes=axes,
        kind=FrameKind.LOCAL_ORBITAL,
        reference=reference,
    )


# import astropy.coordinates as coord
# ntw_position = demoa.prop.ntw_4x4hpB01_to_4x4[12]['NTW position']
# cart12 = coord.CartesianRepresentation(ntw_position)
# coord.SphericalRepresentation.from_cartesian(cart12)
# [2025-09-12 Fri 13:42] .spherical method for PVT will do this, want
# to cyclic permute to WNT?
