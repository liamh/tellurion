import collections.abc

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.attitudes import LofOffset
from org.orekit.frames import LOFType
from org.orekit.orbits import CartesianOrbit
from org.orekit.propagation import SpacecraftState

from tellurion.astro import units
from tellurion.ork import convert, prop


def ntw(relsc, refsc, unitlookup=units.prefunits):
    """Find the NTW (normal, in-track, cross-track) relative
    coordinates of relsc with respect to refsc; each is a
    PositionVelocityTime or org.orekit.propagation.SpacecraftState (as
    from the `spacecraftstate` argument of `prop()` == `True`)
    """
    rf = _relframe(LOFType.NTW, relsc, refsc, unitlookup)
    rf.name = "NTW"
    return rf


def lvlh(relsc, refsc, unitlookup=units.prefunits):
    """Find the LVLH or RSW (radial, along-track, cross-track)
    relative coordinates of relsc with respect to refsc; each is a
    PositionVelocityTime or org.orekit.propagation.SpacecraftState (as
    from the `spacecraftstate` argument of `prop()` == `True`)
    """
    rf = _relframe(LOFType.LVLH, relsc, refsc, unitlookup)
    rf.name = "LVLH"
    return rf


def inv_ntw(relpvt, refsc):
    """Inverse of ntw(): convert NTW relative PVT + inertial reference
    SpacecraftState to secondary inertial SpacecraftState.
    """
    return _inv_relframe(LOFType.NTW, relpvt, refsc)


def inv_lvlh(relpvt, refsc):
    """Inverse of lvlh(): convert LVLH relative PVT + inertial reference
    SpacecraftState to secondary inertial SpacecraftState.
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


# import astropy.coordinates as coord
# ntw_position = demoa.prop.ntw_4x4hpB01_to_4x4[12]['NTW position']
# cart12 = coord.CartesianRepresentation(ntw_position)
# coord.SphericalRepresentation.from_cartesian(cart12)
# [2025-09-12 Fri 13:42] .spherical method for PVT will do this, want
# to cyclic permute to WNT?
