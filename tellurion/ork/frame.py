"""Orekit frame lookup and state transformations."""

from tellurion.ork import ensure_orekit_initialized

ensure_orekit_initialized()

from org.orekit.frames import FramesFactory

from tellurion.core.frame import FrameDescriptor, FrameKind
from tellurion.core.posvel import PositionVelocityT

from . import convert


_EARTH = {"earth"}
_SSB = {"solar system barycenter", "solar-system barycenter", "ssb"}


def _key(value):
    return " ".join(value.casefold().replace("_", " ").split())


def resolve_orekit_frame(descriptor):
    """Resolve a supported inertial frame descriptor to an Orekit ``Frame``.

    Currently supported descriptors are Earth-centered EME2000 and
    solar-system-barycentric ICRF. Other origins, axes, and frame kinds are
    rejected rather than guessed.
    """
    if not isinstance(descriptor, FrameDescriptor):
        raise TypeError("descriptor must be a FrameDescriptor")
    if descriptor.kind is not FrameKind.INERTIAL:
        raise ValueError(
            f"Orekit resolution is not supported for {descriptor.kind.value!r} frames"
        )

    origin = _key(descriptor.origin)
    axes = _key(descriptor.axes)
    if origin in _EARTH and axes in {"eme2000", "j2000"}:
        return FramesFactory.getEME2000()
    if origin in _SSB and axes == "icrf":
        return FramesFactory.getICRF()
    raise ValueError(
        "Unsupported inertial frame descriptor "
        f"(origin={descriptor.origin!r}, axes={descriptor.axes!r})"
    )


def transform_pvt(pvt, target_frame):
    """Transform a PVT state or series to a supported inertial frame.

    Orekit's date-dependent frame transform is applied to the full position
    and velocity state, including origin translation and transport velocity.
    The input PVT must have a frame descriptor and an epoch.
    """
    if not isinstance(pvt, PositionVelocityT):
        raise TypeError("pvt must be a PositionVelocityT")
    if pvt.frame is None:
        raise ValueError("pvt must have a frame descriptor")
    if pvt.time is None:
        raise ValueError("pvt must have an epoch")
    if not isinstance(target_frame, FrameDescriptor):
        raise TypeError("target_frame must be a FrameDescriptor")

    source = resolve_orekit_frame(pvt.frame)
    target = resolve_orekit_frame(target_frame)

    def transform_one(state):
        date = convert._abstime_to_okad(state.time)
        transform = source.getTransformTo(target, date)
        transformed = transform.transformPVCoordinates(convert._tspvc(state))
        result = convert._pvt(transformed)
        result.frame = target_frame
        result.aux = state.aux.copy()
        if hasattr(state, "name"):
            result.name = state.name
        return result

    if pvt.isscalar:
        return transform_one(pvt)
    if len(pvt) == 0:
        raise ValueError("Cannot transform an empty PVT series")

    result = transform_one(pvt[0])
    for index in range(1, len(pvt)):
        result.concatenate(transform_one(pvt[index]))
    if hasattr(pvt, "name"):
        result.name = pvt.name
    return result
