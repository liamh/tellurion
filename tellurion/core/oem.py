"""Optional CCSDS Orbit Ephemeris Message support."""

import importlib

import astropy.units as u
import numpy as np

from tellurion.core.frame import FrameDescriptor, FrameKind
from tellurion.core.posvel import pvtcart


_OEM_UNITS = {
    "time": u.s,
    "length": u.km,
    "speed": u.km / u.s,
    "angle": u.deg,
    "angular speed": u.rad / u.s,
    "dimensionless": u.dimensionless_unscaled,
}

_INERTIAL_AXES = {
    "earthicrf",
    "earthmj2000ec",
    "earthmj2000eq",
    "eme2000",
    "gcrf",
    "icrf",
    "j2000",
    "lunamj2000eq",
    "m50",
    "mod",
    "moonmj2000eq",
    "teme",
    "tod",
}
_LOCAL_AXES = {"lvlh", "ntw", "rsw", "rtn", "tnw", "vnc"}


def _frame_kind(axes):
    key = " ".join(axes.casefold().replace("_", " ").replace("-", " ").split())
    if key in _INERTIAL_AXES:
        return FrameKind.INERTIAL
    if key in _LOCAL_AXES:
        return FrameKind.LOCAL_ORBITAL
    if key.startswith(("itrf", "iau ", "pck")) or key == "tdr":
        return FrameKind.BODY_FIXED
    return FrameKind.CUSTOM


def _center_name(center):
    key = " ".join(
        center.casefold().replace("_", " ").replace("-", " ").split()
    )
    aliases = {
        "earth": "Earth",
        "earth barycenter": "Earth Barycenter",
        "moon": "Moon",
        "luna": "Moon",
        "sun": "Sun",
        "ssb": "Solar System Barycenter",
        "solar system barycenter": "Solar System Barycenter",
        "earth moon barycenter": "Earth-Moon Barycenter",
    }
    return aliases.get(key, center)


def oempvt(file):
    """Read a CCSDS OEM file as frame-tagged PVT state(s).

    The optional ``oem`` package is imported only when this function is
    called; install it with ``tellurion[oem]``. OEM positions and velocities
    use km and km/s. Each segment's ``CENTER_NAME`` and ``REF_FRAME`` form
    its :class:`~tellurion.FrameDescriptor`.

    Returns a single PVT when all adjacent segments share a frame, or a
    tuple of PVT objects when the file changes frame. Each contiguous run
    is kept separate so no states from different frames are concatenated.
    Unknown origins and reference frames are preserved in custom descriptors.

    Parameters
    ----------
    file : str or path-like
        Path to a KVN or XML OEM file.

    Returns
    -------
    PositionVelocityT or tuple of PositionVelocityT
        The ephemeris data, with each PVT carrying its frame descriptor.
    """
    try:
        oem = importlib.import_module("oem")
    except ModuleNotFoundError as exc:
        if exc.name != "oem":
            raise
        raise ImportError(
            "OEM reading requires the optional 'oem' package; "
            "install it with `pip install tellurion[oem]`"
        ) from exc

    message = oem.OrbitEphemerisMessage.open(file)
    segments = message.segments
    if not segments:
        raise ValueError("OEM file contains no ephemeris segments")

    result = []
    current_frame = None
    current_pvt = None
    for segment in segments:
        center = segment.metadata["CENTER_NAME"]
        axes = segment.metadata["REF_FRAME"]
        frame = FrameDescriptor(
            _center_name(center), axes, kind=_frame_kind(axes)
        )
        states = list(segment.states)
        if not states:
            raise ValueError("OEM segment contains no states")
        cartesian = np.asarray(
            [np.concatenate((state.position, state.velocity)) for state in states]
        )
        times = [state.epoch for state in states]
        pvt = pvtcart(cartesian, times, _OEM_UNITS, frame=frame)

        if current_frame == frame:
            current_pvt.concatenate(pvt)
        else:
            if current_pvt is not None:
                result.append(current_pvt)
            current_frame = frame
            current_pvt = pvt
    result.append(current_pvt)
    return result[0] if len(result) == 1 else tuple(result)
