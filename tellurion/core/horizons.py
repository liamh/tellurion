"""Optional JPL Horizons state-vector support."""

import importlib
import re

import astropy.units as u
import numpy as np
from astropy.time import Time

from tellurion.astro.time import abstime
from tellurion.core.frame import FrameDescriptor, FrameKind
from tellurion.core.posvel import pvtcart


_CENTER_NAMES = {
    "0": "Solar System Barycenter",
    "3": "Earth-Moon Barycenter",
    "10": "Sun",
    "199": "Mercury",
    "299": "Venus",
    "399": "Earth",
    "301": "Moon",
    "499": "Mars",
    "599": "Jupiter",
    "699": "Saturn",
    "799": "Uranus",
    "899": "Neptune",
    "999": "Pluto",
}


def _center_origin(center):
    match = re.fullmatch(r"\s*@?(\d+)\s*", str(center))
    if match:
        return _CENTER_NAMES.get(match.group(1), str(center).strip())
    return str(center).strip()


def _iso8601(time):
    """ISO 8601 string, in TDB as Horizons expects for vectors."""
    tdb = abstime(time).tdb
    tdb.precision = 3
    return tdb.isot


def _frame(center, refplane):
    origin = _center_origin(center)
    if refplane == "earth":
        axes = "ICRF/Earth-equator"
        kind = FrameKind.INERTIAL
    elif refplane == "ecliptic":
        axes = "ICRF/Ecliptic"
        kind = FrameKind.INERTIAL
    else:
        axes = "ICRF/Body-equator-of-date"
        kind = FrameKind.CUSTOM
    return FrameDescriptor(origin, axes, kind=kind)


def horizons(
    target_id="-1024",
    center="@399",
    start_time="2026-04-02T02:00:00",
    stop_time="2026-04-05",
    step="10m",
    refplane="earth",
):
    """Fetch geometric state vectors from JPL Horizons as a PVT series.

    Requires the optional ``astroquery`` package, installed with
    ``pip install "tellurion[horizons]"``. The dependency is imported only
    when this function is called.

    Parameters
    ----------
    target_id : str, optional
        Horizons target identifier. Defaults to ``"-1024"`` (Artemis II Orion).
    center : str, optional
        Horizons vector center/location, e.g. ``"@399"`` for Earth or
        ``"@301"`` for the Moon.
    start_time, stop_time : optional
        Start and stop epochs: anything `tellurion.abstime` accepts (ISO 8601
        string, `~astropy.time.Time`, `datetime`, ...). They are converted to
        ISO 8601 strings in TDB, the time scale Horizons uses for vectors.
    step : str, optional
        Horizons sampling step, such as ``"10m"`` or ``"1h"``.
    refplane : {"earth", "ecliptic", "body"}, optional
        Output reference plane. Earth and ecliptic planes are described as
        inertial ICRF-oriented axes; the body-equator-of-date option is
        tagged custom because its axes vary with date. Defaults to "earth".

    Returns
    -------
    PositionVelocityT
        State-vector time series in km and km/s, with frame metadata.
        Horizons epochs are interpreted as TDB.
    """
    try:
        horizons_module = importlib.import_module("astroquery.jplhorizons")
    except ModuleNotFoundError as exc:
        if exc.name != "astroquery":
            raise
        raise ImportError(
            "JPL Horizons queries require the optional 'astroquery' package; "
            "install it with `pip install tellurion[horizons]`"
        ) from exc

    refplane = str(refplane).casefold()
    if refplane not in {"earth", "ecliptic", "body"}:
        raise ValueError("refplane must be 'earth', 'ecliptic', or 'body'")

    obj = horizons_module.Horizons(
        id=str(target_id),
        location=str(center),
        epochs={
            "start": _iso8601(start_time),
            "stop": _iso8601(stop_time),
            "step": step,
        },
    )
    vectors = obj.vectors(refplane=refplane, aberrations="geometric")
    if len(vectors) == 0:
        raise ValueError("JPL Horizons returned no state vectors")

    columns = [
        u.Quantity(vectors[key]).to_value(unit)
        for key, unit in (
            ("x", u.km),
            ("y", u.km),
            ("z", u.km),
            ("vx", u.km / u.s),
            ("vy", u.km / u.s),
            ("vz", u.km / u.s),
        )
    ]
    cartesian = np.column_stack(columns)
    times = Time(
        u.Quantity(vectors["datetime_jd"]).to_value(u.day),
        format="jd",
        scale="tdb",
    )
    times.format = "isot"
    result = pvtcart(
        cartesian,
        times,
        {
            "time": u.s,
            "length": u.km,
            "speed": u.km / u.s,
            "angle": u.deg,
            "angular speed": u.rad / u.s,
            "dimensionless": u.dimensionless_unscaled,
        },
        frame=_frame(center, refplane),
    )
    if "targetname" in vectors.colnames:
        result.name = str(vectors["targetname"][0]).split("(")[0].strip()
    return result
