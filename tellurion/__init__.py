"""Top-level Tellurion public API."""

from importlib import import_module as _import_module, util as _importlib_util
from importlib.metadata import PackageNotFoundError, version as _version

import numpy as np

import tellurion.core as _core
from tellurion.astro import (
    abstime,
    change_units,
    from_array,
    fromtime,
    gravconstunits,
    hcat,
    hstack,
    make_quantity,
    normalizeangle,
    orkunits,
    prefunits,
    prefnumabstime,
    quantity_to_array,
    quantity_to_dict,
    rev,
    revolution,
    siunits,
    sifloat,
    striptime,
    tc,
    time_concat,
    timesec,
    to_array,
    tq,
    vstack,
)
from tellurion.core import (
    PVT,
    EarthObservationT,
    ElementSetT,
    MeanElementSetT,
    PositionBase,
    PositionT,
    PositionVelocityT,
    azelrange,
    circular,
    circeltma_names,
    circeltta_names,
    circtimeelements,
    earthloc,
    elev,
    equeltma_names,
    equeltta_names,
    equinoctial,
    equtimeelements,
    gnuserid,
    iscircels,
    isequels,
    iskepels,
    kepeltma_names,
    kepeltta_names,
    kepler,
    location,
    magdiff,
    nullisland,
    obsdict1,
    observer_location,
    pvtcart,
    satdata,
    sfdict,
    spacetrack_latest,
    statefnval,
    stscdata,
    timeelements,
    userid,
)

np.set_printoptions(suppress=True, precision=4, linewidth=np.inf)

try:
    __version__ = _version("tellurion")
except PackageNotFoundError:
    __version__ = "0+unknown"

def _orekit_bridge_available():
    """Return True when orekit-jpype bridge modules are importable."""
    if _importlib_util.find_spec("orekit_jpype") is None:
        return False
    try:
        _import_module("orekit_jpype")
        _import_module("orekit_jpype.pyhelpers")
    except Exception:
        return False
    return True


# Whether Orekit bridge modules appear importable; does not initialize JVM.
_ORK_AVAILABLE = _orekit_bridge_available()


def _is_missing_orekit_bridge_error(exc):
    current = exc
    seen = set()
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        if isinstance(current, ModuleNotFoundError):
            if (current.name or "").startswith("orekit_jpype"):
                return True
        elif isinstance(current, ImportError):
            if (getattr(current, "name", "") or "").startswith("orekit_jpype"):
                return True
        current = current.__cause__ or current.__context__
    return False


def _is_orekit_init_error(exc):
    text = str(exc)
    return (
        "Failed to initialize Orekit JVM/data." in text
        or "Orekit initialization previously failed." in text
    )


def _raise_orekit_runtime_error(exc):
    raise RuntimeError(
        "Orekit support is not available or failed to initialize. "
        "Install/configure Java + orekit-jpype, then retry."
    ) from exc


def _ork_attr(name):
    """Fetch a lazy Orekit-backed symbol from tellurion.ork."""
    if not _ORK_AVAILABLE:
        _raise_orekit_runtime_error(RuntimeError("orekit_jpype not available"))
    try:
        value = getattr(_import_module("tellurion.ork"), name)
    except RuntimeError as exc:
        if _is_missing_orekit_bridge_error(exc) or _is_orekit_init_error(exc):
            _raise_orekit_runtime_error(exc)
        raise
    except ModuleNotFoundError as exc:
        if _is_missing_orekit_bridge_error(exc):
            _raise_orekit_runtime_error(exc)
        raise
    if callable(value):
        def _wrapped(*args, **kwargs):
            try:
                return value(*args, **kwargs)
            except RuntimeError as exc:
                if _is_missing_orekit_bridge_error(exc) or _is_orekit_init_error(exc):
                    _raise_orekit_runtime_error(exc)
                raise
        return _wrapped
    return value


def init_orekit():
    """Initialize JVM + Orekit data explicitly."""
    return _ork_attr("init_orekit")()


def orekit_available():
    """Return True if Orekit initialization succeeds, otherwise False."""
    if not _ORK_AVAILABLE:
        return False
    try:
        ork = _import_module("tellurion.ork")
        return bool(getattr(ork, "orekit_available")())
    except RuntimeError:
        return False


def prepare(*args, **kwargs):
    """Prepare an Orekit-backed propagation."""
    return _ork_attr("prepare")(*args, **kwargs)


def propagate(*args, **kwargs):
    """Execute a prepared Orekit-backed propagation."""
    return _ork_attr("propagate")(*args, **kwargs)


def setgravity(*args, **kwargs):
    """Create an Orekit force environment."""
    return _ork_attr("setgravity")(*args, **kwargs)


def dragforce(*args, **kwargs):
    """Add atmospheric drag to an Orekit force environment."""
    return _ork_attr("dragforce")(*args, **kwargs)


def siderealtime(*args, **kwargs):
    """Compute sidereal time using Orekit reference frames."""
    return _ork_attr("siderealtime")(*args, **kwargs)


def elementval(*args, **kwargs):
    """Compute orbital element values from a state representation."""
    return _ork_attr("elementval")(*args, **kwargs)


def pvt(*args, **kwargs):
    """Convert an Orekit/object representation to PositionVelocityT."""
    return _ork_attr("pvt")(*args, **kwargs)


def allplane(*args, **kwargs):
    """Convert between equivalent orbital element input pairs."""
    return _ork_attr("allplane")(*args, **kwargs)


def sma(*args, **kwargs):
    """Compute semimajor axis from related orbit parameters."""
    return _ork_attr("sma")(*args, **kwargs)


def kepleranalytic(*args, **kwargs):
    """Create a Kepler-only Orekit force environment."""
    return _ork_attr("kepleranalytic")(*args, **kwargs)


def tselements(*args, **kwargs):
    """Build a TimeSeries of orbital elements from one or more states."""
    return _ork_attr("tselements")(*args, **kwargs)


def __getattr__(name):
    """Fallback lazy Orekit symbol exposure for advanced/internal use."""
    if name in ("element_hdf5", "posvel_hdf5", "spacetrack_hdf5"):
        if hasattr(_core, name):
            value = getattr(_core, name)
            globals()[name] = value
            return value
        raise AttributeError(f"module 'tellurion' has no attribute '{name}'")

    try:
        value = _ork_attr(name)
    except AttributeError as exc:
        raise AttributeError(f"module 'tellurion' has no attribute '{name}'") from exc
    globals()[name] = value
    return value


__all__ = [
    "__version__",
    "ElementSetT",
    "EarthObservationT",
    "MeanElementSetT",
    "PVT",
    "PositionBase",
    "PositionT",
    "PositionVelocityT",
    "abstime",
    "allplane",
    "azelrange",
    "change_units",
    "circular",
    "circeltma_names",
    "circeltta_names",
    "circtimeelements",
    "dragforce",
    "earthloc",
    "elev",
    "elementval",
    "equeltma_names",
    "equeltta_names",
    "equinoctial",
    "equtimeelements",
    "from_array",
    "fromtime",
    "gnuserid",
    "gravconstunits",
    "hcat",
    "hstack",
    "init_orekit",
    "iscircels",
    "isequels",
    "iskepels",
    "kepeltma_names",
    "kepeltta_names",
    "kepler",
    "kepleranalytic",
    "location",
    "magdiff",
    "make_quantity",
    "normalizeangle",
    "orkunits",
    "obsdict1",
    "prefunits",
    "nullisland",
    "observer_location",
    "orekit_available",
    "prefnumabstime",
    "prepare",
    "propagate",
    "pvt",
    "pvtcart",
    "quantity_to_array",
    "quantity_to_dict",
    "rev",
    "revolution",
    "setgravity",
    "siderealtime",
    "satdata",
    "sfdict",
    "siunits",
    "sifloat",
    "sma",
    "spacetrack_latest",
    "statefnval",
    "stscdata",
    "striptime",
    "tc",
    "timeelements",
    "time_concat",
    "timesec",
    "to_array",
    "tselements",
    "tq",
    "userid",
    "vstack",
]

for _opt_name in ("element_hdf5", "posvel_hdf5", "spacetrack_hdf5"):
    if hasattr(_core, _opt_name):
        __all__.append(_opt_name)
