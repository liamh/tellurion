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
    earthloc,
    equinoctial,
    iscircels,
    isequels,
    iskepels,
    kepler,
    magdiff,
    nullisland,
    observer_location,
    pvtcart,
    spacetrack_latest,
)

np.set_printoptions(suppress=True, precision=4, linewidth=np.inf)

try:
    __version__ = _version("tellurion")
except PackageNotFoundError:
    __version__ = "0+unknown"

# Whether Orekit dependency appears importable; does not initialize JVM.
_ORK_AVAILABLE = _importlib_util.find_spec("orekit_jpype") is not None


def _ork_attr(name):
    """Fetch a lazy Orekit-backed symbol from tellurion.ork."""
    if not _ORK_AVAILABLE:
        raise RuntimeError(
            "Orekit support is not available in this environment. "
            "Install `orekit-jpype` and Orekit data to use Orekit-backed APIs."
        )
    try:
        return getattr(_import_module("tellurion.ork"), name)
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Orekit support is not available in this environment. "
            "Install `orekit-jpype` and Orekit data to use Orekit-backed APIs."
        ) from exc


def init_orekit():
    """Initialize JVM + Orekit data explicitly."""
    return _ork_attr("init_orekit")()


def orekit_available():
    """Return True if Orekit initialization succeeds, otherwise False."""
    if not _ORK_AVAILABLE:
        return False
    try:
        return _ork_attr("orekit_available")()
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

    if not _ORK_AVAILABLE:
        raise AttributeError(f"module 'tellurion' has no attribute '{name}'")

    ork = _import_module("tellurion.ork")
    try:
        value = getattr(ork, name)
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
    "dragforce",
    "earthloc",
    "elementval",
    "equinoctial",
    "from_array",
    "fromtime",
    "gravconstunits",
    "hcat",
    "hstack",
    "init_orekit",
    "iscircels",
    "isequels",
    "iskepels",
    "kepler",
    "kepleranalytic",
    "magdiff",
    "make_quantity",
    "normalizeangle",
    "orkunits",
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
    "siunits",
    "sifloat",
    "sma",
    "spacetrack_latest",
    "striptime",
    "tc",
    "time_concat",
    "timesec",
    "to_array",
    "tselements",
    "tq",
    "vstack",
]

for _opt_name in ("element_hdf5", "posvel_hdf5", "spacetrack_hdf5"):
    if hasattr(_core, _opt_name):
        __all__.append(_opt_name)
