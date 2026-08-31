"""Public exports for tellurion.core."""

from importlib import import_module

from tellurion.core.element import (
    ElementSetT,
    circeltma_names,
    circeltta_names,
    circtimeelements,
    circular,
    equeltma_names,
    equeltta_names,
    equinoctial,
    equtimeelements,
    iscircels,
    isequels,
    iskepels,
    kepeltma_names,
    kepeltta_names,
    kepler,
    sfdict,
    statefnval,
    timeelements,
)
from tellurion.core.geog import earthloc, nullisland, observer_location
from tellurion.core.geonames import elev, gnuserid, location, userid
from tellurion.core.obs import EarthObservationT, azelrange, obsdict1
from tellurion.core.posvel import (
    PVT,
    PositionBase,
    PositionT,
    PositionVelocityT,
    pvtcart,
)
from tellurion.core.spacetrack import (
    MeanElementSetT,
    satdata,
    spacetrack_latest,
    stscdata,
)
from tellurion.core.util import magdiff

__all__ = [
    "ElementSetT",
    "EarthObservationT",
    "MeanElementSetT",
    "PVT",
    "PositionBase",
    "PositionT",
    "PositionVelocityT",
    "azelrange",
    "circular",
    "circeltma_names",
    "circeltta_names",
    "circtimeelements",
    "earthloc",
    "elev",
    "equeltma_names",
    "equeltta_names",
    "equinoctial",
    "equtimeelements",
    "gnuserid",
    "iscircels",
    "isequels",
    "iskepels",
    "kepeltma_names",
    "kepeltta_names",
    "kepler",
    "location",
    "magdiff",
    "nullisland",
    "obsdict1",
    "observer_location",
    "pvtcart",
    "satdata",
    "sfdict",
    "spacetrack_latest",
    "statefnval",
    "stscdata",
    "timeelements",
    "userid",
]

# Optional HDF5 serialization modules are exposed when installed.
try:
    import astropy_hdf5io  # noqa: F401

    element_hdf5 = import_module("tellurion.core.element_hdf5")
    posvel_hdf5 = import_module("tellurion.core.posvel_hdf5")
    spacetrack_hdf5 = import_module("tellurion.core.spacetrack_hdf5")

    __all__.extend(["element_hdf5", "posvel_hdf5", "spacetrack_hdf5"])
except ImportError:
    pass
