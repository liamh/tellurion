# ruff: noqa: F403, E402

from importlib import import_module as _import_module, util as _importlib_util
from importlib.metadata import PackageNotFoundError, version as _version

import numpy as np

np.set_printoptions(suppress=True, precision=4, linewidth=np.inf)

# Import everything from astro (AstroPy extensions) and core
from tellurion.astro import *
from tellurion.core import *

try:
    __version__ = _version("tellurion")
except PackageNotFoundError:
    __version__ = "0+unknown"

# Whether Orekit dependency appears importable; does not initialize JVM.
_ORK_AVAILABLE = _importlib_util.find_spec("orekit_jpype") is not None


def __getattr__(name):
    """Lazily expose Orekit symbols from tellurion.ork on first use."""
    ork = _import_module("tellurion.ork")
    try:
        value = getattr(ork, name)
    except AttributeError as exc:
        raise AttributeError(f"module 'tellurion' has no attribute '{name}'") from exc
    globals()[name] = value
    return value
