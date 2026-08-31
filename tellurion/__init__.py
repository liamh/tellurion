from importlib import import_module, util as importlib_util
from importlib.metadata import PackageNotFoundError, version

import numpy as np

np.set_printoptions(suppress=True, precision=4, linewidth=np.inf)

# Import everything from astro (AstroPy extensions) and core
from tellurion.astro import *
from tellurion.core import *

try:
    __version__ = version("tellurion")
except PackageNotFoundError:
    __version__ = "0+unknown"

# Whether Orekit dependency appears importable; does not initialize JVM.
_ORK_AVAILABLE = importlib_util.find_spec("orekit_jpype") is not None


def __getattr__(name):
    """Lazily expose Orekit symbols from tellurion.ork on first use."""
    ork = import_module("tellurion.ork")
    try:
        value = getattr(ork, name)
    except AttributeError as exc:
        raise AttributeError(f"module 'tellurion' has no attribute '{name}'") from exc
    globals()[name] = value
    return value

# Optionally expose availability status
__all__ = ["_ORK_AVAILABLE", "__version__"]
