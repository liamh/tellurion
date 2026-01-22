import numpy as np
np.set_printoptions(suppress=True, precision=4, linewidth=np.inf)

# Import everything from astro (AstroPy extensions) and core
from tellurion.astro import *
from tellurion.core import *

# Try to import ork modules if Java/jpype is available
try:
    from tellurion.ork import *
    _ORK_AVAILABLE = True
except (ImportError, OSError) as e:
    # OSError catches jpype-specific errors when JVM fails to start
    _ORK_AVAILABLE = False
    import warnings
    warnings.warn(
        f"tellurion.ork modules are not available: {e}.  "
        "Only core functionality will be loaded.",
        ImportWarning
    )

# Optionally expose availability status
__all__ = ['_ORK_AVAILABLE']
