# Import all contents from each module
from tellurion.core.util import *
from tellurion.core.posvel import *
from tellurion.core.element import *
from tellurion.core.geog import *
from tellurion.core.obs import *
from tellurion.core.geonames import *
from tellurion.core.spacetrack import *

# Optional HDF5 serialization — only available if tellurion[hdf5] is installed
try:
    import astropy_hdf5io  # noqa: F401 - registers AstroPy serializers as side effect
    from tellurion.core.posvel_hdf5 import *
    from tellurion.core.element_hdf5 import *
    from tellurion.core.spacetrack_hdf5 import *
except ImportError:
    pass