import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

from . import (
    force,
    element,
    prop,
    relative,
    geog,
    obs
)

from tellurion.ork.force import *
from tellurion.ork.element import *
from tellurion.ork.prop import *
from tellurion.ork.relative import *
from tellurion.ork.geog import *
from tellurion.ork.obs import *

# Define the __all__ variable
__all__ = ["force", "element", "prop", "relative", "geog", "obs"]
