import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

from . import (
     posvel,
     element,
     tle,
     force,
     prop,
     geog,
)

from .posvel import *
from .element import *
from .tle import *
from .force import *
from .prop import *
from .geog import *

# Define the __all__ variable
__all__ = ["init", "posvel", "element", "tle", "force", "prop", "geog"]
