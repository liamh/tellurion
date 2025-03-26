import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

from . import (
    force,
    element,
    tle,
    prop,
    geog
)

from .force import *
from .element import *
from .tle import *
from .prop import *
from .geog import *

# Define the __all__ variable
__all__ = ["force", "element", "tle", "prop", "geog"]
