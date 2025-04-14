import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

from . import (
    force,
    element,
    prop,
    relative,
    geog
)

from .force import *
from .element import *
from .prop import *
from .relative import *
from .geog import *

# Define the __all__ variable
__all__ = ["force", "element", "prop", "relative", "geog"]
