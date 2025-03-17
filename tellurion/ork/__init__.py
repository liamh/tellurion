import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

from . import (
    force,
    convert,
    element,
    tle,
    prop,
    geog
)

from .force import *
from .convert import *
from .element import *
from .tle import *
from .prop import *
from .geog import *

# Define the __all__ variable
__all__ = ["force", "convert", "element", "tle", "prop", "geog"]
