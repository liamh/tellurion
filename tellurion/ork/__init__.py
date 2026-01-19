import orekit_jpype as orekit
orekit.initVM()
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data()

# Import all contents from each module
from tellurion.ork.force import *
from tellurion.ork.element import *
from tellurion.ork.prop import *
from tellurion.ork.relative import *
from tellurion.ork.geog import *
from tellurion.ork.obs import *
