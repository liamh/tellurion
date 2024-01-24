# Set up for orbit computations from Orekit
# Load with: from pyork import *
import sys

# Easy way to see what is defined for an object
from inspect import getmembers

# Orekit setup
from math import radians, degrees
import orekit
vm = orekit.initVM()
print ('Python version:',sys.version)
print ('Java version:',vm.java_version)
print ('Orekit version:', orekit.VERSION)
from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
setup_orekit_curdir()

# Orbital elements and PVT
from org.orekit.utils import PVCoordinates, TimeStampedPVCoordinates
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.orbits import CartesianOrbit, KeplerianOrbit, PositionAngleType, OrbitType
from org.orekit.frames import FramesFactory
from org.orekit.utils import Constants

# Propagation and ephemeris
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState, EphemerisGenerator
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.utils import IERSConventions
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from orekit import JArray_double


utc = TimeScalesFactory.getUTC()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earthrad = Constants.IERS2010_EARTH_EQUATORIAL_RADIUS
earthangspd = Constants.IERS2010_EARTH_ANGULAR_VELOCITY
earthJ2 = -Constants.IERS2010_EARTH_C20
earthmu = Constants.IERS2010_EARTH_MU


#def nowutc():
#    now = datetime.datetime.now(datetime.UTC)
