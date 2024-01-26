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
from org.orekit.propagation import Propagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory
from org.orekit.utils import IERSConventions
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from orekit import JArray_double
from org.orekit.forces.drag import AbstractDragForceModel, DragForce
from org.orekit.models.earth.atmosphere import Atmosphere, HarrisPriester, DTM2000, NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.forces.drag import IsotropicDrag

# Orekit constants
okct = {'utc': TimeScalesFactory.getUTC(),
        'earthframe': FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        'celestialframe': FramesFactory.getEME2000(),
        'earthrad': Constants.IERS2010_EARTH_EQUATORIAL_RADIUS,
        'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY,
        'earthJ2': -Constants.IERS2010_EARTH_C20,
        'earthmu': Constants.IERS2010_EARTH_MU,
        'earthflat': Constants.IERS2010_EARTH_FLATTENING,
        'meananom': PositionAngleType.MEAN,
        'trueanom': PositionAngleType.TRUE,
        'cartesian': OrbitType.CARTESIAN,
        'sun': CelestialBodyFactory.getSun(),
        'gravity10x10': GravityFieldFactory.getNormalizedProvider(10, 10), # 10x10
        'gravity0x0': GravityFieldFactory.getNormalizedProvider(0, 0), # 0x0
        'swdata': CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt"),
        'stopalt': 125.0e3  # Altitude at which propagation should stop
}
okct['gravity']=okct['gravity0x0']
okct['earth'] = OneAxisEllipsoid(okct['earthrad'], okct['earthflat'],  okct['earthframe'])
okct['sphearth'] = OneAxisEllipsoid(okct['earthrad'], 0.0,  okct['earthframe'])
# Atmospheric density models
okct['hp'] = HarrisPriester(okct['sun'], okct['earth']) # Harris-Priester atmospheric density model
okct['dtm'] = DTM2000(okct['swdata'], okct['sun'], okct['earth']) # DTM2000 atmospheric density model
okct['msis'] = NRLMSISE00(okct['swdata'], okct['sun'], okct['earth'])

# Create an example spacecraft with B = C_D A/m = 0.01 m^2/kg
scB010 = {'mass': 100.0,  # The models need a spacecraft mass, unit kg.
        'dragarea': 1.0, # Cross-sectional area perpendicular to atmosphere direction, m^2
        'dragcoef': 1.0 # Coefficient of drag
        }
scB010['drag'] = IsotropicDrag(scB010['dragarea'], scB010['dragcoef'])
scB010['atmdens'] = okct['hp']
scB010['dragforce'] = DragForce(scB010['atmdens'], scB010['drag']);

def atmdens(location, time, model = 'hp'):
    return(okct[model].getDensity(time, location, okct['celestialframe']))

import datetime

def datm(year, month, day, hour=12, minute=0, second=0.0, microsecond=0):
    return(AbsoluteDate(year, month, day, hour, minute, float(second), okct['utc']))

def nowutc():
    now = datetime.datetime.now(datetime.UTC)
    return(datm(now.year, now.month, now.day, now.hour, now.minute, now.second+1.0e-6*now.microsecond))

hour = 3600.0
day = 24.0*hour
