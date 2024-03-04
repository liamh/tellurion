from orkinit import *

from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from orekit import JArray_double
from org.orekit.forces.drag import AbstractDragForceModel, DragForce
from org.orekit.models.earth.atmosphere import Atmosphere, HarrisPriester, DTM2000, NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.forces.drag import IsotropicDrag

from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory
OneAxisEllipsoid.__repr__ = \
    lambda self: f"<Near-spherical body equatorial radius {self.equatorialRadius}m, polar radius difference {-self.equatorialRadius*self.flattening}m >"

# Default celestial environment constants
celestdflt = {
        'earthframe': FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        'celestialframe': FramesFactory.getEME2000(),
        'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY,
        'earthflat': Constants.IERS2010_EARTH_FLATTENING,
        'sun': CelestialBodyFactory.getSun(),
}

def setgravity(degree, order, mass = 100.0):
    eg = GravityFieldFactory.getNormalizedProvider(degree, order)
    force = {'gravity': eg,
             'gravity-degree-order': [eg.maxDegree, eg.maxOrder],
             'earthrad': eg.ae,
             'earthmu': eg.mu,
             # This isn't necessarily what it's using, how to get actual coefficient used?
             'zonals': [-Constants.IERS2010_EARTH_C20],
             'mass': mass} # Needed by several forces but not gravity
    # Define spherical altitude for convience, not specifically force related, but uses the definitions
    global sphalt
    sphalt = OneAxisEllipsoid(force['earthrad'], 0.0, celestdflt['earthframe'])
    global forcedflt
    forcedflt = celestdflt | force

# Set the default gravity 0x0
setgravity(0, 0)

swdata = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

# Optionally add a drag force
# Default spacecraft parameters gives B = C_D A/m = 0.01 m^2/kg
def dragforce(atmdensname = None, dragcoef = 1.0, dragarea = 1.0):
    earth = OneAxisEllipsoid(forcedflt['earthrad'], forcedflt['earthflat'],  forcedflt['earthframe'])
    match atmdensname:
        case 'hp': # Harris-Priester atmospheric density model
            atmdens = HarrisPriester(forcedflt['sun'], earth)
        case 'dtm': # DTM2000 atmospheric density model
            atmdens = DTM2000(swdata, forcedflt['sun'], earth)
        case 'msis':
            atmdens = NRLMSISE00(swdata, forcedflt['sun'], earth)
        case None:
            return None
    force = {'mass': mass,  # The models need a spacecraft mass, unit kg.
             'dragarea': dragarea, # Cross-sectional area perpendicular to atmosphere direction, m^2
             'dragcoef': dragcoef, # Coefficient of drag
             'B': dragarea*dragcoef/mass,
             'atmdens': atmdens,
             'dragforce': DragForce(atmdens, IsotropicDrag(dragarea, dragcoef))}
    forcedflt.update(force)

# def atmdens(location, time, model = 'hp'):
#     return(envct[model].getDensity(time, location, envct['celestialframe']))
