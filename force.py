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
OneAxisEllipsoid.__repr__ = lambda self: f"<Spherical body of radius {self.equatorialRadius}m>"

# Default celestial environment constants
celestdflt = {
        'earthframe': FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        'celestialframe': FramesFactory.getEME2000(),
        'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY,
        'earthflat': Constants.IERS2010_EARTH_FLATTENING,
        'sun': CelestialBodyFactory.getSun(),
}

# # Atmospheric density models
# envct['hp'] = HarrisPriester(envct['sun'], envct['earth']) # Harris-Priester atmospheric density model
# envct['dtm'] = DTM2000(envct['swdata'], envct['sun'], envct['earth']) # DTM2000 atmospheric density model
# envct['msis'] = NRLMSISE00(envct['swdata'], envct['sun'], envct['earth'])

def setgravity(degree, order):
    eg = GravityFieldFactory.getNormalizedProvider(degree, order)
    force = {'gravity': eg,
             'gravity-degree-order': [eg.maxDegree, eg.maxOrder],
             'earthrad': eg.ae,
             'earthmu': eg.mu,
             # This isn't necessarily what it's using, how to get actual coefficient used?
             'zonals': [-Constants.IERS2010_EARTH_C20]}
    # Define spherical altitude for convience, not specifically force related, but uses the definitions
    global sphalt
    sphalt = OneAxisEllipsoid(force['earthrad'], 0.0, celestdflt['earthframe'])
    return celestdflt | force

# Set the default gravity 0x0
forcedflt = setgravity(0, 0)

# Create an example spacecraft with B = C_D A/m = 0.01 m^2/kg
# scB010 = {'mass': 100.0,  # The models need a spacecraft mass, unit kg.
#         'dragarea': 1.0, # Cross-sectional area perpendicular to atmosphere direction, m^2
#         'dragcoef': 1.0 # Coefficient of drag
#         }
# scB010['drag'] = IsotropicDrag(scB010['dragarea'], scB010['dragcoef'])
# scB010['atmdens'] = envct['hp']
# scB010['dragforce'] = DragForce(scB010['atmdens'], scB010['drag']);

# def atmdens(location, time, model = 'hp'):
#     return(envct[model].getDensity(time, location, envct['celestialframe']))
