import orkinit
import dttm

from org.orekit.utils import Constants
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

def setgravity(degree, order, mass = 100.0):
    # Default celestial environment constants
    celestdflt = {
        'earthframe': FramesFactory.getITRF(IERSConventions.IERS_2010, True),
        'celestialframe': FramesFactory.getEME2000(),
        'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY,
        'earthflat': Constants.IERS2010_EARTH_FLATTENING,
        'sun': CelestialBodyFactory.getSun(),
    }
    eg = GravityFieldFactory.getNormalizedProvider(degree, order)
    force = {'gravity': eg,
             'gravity-degree-order': [eg.maxDegree, eg.maxOrder],
             'earthrad': eg.ae,
             'earthmu': eg.mu,
             'mass': mass} # Needed by several forces but not gravity
    # Define spherical altitude for convience, not specifically force related, but uses the definitions
    return celestdflt | force | \
        {'earth': OneAxisEllipsoid(force['earthrad'],
                                   celestdflt['earthflat'],  celestdflt['earthframe']),
         'sphalt': OneAxisEllipsoid(force['earthrad'], 0.0, celestdflt['earthframe'])}

swdata = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

# This will find the unnormalized coefficients independent of any simulation.
# Returns three arrays: zonals (J2,...), Cnm, Snm
def unnormcoef(degree, order, when = dttm.nowutc(True)):
    if type(when) is Quantity and get_physical_type(when) == 'time':
        when = dttm.to_okad(dttm.nowutc() + when)
    provider = GravityFieldFactory.getUnnormalizedProvider(degree, order)
    deg = provider.maxDegree
    ord = provider.maxOrder
    if deg > 0:
        if ord > 0:
            cnm = np.zeros((deg-1, ord))
            with np.nditer(cnm, flags=['multi_index'], op_flags=['writeonly']) as it:
                for x in it:
                    n = it.multi_index[0]+2
                    m = it.multi_index[1]
                    if n >= m:
                        x[...] = provider.onDate(when).getUnnormalizedCnm(n, m)
            snm = np.zeros((deg-1, ord))
            with np.nditer(snm, flags=['multi_index'], op_flags=['writeonly']) as it:
                for x in it:
                    n = it.multi_index[0]+2
                    m = it.multi_index[1]
                    if n >= m:
                        x[...] = provider.onDate(when).getUnnormalizedSnm(n, m)
        else:
            cnm = np.zeros(deg-1)
            with np.nditer(cnm, flags=['multi_index'], op_flags=['writeonly']) as it:
                for x in it:
                    x[...] = provider.onDate(when).getUnnormalizedCnm(it.multi_index[0]+2, 0)
            snm = np.zeros(deg-1)
            with np.nditer(snm, flags=['multi_index'], op_flags=['writeonly']) as it:
                for x in it:
                    x[...] = provider.onDate(when).getUnnormalizedSnm(it.multi_index[0]+2, 0)
        znl = -cnm[:,0]
    else:
        znl = None
        cnm = None
        snm = None
    return (znl, cnm, snm)

# Optionally add a drag force
# Default spacecraft parameters gives B = C_D A/m = 0.01 m^2/kg
default_atmdens = 'hp'
def dragforce(force, atmdensname = default_atmdens, dragcoef = 1.0, dragarea = 1.0):
    match atmdensname:
        case 'hp': # Harris-Priester atmospheric density model
            atmdens = HarrisPriester(force['sun'], force['earth'])
        case 'dtm': # DTM2000 atmospheric density model
            atmdens = DTM2000(swdata, force['sun'], force['earth'])
        case 'msis':
            atmdens = NRLMSISE00(swdata, force['sun'], force['earth'])
    mass = force['mass']  # The models need a spacecraft mass, unit kg.
    drforce = {'dragarea': dragarea, # Cross-sectional area perpendicular to atmosphere direction, m^2
               'dragcoef': dragcoef, # Coefficient of drag
               'B': dragarea*dragcoef/mass,
               'atmdens': atmdens,
               'dragforce': DragForce(atmdens, IsotropicDrag(dragarea, dragcoef))}
    return force | drforce

# def atmdens(location, time, model = 'hp'):
#     return(envct[model].getDensity(time, location, envct['celestialframe']))
