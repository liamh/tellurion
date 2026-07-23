import astropy.units as u

#from .
from org.orekit.utils import Constants
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.drag import AbstractDragForceModel, DragForce
from org.orekit.models.earth.atmosphere import Atmosphere, HarrisPriester, DTM2000, NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.forces.drag import IsotropicDrag

from org.orekit.frames import FramesFactory
from org.orekit.time import TimeScalesFactory
from org.orekit.utils import IERSConventions
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory

OneAxisEllipsoid.__repr__ =    lambda self: f"<Near-spherical body equatorial radius {self.getEquatorialRadius()}m, polar radius difference {-self.getEquatorialRadius()*self.getFlattening()}m >"

def celestial(force):
    efr = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    efl = Constants.IERS2010_EARTH_FLATTENING
    return {'earthframe': efr, \
            'celestialframe': FramesFactory.getEME2000(), \
            'earthangspd': Constants.IERS2010_EARTH_ANGULAR_VELOCITY, \
            'earthflat': efl, \
            'sun': CelestialBodyFactory.getSun(), \
            'sunrad': Constants.SUN_RADIUS*u.m, \
            'gmst': IERSConventions.IERS_2010.getGMSTFunction(TimeScalesFactory.getUT1(IERSConventions.IERS_2010, True)),
            'earth': OneAxisEllipsoid(force['earthrad'].si.value, efl,  efr),
            'sphalt': OneAxisEllipsoid(force['earthrad'].si.value, 0.0, efr),
            'moon': OneAxisEllipsoid(Constants.MOON_EQUATORIAL_RADIUS, 0.0012, \
                                     CelestialBodyFactory.getMoon().getBodyOrientedFrame())}


def setgravity(degree, order, mass = 100.0):
    """Set the environmental constants such as reference frame and
    planetary properties.  Arguments are the degree and order of the
    gravitational model to use, and the mass of the spacecraft.
    """
    eg = GravityFieldFactory.getNormalizedProvider(degree, order)
    force = {'gravity': eg,
             'gravity-unnorm':  # For Brouwer-Lyddane
             GravityFieldFactory.getUnnormalizedProvider(degree, order),
             'gravity-degree-order': [eg.getMaxDegree(), eg.getMaxOrder()],
             'earthrad': eg.getAe()*u.m,
             'earthmu': eg.getMu()*u.m**3/u.s**2,
             'mass': mass}  # Needed by several forces but not gravity
    return celestial(force) | force


def kepleranalytic():
    eg = GravityFieldFactory.getNormalizedProvider(0,0)
    force = {'earthrad': eg.getAe()*u.m,
             'earthmu': eg.getMu()*u.m**3/u.s**2}
    return  celestial(force) | force

# Default force & environment
deffe = setgravity(0,0)

_swdata = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

# This will find the unnormalized coefficients independent of any simulation.
# Returns three arrays: zonals (J2,...), Cnm, Snm
def unnormcoef(degree, order, when = 0.0*u.s):
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

default_atmdens = 'hp' # Harris-Priester
def dragforce(force, atmdensname = default_atmdens, dragcoef = 1.0, dragarea = 1.0):
    """Add atmospheric drag force to the force model. Default
    spacecraft parameters gives B = C_D A/m = 0.01 m^2/kg
    `atmdens` is one of `'hp'` for Harris-Priester (default), `'dtm'` for DTM2000, or `'msis'` for NRLMSIS.

    The parameters `dragcoef` and `drarea` should be either a
    numerical value, or a list of numerical value and a boolean. The
    boolean specifies whether the Jacobian matrix should be computed;
    if only a number is given, then the Jacobian will not be computed
    for that parameter.
    """
    ### There are no units given for dragarea; they should be SI
    match atmdensname:
        case 'hp': # Harris-Priester atmospheric density model
            atmdens = HarrisPriester(force['sun'], force['earth'])
        case 'dtm': # DTM2000 atmospheric density model
            atmdens = DTM2000(_swdata, force['sun'], force['earth'])
        case 'msis':
            atmdens = NRLMSISE00(_swdata, force['sun'], force['earth'])
    mass = force['mass']  # The models need a spacecraft mass, unit kg.
    if type(dragarea) is list:
        da = dragarea[0]
        dajac = dragarea[1]
    else:
        da = dragarea
        dajac = False
    if type(dragcoef) is list:
        dc = dragcoef[0]
        dcjac = dragcoef[1]
    else:
        dc = dragcoef
        dcjac = False
    df = DragForce(atmdens, IsotropicDrag(da, dc))
    df.getParametersDrivers()[0].setSelected(dajac)
    df.getParametersDrivers()[1].setSelected(dcjac)
    drforce = {'dragarea': da, # Cross-sectional area perpendicular to atmosphere direction, m^2
               'dragcoef': dc, # Coefficient of drag
               'B': da*dc/mass,
               'atmdens': atmdens,
               'dragforce': df}
    return force | drforce

def compute_drag_pjac(forceenv):
    """Boolean indicating whether the atmospheric drag parameters Jacobian should be calculated."""
    return forceenv.get('dragforce') \
        and (forceenv['dragforce'].getParametersDrivers()[0].isSelected() or
             forceenv['dragforce'].getParametersDrivers()[1].isSelected())




# def atmdens(location, time, model = 'hp'):
#     return(envct[model].getDensity(time, location, envct['celestialframe']))
