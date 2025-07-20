"""Propagate a state in two steps: `prepare()` and `propagate()`

Events (e.g. eclipse, visibility) should supply three functions
`_mkdetlog()` to add the event detector(s)/logger(s), `_gentrans` to
generate the event transitions, and `_statechar()` to generate a
character indicating each state.

"""

import collections
import collections.abc
import numpy as np
import astropy.units as u
import astropy.time
import astropy.table
from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit, KeplerianOrbit
from org.orekit.propagation.numerical import NumericalPropagator
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import Propagator, BoundedPropagator, SpacecraftState, EphemerisGenerator
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.utils import AbsolutePVCoordinates, TimeStampedPVCoordinates, PVCoordinatesProvider
import org.orekit.forces.gravity as okgrav

from ..core import astro
from ..core import posvel
from ..core.posvel import PVT
from ..core import element
from ..core.spacetrack import MeanElementSetT
from . import force
from . import element as oelement
from . import convert
from . import event

defev = {'altitude': 125.0*u.km, 'eclipse': [], 'visibility': []}

# import astropy.units as u
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# sentprep = tork.SGP4prep(isssent['SENTINEL 3A'], 1*u.day, {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []})
# tork.propagate(sentprep, np.linspace(5.0*u.minute, 60.0*u.minute, 12), True)

def prepare(initstate, proptime, events=defev, forceenv=force.deffe, reftime=None, occluder='earth', \
            output='et'):
    if type(initstate)==MeanElementSetT:
        return SGP4prep(initstate, proptime, events, forceenv, reftime, occluder)
    elif type(initstate)==PVT:
        return niprep(initstate, proptime, events, forceenv, reftime, occluder, output=output)

def SGP4prep(meanels, proptime, events=defev, forceenv=force.deffe, reftime=None, occluder='earth', \
             output='et'):
    '''Propagate mean elements using SGP4'''
    if hasattr(meanels, 'model') and meanels.model == 'SGP4':
        tle = TLE(*meanels.tle)
        propagator = TLEPropagator.selectExtrapolator(tle)
        gendict = _make_gendict(tle, lambda propto: \
                                propagator.getPVCoordinates(propto, forceenv['celestialframe']))
        gendict['pvt0'] = convert._pvt(gendict['propfn'](gendict['epoch']))
        event._add(events, propagator, gendict, proptime, forceenv, reftime, output)
        return gendict
    else:
        raise ValueError("Can only propagate SGP4 mean elements with SGP4")

def niprep(initstate, proptime, events=defev, forceenv=force.deffe, reftime=None, occluder='earth', \
           output='et'):
    """Make a generator for an ephemeris, optionally include eclipse
    information. The result of this function is passed as the first
    argument to `propagate()`.

    Parameters
    ----------
    initstate: a tell.ispvter() representing the initial state

    proptime:  u.Quantity, float
      The maximum time to propagate, numbers are in seconds

    forceenv:  dict
      Forces to use; output of force.setgravity()

    stopalt:   float
      Stop propagation if altitude above spherical earth (m) drops below this threshold

    Returns
    -------
    An Orekit object that is passed to `propagate()` as the first
    argument, if no event other than altitude is included in `events`.

    If eclipse detection is added to `events`, then a dictionary with
    the generator `['ephgen']`, and an ephemeris table `['sun
    transitions']` with xyz positions is returned. The latter has a
    two-character string 'suntrans', and elapsed time from the
    previous transition 'elapsed'. There are also two detectors used
    by `propagate()` in the dictionary.

    """

    ork0 = CartesianOrbit(convert._tspvc(initstate), \
                          forceenv['celestialframe'], forceenv['earthmu'].si.value)

    # Set parameters
    minstep = 0.001
    maxstep = 1000.0
    initStep = 60.0
    positionTolerance = 1.0
    tolerances = NumericalPropagator.tolerances(positionTolerance, ork0, ork0.getType())

    # Initialize the integrator
    integrator = DormandPrince853Integrator(
        minstep,
        maxstep,
        tolerances[0],
        tolerances[1])
    integrator.setInitialStepSize(initStep)

    initialState = SpacecraftState(ork0, forceenv['mass'])
    okprop = NumericalPropagator(integrator)
    okprop.setOrbitType(OrbitType.CARTESIAN)
    okprop.setInitialState(initialState)
    generator = okprop.getEphemerisGenerator()

    # Forces
    okprop.addForceModel(okgrav.HolmesFeatherstoneAttractionModel(forceenv['earthframe'], forceenv['gravity']))
    if 'dragforce' in forceenv:
        okprop.addForceModel(forceenv['dragforce'])

    # Set up event detectors and propagate
    gendict = _make_gendict(ork0, lambda propto: okprop.propagate(gendict['epoch'], propto))
    event._add(events, okprop, gendict, proptime, forceenv, reftime, output)
    gge = generator.getGeneratedEphemeris()
    gendict['mindate'] = gge.getMinDate()
    gendict['maxdate'] = gge.getMaxDate()

    return gendict

def _make_gendict(getdatefrom, propfn):
    gendict = {}
    gendict['event detectors'] = {}
    gendict['epoch'] = getdatefrom.getDate()
    gendict['propfn'] = propfn
    return gendict

def propagate(generator, reltimes, include_init=True, reftime='epoch', output='et'):
    '''From an existing ephemeris generator, propagate to the time(s)
    relative to epoch of the initial state. The relative times must
    satisfy posvel.isreltime(reltimes), and if the size
    prop5m1h.shape[0] > 0, an ephemeris table is returned. If reltimes
    is a single time, then a PVT is returned. If `include`_init is
    true, then include the initial PVT in the ephemeris table.

    Output is specified as one of
      `et`: Ephemeris table (an AstroPy time series; default)
      `pvt`: a PVT
      `ss`: SpacecraftState

    If an eclipse detector has been added, the `'sunlight'` value will
    be one of 'u' (umbra, or total eclipse), 'p' (penumbra, or partial
    eclipse, or 's' (full sun).
    '''

    # Compute the reference time `reft`, an astropy.time.Time
    if reftime=='epoch' and generator.get('epoch'):
        reft = convert._okad(generator.get('epoch'))
    elif reftime=='epoch' and generator.get('mindate'):
        reft = convert._okad(generator.get('mindate'))
    else:
        reft = astro.abstime(reftime)

    atimes = astro.abstime(reltimes, reft)
    if include_init:
        reltimes = np.insert(atimes, 0, reft)
    atscalar = atimes.shape == ()

    if atscalar:
        # This includes the value of the event function
        ss = _to_spacecraft_state(generator['propfn'](convert._okad(atimes)))
        if output=='ss':
            return ss
        else:
            return event._ephemeris(generator, ss)
    else:
        if output=='ss':
            datalist = [propagate(generator, rt, False, reftime, output) for rt in reltimes]
            return datalist
        else:
            retpvt = propagate(generator, reltimes[0], False, reftime, 'pvt')
            for rt in reltimes[1:]:
                retpvt.concatenate(propagate(generator, rt, False, reftime, 'pvt'))
            if output=='pvt':
                return retpvt
            else:
                return retpvt.ephemeris()

def timerange(object):
    '''The time difference between the earliest (usually the initial
    time) and the latest; not always what is requested as atmospheric
    drag can shorten the timespan
    '''
    if object.get('maxdate') and object.get('mindate'):
        return (convert._okad(object.get('maxdate'))-convert._okad(object.get('mindate'))).to(u.s)
    else:
        raise ValueError('Cannot get timerange for this object')

def _to_spacecraft_state(tspvc, forceenv=force.deffe):
    '''From the TimeStampedPVCoordinates, create a SpacecraftState'''
    if type(tspvc)==TimeStampedPVCoordinates:
        apvc = AbsolutePVCoordinates(forceenv['celestialframe'], tspvc)
        return SpacecraftState(apvc)
    else:
        return tspvc
