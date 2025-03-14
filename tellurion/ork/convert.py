#from org.orekit.orbits import CartesianOrbit, OrbitType, Orbit

from ..core import astro
from ..core import element
from . import force
from . import posvel as oposvel
from . import element as oelement

# The following return AstroPy objects
# SpacecraftState.pvt = lambda self, unitlookup=astro.prefunits: oposvel._pvtork(self.pVCoordinates, unitlookup)
# BoundedPropagator.pvt = lambda self: self.initialState.pvt()
# SpacecraftState.kepler = lambda self: self.orbit.kepler()
# BoundedPropagator.kepler  = lambda self: self.initialState.kepler()
# The following return Orekit objects
# SpacecraftState.cartesianorbit = lambda self: CartesianOrbit.cast_(OrbitType.CARTESIAN.convertType(self))
# BoundedPropagator.cartesianorbit  = lambda self: self.initialState.cartesianorbit()
# SpacecraftState.keplerianorbit = lambda self: self.orbit.keplerianorbit()
# BoundedPropagator.keplerianorbit  = lambda self: self.initialState.keplerianorbit()

# The following return AstroPy objects
def pvt(object, unitlookup=astro.prefunits):
    if type(object) is Orbit:
        return oposvel._pvtork(object.pVCoordinates, unitlookup)
    elif type(object) is SpacecraftState:
        return oposvel._pvtork(object.pVCoordinates, unitlookup)
    elif type(object) is BoundedPropagator:
        return pvt(object.initialState)
    elif type(object) is TimeStampedPVCoordinates:
        return oposvel._pvtork(self, unitlookup)

def kepler(object, units=(astro.prefunits['length'], astro.prefunits['angle']), forceenv=force.deffe):
    if type(object) is Orbit:
        return oelement._kepler(object, units)
    elif type(object) is SpacecraftState:
        return kepler(object.orbit, units)
    elif type(object) is BoundedPropagator:
        return kepler(object.initialState, units)
    elif type(object) is TimeStampedPVCoordinates:
        return kepler(cartesianorbit(object, forceenv), forceenv, units)

# The following return Orekit objects
def keplerianorbit(object):
    if type(object) is Orbit:
        return OrbitType.KEPLERIAN.convertType(object)
    elif type(object) is SpacecraftState:
        return keplerianorbit(object.orbit)
    elif type(object) is BoundedPropagator:
        return keplerianorbit(object.initialState)

def cartesianorbit(object):
    if type(object) is Orbit:
        return OrbitType.CARTESIAN.convertType(object)
    elif type(object) is SpacecraftState:
        return cartesianorbit(object.orbit)
    elif type(object) is BoundedPropagator:
        return cartesianorbit(object.initialState)
