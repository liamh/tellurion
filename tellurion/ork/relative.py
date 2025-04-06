import collections.abc

from org.orekit.attitudes import LofOffset
from org.orekit.frames import LOFType
from org.orekit.propagation import SpacecraftState

from ..core import astro
from ..core import posvel
from . import convert


def ntw(relsc, refsc, unitlookup=astro.prefunits):
    '''Find the NTW (normal, in-track, cross-track) relative
    coordinates of relsc with respect to refsc; each is a
    org.orekit.propagation.SpacecraftState (set the `spacecraftstate`
    argument of `prop()` to `True`)
    '''
    rf =  _relframe(LOFType.NTW, relsc, refsc, unitlookup)
    return posvel.tsephem(rf, None, ['time', 'NTW position', 'NTW velocity'])

def lvlh(relsc, refsc, unitlookup=astro.prefunits):
    '''Find the LVLH or RSW (radial, along-track, cross-track) relative
    coordinates of relsc with respect to refsc; each is a
    org.orekit.propagation.SpacecraftState (set the `spacecraftstate`
    argument of `prop()` to `True`)
    '''
    rf = _relframe(LOFType.LVLH, relsc, refsc, unitlookup)
    return posvel.tsephem(rf, None, ['time', 'LVLH position', 'LVLH velocity'])

def _relframe(frame, relsc, refsc, unitlookup=astro.prefunits):
    if isinstance(refsc, collections.abc.Iterable):
        return [_relframe(frame, rel, ref, unitlookup) for (rel, ref) in zip(relsc, refsc)]
    # See https://forum.orekit.org/t/1478/2
    relframe = LofOffset(refsc.getFrame(), frame);
    converted = SpacecraftState(refsc.getOrbit(), \
                                relframe.getAttitude(refsc.getOrbit(), \
                                                refsc.getDate(), refsc.getFrame()));
    pvrel = converted.toTransform().transformPVCoordinates(relsc.getPVCoordinates(refsc.getFrame()));
    return convert._pvt(pvrel, unitlookup)

# import astropy.coordinates as coord
# cart12 = coord.CartesianRepresentation(demoa.prop.ntw_4x4hpB01_to_4x4[12]['NTW position'])
# coord.SphericalRepresentation.from_cartesian(cart12)
