import collections.abc

from org.orekit.attitudes import LofOffset
from org.orekit.frames import LOFType
from org.orekit.propagation import SpacecraftState

from tellurion.astro import units
from tellurion.ork import convert

def ntw(relsc, refsc, unitlookup=units.prefunits):
    '''Find the NTW (normal, in-track, cross-track) relative
    coordinates of relsc with respect to refsc; each is a
    org.orekit.propagation.SpacecraftState (set the `spacecraftstate`
    argument of `prop()` to `True`)
    '''
    rf = _relframe(LOFType.NTW, relsc, refsc, unitlookup)
    return rf.ephemeris(['time', 'NTW position', 'NTW velocity'])

def lvlh(relsc, refsc, unitlookup=units.prefunits):
    '''Find the LVLH or RSW (radial, along-track, cross-track) relative
    coordinates of relsc with respect to refsc; each is a
    org.orekit.propagation.SpacecraftState (set the `spacecraftstate`
    argument of `prop()` to `True`)
    '''
    rf = _relframe(LOFType.LVLH, relsc, refsc, unitlookup)
    return rf.ephemeris(['time', 'LVLH position', 'LVLH velocity'])

def _relframe(frame, relsc, refsc, unitlookup=units.prefunits):
    if isinstance(refsc, collections.abc.Iterable):
        rf = [_relframe(frame, rel, ref, unitlookup) for (rel, ref) in zip(relsc, refsc)]
        return rf[0].merge(rf[1:])
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
# [2025-09-12 Fri 13:42] .spherical method for PVT will do this, want to cyclic permute to WNT?
