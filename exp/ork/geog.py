import astropy.units as u

from org.orekit.bodies import GeodeticPoint
from org.orekit.frames import TopocentricFrame
import org.orekit.models.earth as oearth

import geog
import posvel
import cdttm
import ork.force as ofr

def geodpt(earthloc):
    '''Create the Orekit GeodeticPoint from an AstroPy EarthLocation'''
    geod = earthloc.geodetic
    lat_rdn = float(u.Quantity(geod.lat).to(u.radian).value)
    lon_rdn = float(u.Quantity(geod.lon).to(u.radian).value)
    altitude_m = float(u.Quantity(geod.height).to(u.m).value)
    return GeodeticPoint(lat_rdn, lon_rdn, altitude_m)

def sitevec(loc, dttm, forceenv=ofr.deffe):
    '''Find the site vector (in the GCRS frame) of the EarthLocation using Orekit'''
    topoframe = TopocentricFrame(forceenv['earth'], geodpt(loc), "the location")
    posv3d = topoframe.getPVCoordinates(cdttm.okad(dttm), forceenv['celestialframe']).getPosition()
    pos = posvel.makepos(posv3d.quant(u.m))
    return (pos, dttm)
