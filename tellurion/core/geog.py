"""Computation of GeodeticPoint and ECI computation of site vectors
and angles & range observation eciobs()

"""

# Site vector calculation
import astropy.coordinates as coord
import astropy.units as u

from tellurion.astro import units
from tellurion.core import geonames


def observer_location(location, name="Unnamed", minimum_elevation=coord.Angle(10.0, u.deg)):
    """Define an observer's location; location is an
    `astropy.coordinates.earth.EarthLocation`, e.g. the output of
    `earthloc()`. The `minimum_elevation` is used for visibility
    calculation; if a number is specified, angle preferred units
    `prefunits['angle']` are assumed.
    """
    if type(minimum_elevation) is coord.Angle:
        me = minimum_elevation
    else:
        me = coord.Angle(minimum_elevation, units.prefunits["angle"])
    return {"location": location, "name": name, "minelev": me}

# List of sites: astropy.coordinates.EarthLocation.get_site_names()

def earthloc(lon, lat, elevation=None):
    """Make an earth location; if elev=None (default), look it up."""
    longi = coord.Angle(lon)
    latit = coord.Angle(lat)
    if elevation is None:
        try:
            elevation = geonames.elev(latit.value, longi.value)*u.m
        except:
            elevation = 0.0*u.m
    return coord.EarthLocation.from_geodetic(lon=longi, lat=latit, height=elevation.si)

nullisland = earthloc(0*u.deg, 0*u.deg, 0*u.m)
