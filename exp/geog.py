# Site vector calculation
import astropy.coordinates as coord
import cdttm
import posvel

# List of sites: astropy.coordinates.EarthLocation.get_site_names()
def sitevec(loc, dttm=None):
    '''Find the site vector (in the GCRS frame) of the EarthLocation using AstroPy'''
    if dttm==None:
        dttm = posvel.nowutc()
    itrs = coord.ITRS(coord.CartesianRepresentation(x=loc.x, y=loc.y, z=loc.z), obstime=dttm)
    gcrs = itrs.transform_to(coord.GCRS(obstime=dttm)).cartesian
    return (posvel.makepos(gcrs), dttm)
