"""
Geographic points specified by latitude, longitude, altitude/elevation
Lookup geographic points by name; must set geonames.geonames_userid first
See https://gitlab.com/-/snippets/3743091 to get a username
For examples, see observersite.py
"""

# ex1lla = ex1.pvt.lla(ex1.forceenv)

from . import dttm
from . import geonames
from . import astro
from . import pvork  # Vectors and PVT (position, velocity, time) sets in Orekit
from . import orbit
from . import force
import astropy.units as u
import astropy.time
import astropy.timeseries as apts
import astropy.table as aptbl
import skyfield.api
from org.orekit.bodies import GeodeticPoint, FieldGeodeticPoint
from org.orekit.frames import TopocentricFrame
from org.orekit.utils import TimeStampedPVCoordinates


"""
Class of gegraphic points (longitude, latitude, altitude) with or without a timestamp

To make without a timestamp, call with `time` argument empty or
`None`. To make with a timestamp, give the time as the `time`
argument, or `True` to make it now.

import astrodynamics.frames as frames
import astrodynamics.observersite as obsite
kickapoonow = frames.LLA(obsite.kickapoo, True)
# Convert to ECI
kneci = kickapoonow.eci() # Gives the Kickapoo site vector at the current time
# Convert back to LLA
knecilla = kneci.lla().llaq
# Same as the original
knlla = kickapoonow.llaq
"""
class LLA:
    llaq: astro.TQuantity
    ork: GeodeticPoint

    def __init__(self, fromthing, time=None,
                 length_unit=astro.prefunits["length"], angle_unit=astro.prefunits["angle"], latlon=False):
        if time==True:
            time = dttm.nowutc()
        if type(fromthing) == GeodeticPoint or type(fromthing) == FieldGeodeticPoint:
            self.ork = fromthing
            lu = astro.lonlatalt([fromthing.longitude, fromthing.latitude, fromthing.altitude],
                           time, length_unit=u.meter, angle_unit=u.radian)
            self.llaq = lu.convert_units((angle_unit,angle_unit,length_unit))
            # elif type(fromthing) == astropy.table.row.Row:
            # elif type(fromthing) == Orbit:
        elif type(fromthing) == LLA:
            self.llaq = astro.tquant(astro.quant(fromthing.llaq), time)
            self.ork = fromthing.ork
            self.skf = fromthing.skf
            if hasattr(fromthing,"info"):
                self.info = fromthing.info
        elif type(fromthing) == TimeStampedPVCoordinates:
            self.__init__(fromthing.getPosition(), fromthing.getDate())
        else:
            if latlon:
                q = astro.latlonalt(fromthing, time=time, length_unit=length_unit, angle_unit=angle_unit)
            else:
                q = astro.lonlatalt(fromthing, time=time, length_unit=length_unit, angle_unit=angle_unit)
            self.llaq = q
            self.ork = GeodeticPoint(float(q["lat"].to(u.radian).value), float(q["lon"].to(u.radian).value),
                                     float(q["alt"].to(u.meter).value))
            if len(fromthing) > 3:
                self.info = ', '.join(fromthing[3:]) # other information, such as name and location
        q = self.llaq
        self.skf = skyfield.api.wgs84.latlon(q["lat"].to(u.deg).value, q["lon"].to(u.deg).value, q["alt"].to(u.meter).value)
    def __repr__(self):
        if hasattr(self.llaq,"time"):
            if hasattr(self,"info"):
                return f"<{self.info}; epoch {self.llaq.time} (UTC)>"
            else:
                return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) " \
                    f"latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) " \
                    f"altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()}) epoch {self.llaq.time} (UTC)>"
        else:
            if hasattr(self,"info"):
                return f"<{self.info}; timeless>"
            else:
                return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) " \
                    f"latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) " \
                    f"altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()})>"
    def eci(self, forceenv=force.setgravity(0,0)):
        '''
        Find the ECI position from the geographic (LLA) coordinates, return a PVT
        '''
        topoframe = TopocentricFrame(forceenv['earth'], self.ork, "ptfromlla")
        pos = topoframe.getPVCoordinates(dttm.to_okad(self.llaq.time), forceenv['celestialframe']).getPosition()
        orkpt = pvork.orkpvt(pos, pvork.v3dnan, self.llaq.time)
        return orbit.PVT(orkpt)

def llafrompt(position, time, forceenv):
    '''Find the geographic (LLA) coordinates from the ECI position and time'''
    return LLA(forceenv['earth'].transform(position, forceenv['celestialframe'], time), time);

# Transform from ECI to ECEF
# eci_to_ecef(ex1.pvt, ex1.forceenv)
def eci_to_ecef(pvt, forceenv):
    xf = forceenv['celestialframe'].getTransformTo(forceenv['earthframe'], pvt.ork.date)
    ecef = xf.transformPVCoordinates(pvt.ork)
    return ecef

# Not tested
def ecef_to_eci(ecef, forceenv):
    xf = forceenv['earthframe'].getTransformTo(forceenv['celestialframe'],ecef.ork.date)
    return xf.transformPVCoordinates(ecef.ork)

nullisland = LLA([0.0, 0.0, 0.0])

# Lookup geographic points by name; must set geonames.geonames_userid first
def geopt(name, time=None):
    return LLA(geonames.location(name), time, length_unit=u.m, latlon=True)

# Find the local sidereal time for the location
# If time is supplied explicitly, it is used, if it's not, then the LLA's time is used.
# If that doesn't exist, the current time is used.
# If the LLA is not supplied or given as None, then the Greenwich
#  sidereal time (GST) is returned.
# siderealtime(exg)
def siderealtime(lla=None, time=None):
    if lla is None:
       lla = nullisland
    if type(time)==astropy.time.Time:
        tm = time
    elif hasattr(lla.llaq, 'time'):
        tm = lla.llaq.time
    else:
        tm= dttm.nowutc()
    return astropy.coordinates.Angle(lla.skf.lst_hours_at(dttm.to_skftime(tm)), u.hour).to(u.deg)

# Time series of lla
def tslla(ephem, elements):
    return apts.TimeSeries(time=[], data=[])

## Time series of orbital elements
def tselements(ephem, elements):
    return apts.TimeSeries(time=[pvt.pvtq.time for pvt in ephem.pvt],
                      data=[dict(zip(elements, ))
                            for orb in ephem.orbit])
