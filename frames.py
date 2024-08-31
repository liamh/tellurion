# Longitude, latitude, altitude
# Make an lla class that saves the ork as an attribute
# ex1lla = ex1.pvt.lla(ex1.forceenv)

from astro import *
from dttm import *
from astropy.coordinates import Angle
from skyfield.api import wgs84
from org.orekit.bodies import GeodeticPoint, FieldGeodeticPoint
from org.orekit.frames import TopocentricFrame

# Find the geographic (LLA) coordinates from the ECI position and time
# Returns a FieldGeodeticPoint, should return an LLA
# llafrompt(ex1.pvt.ork, ex1.pvt.ork.date, ex1.forceenv)
def llafrompt(position, time, forceenv):
    return forceenv['earth'].transform(position, forceenv['celestialframe'], time);

# Find the ECI position from the geographic (LLA) coordinates
# Returns a Vector3D, should return a PVT
# ptfromlla(ex1.pvt.lla(ex1.forceenv),ex1.forceenv)
def ptfromlla(lla, forceenv):
    topoframe = TopocentricFrame(forceenv['earth'], lla.ork, "ptfromlla")
    cart = topoframe.getPVCoordinates(to_okad(lla.llaq.time), forceenv['celestialframe']).getPosition()
    return cart

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

# Class of gegraphic points (longitude, latitude, altitude)
# exg = LLA([-40.0, 0.1, 250.0], nowutc())
class LLA:
    llaq: TQuantity
    ork: GeodeticPoint

    def __init__(self, fromthing, time=None,
                 length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
        if type(fromthing) == GeodeticPoint or type(fromthing) == FieldGeodeticPoint:
            self.ork = fromthing
            lu = lonlatalt([fromthing.longitude, fromthing.latitude, fromthing.altitude],
                           time, length_unit=u.meter, angle_unit=u.radian)
            self.llaq = lu.convert_units((angle_unit,angle_unit,length_unit))
            # elif type(fromthing) == astropy.table.row.Row:
            # elif type(fromthing) == Orbit:
        else:
            q = lonlatalt(fromthing, time=time, length_unit=length_unit, angle_unit=angle_unit)
            self.llaq = q
            self.ork = GeodeticPoint(float(q["lat"].to(u.radian).value), float(q["lon"].to(u.radian).value),
                                     float(q["alt"].to(u.meter).value))
        q = self.llaq
        self.skf = wgs84.latlon(q["lat"].to(u.deg).value, q["lon"].to(u.deg).value, q["alt"].to(u.meter).value)
    def __repr__(self):
        if hasattr(self,"time"):
            return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()}) epoch {self.llaq.time} (UTC)>"
        else:
            return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()})>"

nullisland = LLA([0.0, 0.0, 0.0])

# Find the local sidereal time for the location
# If time is supplied explicitly, it is used, if it's not, then the LLA's time is used.
# If that doesn't exist, the current time is used.
# If the LLA is not supplied or given as None, then the Greenwich
#  sidereal time (GST) is returned.
# siderealtime(exg)
def siderealtime(lla=None, time=None):
    if lla is None:
       lla = nullisland
    if type(time)==Time:
        tm = time
    elif hasattr(lla.llaq, 'time'):
        tm = lla.llaq.time
    else:
        tm= nowutc()
    return Angle(lla.skf.lst_hours_at(to_skftime(tm)), u.hour).to(u.deg)

# Time series of lla
#def tslla(ephem, elements):
#    return TimeSeries(time=[], data=[])
