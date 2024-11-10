"""
# Geographic points specified by latitude, longitude, altitude/elevation

# Lookup geographic points by name; must set lookup.geonames_username first
# See https://gitlab.com/-/snippets/3743091 to get a username
## Example without time
In [8]: angellhall = geopt('angell hall')
In [10]: angellhall.llaq
Out[10]: <Quantity (-83.74001, 42.2771, 281.) (deg, deg, m)>
In [11]: angellhall.ork
Out[11]: <GeodeticPoint: {lat: 42.2771 deg, lon: -83.74001 deg, alt: 281}>
In [12]: angellhall.skf
Out[12]: <GeographicPosition WGS84 latitude +42.2771 N longitude -83.7400 E elevation 281.0 m>

## Add a time to it
angellhall
<Angell Hall, Michigan, United States; timeless>
angellhallnow = LLA(angellhall, dttm.nowutc())
<Angell Hall, Michigan, United States; epoch 2024-09-03 01:22:38.409828 (UTC)>

## Example with time
In [13]: mcdonald = geopt('mcdonald observatory', dttm.nowutc())
In [14]: mcdonald
Out[14]: <McDonald Observatory, Texas, United States; epoch 2024-09-02 23:14:39.293387 (UTC)>
In [15]: mcdonald.llaq
Out[15]: <TQuantity (-104.02158, 30.67154, 2050.) (deg, deg, m), time=2024-09-02 02:12:52.997075>
In [16]: mcdonald.ork
Out[16]: <GeodeticPoint: {lat: 30.67154 deg, lon: -104.02158 deg, alt: 2,050}>
In [17]: mcdonald.skf
Out[17]: <GeographicPosition WGS84 latitude +30.6715 N longitude -104.0216 E elevation 2050.0 m>

# The old NAVSPASUR main transmitter
# Coordinates from https://www.thelivingmoon.com/45jack_files/04images/Kickapoo/Kickapoo_001.png
# klat = Angle('33°33′08.50″N').value
# klon = Angle('98°45′49.82″W').value
# kickapoo = LLA([klat, klon, elev(klat, klon), 'Historic NAVSPASUR transmitter', 'Texas', 'US'], latlon=True, length_unit=u.m)
"""

# ex1lla = ex1.pvt.lla(ex1.forceenv)

import dttm
import lookup
import astropy.coordinates
import skyfield.api
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
    cart = topoframe.getPVCoordinates(dttm.to_okad(lla.llaq.time), forceenv['celestialframe']).getPosition()
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
# exg = LLA([-40.0, 0.1, 250.0], dttm.nowutc())
class LLA:
    llaq: TQuantity
    ork: GeodeticPoint

    def __init__(self, fromthing, time=None,
                 length_unit=prefunits["length"], angle_unit=prefunits["angle"], latlon=False):
        if type(fromthing) == GeodeticPoint or type(fromthing) == FieldGeodeticPoint:
            self.ork = fromthing
            lu = lonlatalt([fromthing.longitude, fromthing.latitude, fromthing.altitude],
                           time, length_unit=u.meter, angle_unit=u.radian)
            self.llaq = lu.convert_units((angle_unit,angle_unit,length_unit))
            # elif type(fromthing) == astropy.table.row.Row:
            # elif type(fromthing) == Orbit:
        elif type(fromthing) == LLA:
            self.llaq = tquant(quant(fromthing.llaq), time)
            self.ork = fromthing.ork
            self.skf = fromthing.skf
            if hasattr(fromthing,"info"):
                self.info = fromthing.info
        else:
            if latlon:
                q = latlonalt(fromthing, time=time, length_unit=length_unit, angle_unit=angle_unit)
            else:
                q = lonlatalt(fromthing, time=time, length_unit=length_unit, angle_unit=angle_unit)
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
                return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()}) epoch {self.llaq.time} (UTC)>"
        else:
            if hasattr(self,"info"):
                return f"<{self.info}; timeless>"
            else:
                return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()})>"

nullisland = LLA([0.0, 0.0, 0.0])

# Lookup geographic points by name; must set lookup.geonames_username first
def geopt(name, time=None):
    return LLA(lookup.location(name), time, length_unit=u.m, latlon=True)

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
        tm= dttm.nowutc()
    return astropy.coordinates.Angle(lla.skf.lst_hours_at(dttm.to_skftime(tm)), u.hour).to(u.deg)

# Time series of lla
#def tslla(ephem, elements):
#    return TimeSeries(time=[], data=[])
