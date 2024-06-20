# Longitude, latitude, altitude
# Make an lla class that saves the ork as an attribute
# ex1lla = ex1.pvt.lla(ex1.forceenv)

from astro import *
from org.orekit.bodies import GeodeticPoint

def llafrompt(position, time, forceenv):
    return forceenv['earth'].transform(position, forceenv['celestialframe'], time);

# Class of gegraphic points (longitude, latitude, altitude)
# LLA([-40.0, 0.1, 250.0], nowutc())
class LLA:
    llaq: TQuantity
    ork: GeodeticPoint

    def __init__(self, fromthing, time=None,
                 length_unit=prefunits["length"], angle_unit=prefunits["angle"]):
        if type(fromthing) == GeodeticPoint:
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
    def __repr__(self):
        if hasattr(self,"time"):
            return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()}) epoch {self.llaq.time} (UTC)>"
        else:
            return f"<LLA longitude: {self.llaq['lon'].value} ({self.llaq.unit[0].to_string()}) latitude: {self.llaq['lat'].value} ({self.llaq.unit[1].to_string()}) altitude: {self.llaq['alt'].value} ({self.llaq.unit[2].to_string()})>"
