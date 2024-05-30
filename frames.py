# Longitude, latitude, altitude
# Make an lla class that saves the ork as an attribute
# ex1lla = ex1.pvt.lla(ex1.forceenv)
def lla(position, time, forceenv):
    return forceenv['earth'].transform(position, forceenv['celestialframe'], time);
