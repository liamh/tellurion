"""
Define observation sites

 from astrodynamics import observersite

Defines: mcdonald, angellhall, kickapoo

To get different formats,
 .llaq AstrodynamicsPy
 .ork  Orekit
 .skf  Skyfield
General information is in .info.

# Add time
 from astrodynamics import dttm
 from astrodynamics import frames
 mcdonaldnow = frames.LLA(observersite.mcdonald, time=dttm.nowutc())

"""


from . import frames
from . import geonames
import astropy.coordinates as apc
import astropy.units as u



# Optical
mcdonald = frames.geopt('mcdonald observatory')
angellhall = frames.geopt('angell hall')


# Historic NAVSPASUR main transmitter
klat = apc.Angle('33°33′08.50″N').value
klon = apc.Angle('98°45′49.82″W').value
kickapoo = frames.LLA([klat, klon, geonames.elev(klat, klon), 'Historic NAVSPASUR transmitter', 'Texas', 'US'], latlon=True, length_unit=u.m)
