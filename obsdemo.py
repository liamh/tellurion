# Geographic locations

from propdemo import *

import astropy.coordinates as coord
import tellurion.geonames as geonames
import tellurion.ork.geog as ogeog
import tellurion.geog as geog

############### Earth locations

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = geog.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcd = astro.hcat(demoa.prop.ephem, eloc.mcdsv) # Ephemeris table with additional column for McDonald site vector
eloc.mcdsvork = ogeog.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcdork = astro.hcat(demoa.prop.ephem, eloc.mcdsvork) # Ephemeris table with additional column for McDonald site vector

# Difference between AstroPy and Orekit
eloc.mcdsv_apy = posvel.makepos(geog.eciobs(eloc.mcd, newyear))
eloc.mcdsv_ork = posvel.makepos(ogeog.eciobs(eloc.mcd, newyear))
eloc.mcdsv_apy_ork_dist = np.linalg.norm(eloc.mcdsv_ork - eloc.mcdsv_apy).si

eloc.kickapoo = geog.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')
# From geonames.location('carbarn') - integrate earthloc with location()?
eloc.carbarn = geog.earthloc(lat=38.87206*u.deg, lon=-77.01748*u.deg, elevation=13.0*u.m)
eloc.carbarn_lst_newyear = geog.siderealtime(newyear, eloc.carbarn)

############### Earth observations

eobs = Munch()

# A simulated observation from carbarn
eobs.ob = geog.azelrange(255*u.deg, 80*u.deg, 1455*u.km, eloc.carbarn, newyear)
# ECI Cartesian coordinates
eobs.obeci = eobs.ob.transform_to(coord.GCRS)
eobs.obeci_cart = eobs.obeci.cartesian
eobs.obeci_radec = eobs.obeci.spherical

# An observation from McDonald
eobs.obmcd = eobs.ob.transform_to(coord.AltAz(location=eloc.mcd, obstime=newyear))

# Angles only
# obao = coord.SkyCoord(coord.AltAz(az=63*u.deg, alt=80*u.deg, location=carbarn, obstime=newyear))
# obsloc_carbarn = coord.AltAz(location=carbarn, obstime=newyear)

eobs.ptobork = ogeog.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ptobapy = geog.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ob_ork_apy_dist = np.linalg.norm(eobs.ptobork.cartesian.xyz - eobs.ptobapy.cartesian.xyz).si
