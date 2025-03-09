# Geographic locations and observations

import astropy.coordinates as coord
from propdemo import *

############### Earth locations

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = tell.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcd = tell.hcat(demoa.prop.ephem, eloc.mcdsv) # Ephemeris table with additional column for McDonald site vector
eloc.mcdsvork = tork.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcdork = tell.hcat(demoa.prop.ephem, eloc.mcdsvork) # Ephemeris table with additional column for McDonald site vector

# Difference between AstroPy and Orekit
eloc.mcdsv_apy = tell.makepos(tell.eciobs(eloc.mcd, newyear))
eloc.mcdsv_ork = tell.makepos(tork.eciobs(eloc.mcd, newyear))
eloc.mcdsv_apy_ork_dist = np.linalg.norm(eloc.mcdsv_ork - eloc.mcdsv_apy).si

eloc.kickapoo = tell.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')
# From tork.location('carbarn') - integrate earthloc with location()?
eloc.carbarn = tell.earthloc(lat=38.87206*u.deg, lon=-77.01748*u.deg, elevation=13.0*u.m)
eloc.carbarn_lst_newyear = tell.siderealtime(newyear, eloc.carbarn)

############### Earth observations

eobs = Munch()

# A simulated observation from carbarn
eobs.ob = tell.azelrange(255*u.deg, 80*u.deg, 1455*u.km, eloc.carbarn, newyear)
# ECI Cartesian coordinates
eobs.obeci = eobs.ob.transform_to(coord.GCRS)
eobs.obeci_cart = eobs.obeci.cartesian
eobs.obeci_radec = eobs.obeci.spherical

# An observation from McDonald
eobs.obmcd = eobs.ob.transform_to(coord.AltAz(location=eloc.mcd, obstime=newyear))

# Angles only
# obao = coord.SkyCoord(coord.AltAz(az=63*u.deg, alt=80*u.deg, location=carbarn, obstime=newyear))
# obsloc_carbarn = coord.AltAz(location=carbarn, obstime=newyear)

eobs.ptobork = tork.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ptobapy = tell.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ob_ork_apy_dist = np.linalg.norm(eobs.ptobork.cartesian.xyz - eobs.ptobapy.cartesian.xyz).si
