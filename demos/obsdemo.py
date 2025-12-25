# Geographic locations and observations

import astropy.coordinates as coord
from demos.propdemo import *

############### Earth locations

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = tork.sitevec(eloc.mcd, demoa.prop.ephem.time, "McDonald")

# Name doesn't transfer, too many columns
eloc.cmcd = tell.hcat(demoa.prop.ephem, eloc.mcdsv.ephemeris()) # Ephemeris table with additional column for McDonald site vector

# Fails after this point
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

############### Satellite visibility from ground observers

demod = Munch() # Example with SGP4 propagation
# tell.mestrt(tell.spacetrack_latest(stclient, 41335)) # So that we can reproduce the results
demod.sentinel3A = tell.makemest(((7180.806, 0.0001433, 98.6245, 260.8032, 90.6484, 269.4861,
                                    14.26733945, 6.4e-07, 0.0, 100.93, 801.642, 803.7, 0.0005667336263999999),
                                   ('sma', 'ecc', 'inc', 'raan', 'argper', 'ma', 'memo', 'memod', 'memodd',
                                    'period', 'peralt', 'apoalt', 'B'),
                                   ('km', '', 'deg', 'deg', 'deg', 'deg', 'revolution / d',
                                    'revolution / d2', 'revolution / d3', 'min', 'km', 'km', 'm2 / kg')),
                                  '2025-07-12 13:51:37.099',
                                  ('1 41335U 16011A   25193.57751272  .00000064  00000-0  44479-4 0  9993',
                                   '2 41335  98.6245 260.8032 0001433  90.6484 269.4861 14.26733945489635'),
                                  'SGP4',
                                  {'name': 'SENTINEL 3A',
                                   'type': 'PAYLOAD',
                                   'catid': 41335,
                                   'intldes': '2016-011A'})
demod.mcd_cb_vis = Munch() # Visibility from McDonald and Carbarn
demod.mcd_cb_vis.ev = {'altitude': 125.0*u.km, 'eclipse': [], \
                       'visibility': [tell.observer_location(eloc.mcd, "McDonald"), \
                                      tell.observer_location(eloc.carbarn, "Carbarn", 30.0)]}
demod.mcd_cb_vis.prep = tork.prepare(demod.sentinel3A, 1*u.day, demod.mcd_cb_vis.ev)
demod.mcd_cb_vis.mcdonald = demod.mcd_cb_vis.prep['visibility']['McDonald'] # McDonald visibility transitions
demod.mcd_cb_vis.carbarn = demod.mcd_cb_vis.prep['visibility']['Carbarn'] # Carbarn visibility transitions
demod.mcd_cb_vis.ephempvt = tork.propagate(demod.mcd_cb_vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32), output='pvt')
demod.mcd_cb_vis.ephem = tork.propagate(demod.mcd_cb_vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))
demod.mcd_ecl = Munch()  # Visibility from McDonald and eclipses
demod.mcd_ecl.ev = {'altitude': 125.0*u.km, 'eclipse': [True, True], \
                    'visibility': [tell.observer_location(eloc.mcd, "McDonald")]}
demod.mcd_ecl.prep = tork.prepare(demod.sentinel3A, 1*u.day, demod.mcd_ecl.ev)
demod.mcd_ecl.mcdonald = demod.mcd_ecl.prep['visibility']['McDonald'] # McDonald visibility transitions
demod.mcd_ecl.suntrans = demod.mcd_ecl.prep['sun transition']
demod.mcd_ecl.ephem = tork.propagate(demod.mcd_ecl.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))

############ Upward transitions not present

demoa.prop.vis = Munch() # Example with numerical propagation defined in demoa
demoa.prop.vis.events = {'altitude': 125.0*u.km, 'eclipse': [], \
                         'visibility': [tell.observer_location(eloc.mcd, "McDonald")]}
demoa.prop.vis.prep = tork.prepare(demoa.init.pvt, 8*u.hour, demoa.prop.vis.events)
demoa.prop.vis.ephem = tork.propagate(demoa.prop.vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))
demoa.prop.vis.mcdonald = demoa.prop.vis.prep['visibility']['McDonald']
