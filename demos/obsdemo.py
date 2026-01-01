# Geographic locations and observations

import astropy.coordinates as coord
from demos.propdemo import *

############### Earth locations

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = tork.sitevec(eloc.mcd, demoa.propn.ephem.time, "McDonald")
eloc.mcdsvts = eloc.mcdsv.ephemeris(False, None)   # Time series of site vectors for McDonald
eloc.cmcd = tell.hcat(demoa.propn.ephem, eloc.mcdsvts) # Ephemeris table with additional columns for McDonald site vector

eloc.kickapoo = tell.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')
# From tork.location('carbarn') - integrate earthloc with location()?
eloc.carbarn = tell.earthloc(lat=38.87206*u.deg, lon=-77.01748*u.deg, elevation=13.0*u.m)
eloc.carbarn_lst_newyear = tork.siderealtime(newyear, eloc.carbarn)

############### Earth observations

eobs = Munch()

# A simulated observation from carbarn
eobs.ob = tell.azelrange(255*u.deg, 80*u.deg, 1455*u.km, eloc.carbarn, newyear)
eobs.obpt = tork.eciaer(eobs.ob)
eobs.obret = tork.aereci(eobs.obpt, eloc.carbarn) # This is the same as eobs.ob

# Simulated observations from Carbarn of demoa, they are all below the horizon
demoa.propn.obs_carbarn = tork.aereci(demoa.propn.pvt.position, eloc.carbarn)

############### Satellite visibility from ground observers

demod = Munch() # Example with SGP4 propagation
# tell.mestrt(tell.spacetrack_latest(stclient, 41335)) # So that we can reproduce the results
demod.sentinel3A =\
    tell.makemest({'sma': (np.float64(7180.799), None, 'km'),
                    'ecc': (np.float64(8.89e-05), None, ''),
                    'inc': (np.float64(98.6296), None, 'deg'),
                    'raan': (np.float64(65.5838), None, 'deg'),
                    'argper': (np.float64(97.7597), None, 'deg'),
                    'ma': (np.float64(262.3685), None, 'deg'),
                    'memo': (np.float64(14.26736057), None, 'revolution / d'),
                    'memod': (np.float64(1.43e-06), None, 'revolution / d2'),
                    'memodd': (np.float64(0.0), None, 'revolution / d3'),
                    'period': (np.float64(100.93), None, 'min'),
                    'peralt': (np.float64(802.025), None, 'km'),
                    'apoalt': (np.float64(803.302), None, 'km'),
                    'B': (np.float64(0.0009831673392), None, 'm2 / kg')},
                   '2025-12-26 18:24:37.422',
                   ('1 41335U 16011A   25360.76709979  .00000143  00000-0  77162-4 0  9990',
                    '2 41335  98.6296  65.5838 0000889  97.7597 262.3685 14.26736057513475'),
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

demoa.propn.vis = Munch() # Example with numerical propagation defined in demoa
demoa.propn.vis.events = {'altitude': 125.0*u.km, 'eclipse': [], \
                         'visibility': [tell.observer_location(eloc.mcd, "McDonald")]}
demoa.propn.vis.prep = tork.prepare(demoa.init.pvt, 8*u.hour, demoa.propn.vis.events)
demoa.propn.vis.ephem = tork.propagate(demoa.propn.vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))
demoa.propn.vis.mcdonald = demoa.propn.vis.prep['visibility']['McDonald']
