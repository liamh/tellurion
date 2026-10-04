# Geographic locations and observations

import pathlib
import astropy.coordinates as coord
from demos.propdemo import *
from importlib.resources import files
import astropy_hdf5io  # registers all serializers including MeanElementSetT

_MODULE_DIR = pathlib.Path(__file__).parent # To find pre-recorded file sentinel3a.h5 in this directory

############################################################
## Observation demos
############################################################

############### Earth locations

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = tell.sitevec(eloc.mcd, demoa.propn.ephem.time, "McDonald")
eloc.mcdsvts = eloc.mcdsv.ephemeris(False, None)   # Time series of site vectors for McDonald
eloc.cmcd = tell.hcat(demoa.propn.ephem, eloc.mcdsvts) # Ephemeris table with additional columns for McDonald site vector

eloc.kickapoo = tell.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')
# From tell.location('carbarn') - integrate earthloc with location()?
eloc.carbarn = tell.earthloc(lat=38.87206*u.deg, lon=-77.01748*u.deg, elevation=13.0*u.m)
eloc.carbarn_lst_newyear = tell.siderealtime(newyear, eloc.carbarn)

############### Earth observations

eobs = Munch()

# A simulated observation from carbarn
eobs.ob = tell.azelrange(255*u.deg, 80*u.deg, 1455*u.km, eloc.carbarn, newyear)
eobs.obpt = tell.eciaer(eobs.ob)
eobs.obret = tell.aereci(eobs.obpt, eloc.carbarn) # This is the same as eobs.ob

# Simulated observations from Carbarn of demoa, they are all below the horizon
demoa.propn.obs_carbarn = tell.aereci(demoa.propn.pvt.position, eloc.carbarn)

############### Satellite visibility from ground observers

demod = Munch() # Example with SGP4 propagation
demod.sentinel3A = load(_MODULE_DIR / 'sentinel3a.h5')
demod.mcd_cb_vis = Munch() # Visibility from McDonald and Carbarn
demod.mcd_cb_vis.ev = {'altitude': 125.0*u.km, 'eclipse': [], \
                       'visibility': [tell.observer_location(eloc.mcd, "McDonald"), \
                                      tell.observer_location(eloc.carbarn, "Carbarn", 30.0)]}
demod.mcd_cb_vis.prep = tell.prepare(demod.sentinel3A, 1*u.day, demod.mcd_cb_vis.ev)
demod.mcd_cb_vis.mcdonald = demod.mcd_cb_vis.prep['visibility']['McDonald'] # McDonald visibility transitions
demod.mcd_cb_vis.carbarn = demod.mcd_cb_vis.prep['visibility']['Carbarn'] # Carbarn visibility transitions
demod.mcd_cb_vis.ephempvt = tell.propagate(demod.mcd_cb_vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32), output='pvt')
demod.mcd_cb_vis.ephem = tell.propagate(demod.mcd_cb_vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))
demod.mcd_ecl = Munch()  # Visibility from McDonald and eclipses
demod.mcd_ecl.ev = {'altitude': 125.0*u.km, 'eclipse': [True, True], \
                    'visibility': [tell.observer_location(eloc.mcd, "McDonald")]}
demod.mcd_ecl.prep = tell.prepare(demod.sentinel3A, 1*u.day, demod.mcd_ecl.ev)
demod.mcd_ecl.mcdonald = demod.mcd_ecl.prep['visibility']['McDonald'] # McDonald visibility transitions
demod.mcd_ecl.suntrans = demod.mcd_ecl.prep['sun transition']
demod.mcd_ecl.ephem = tell.propagate(demod.mcd_ecl.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))

############ Upward transitions not present

demoa.propn.vis = Munch() # Example with numerical propagation defined in demoa
demoa.propn.vis.events = {'altitude': 125.0*u.km, 'eclipse': [], \
                         'visibility': [tell.observer_location(eloc.mcd, "McDonald")]}
demoa.propn.vis.prep = tell.prepare(demoa.init.pvt, 8*u.hour, demoa.propn.vis.events)
demoa.propn.vis.ephem = tell.propagate(demoa.propn.vis.prep, np.linspace(0.25*u.hour, 8*u.hour, 32))
demoa.propn.vis.mcdonald = demoa.propn.vis.prep['visibility']['McDonald']
