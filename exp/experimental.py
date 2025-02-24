""" Demonstration of state, element, propagation

See
demoa.keys()
demoa.init.keys()
demoa.prop.keys()
"""

import numpy as np
import astropy.units as u
from munch import Munch
import posvel
import astro
import element
import ork.init
import ork.posvel
import ork.element
import ork.force
import ork.prop


newyear = posvel.dttm('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour

################ State and elements

demoa = Munch()
demoa.init = Munch()
demoa.init.check = Munch()

demoa.init.pos = [5740.13268349, 3314.06715   ,    0.]
demoa.init.vel = [-2.75082684,  4.76457184,  5.50165367]
demoa.init.dttm = newyear
demoa.init.pvt = posvel.pvt((demoa.init.pos,demoa.init.vel,demoa.init.dttm))  # A tuple (Quantity, Time)
demoa.init.kep = ork.element.kepler(demoa.init.pvt)  # Convert PVT to Kepler elements
demoa.init.seekep = astro.splitsq(demoa.init.kep[0])  # Easier to read Kepler elements

# Orekit representations
demoa.init.ork = Munch()
demoa.init.ork.pvt = ork.posvel.tspvc(*demoa.init.pvt)           # org.orekit.utils.TimeStampedPVCoordinates
demoa.init.ork.cartorb = demoa.init.ork.pvt.cartesianorbit()  # org.orekit.orbits.CartesianOrbit
demoa.init.ork.keporb = demoa.init.ork.pvt.keplerianorbit()   # org.orekit.orbits.KeplerianOrbit

# AstroPy representations
demoa.init.check.kepap = demoa.init.ork.pvt.kepler()

# Get the pvt back
demoa.init.check.opvt0 = demoa.init.ork.pvt.pvt()
demoa.init.check.ocartorb0 = demoa.init.ork.cartorb.pvt()
demoa.init.check.okeporb0 = demoa.init.ork.keporb.pvt()

################ Propagation two-body

demoa.prop = Munch()
demoa.prop.check = Munch()

# The generator
demoa.prop.gen = ork.prop.generate(demoa.init.pvt, 1*u.day) # Use generator for any propagation up to 1 day

# The example pvt as a CartesianOrbit

demoa.prop.ephem = ork.prop.propagate(demoa.prop.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
demoa.prop.pvt12h = ork.prop.propagate(demoa.prop.gen, 12*u.hour) # Propagate to a single time, as a PVT

demoa.prop.ckep0 = demoa.prop.gen.kepler() # Convert pvt0 initial state directly from generator
demoa.prop.cpvt0 = demoa.prop.gen.pvt() # Convert pvt0 initial state directly from generator

demoa.prop.altperapo = ork.element.tselements(demoa.prop.ephem, ["altper","altapo"]) # TimeTable of altitudes of perigee and apogee
demoa.prop.pvapa = astro.hcat(demoa.prop.ephem, demoa.prop.altperapo) # Ephemeris table with additional columns for perige and apogee altitude

demoa.prop.check.eph_has_pvt0 = posvel.pvt(demoa.prop.ephem, 0) == demoa.init.pvt # Check that initial state is in the ephemeris

# orb1h_ork = ork.posvel.tspvc(*pvt1h).cartesianorbit()
# pvt1h, ork.posvel._pvtork(kep1h) are equal but can't be compared, u.allclose does not work on the pv part
# https://github.com/astropy/astropy/issues/17602, it is a numpy issue https://github.com/numpy/numpy/issues/28104

# Eclipsing
demoa.prop.eclipse = Munch()
demoa.prop.eclipse.events = {'altitude': 125.0*u.km, 'eclipse': ['umbra'], 'visibility': []}
demoa.prop.eclipse.genev = ork.prop.generate(demoa.init.pvt, 1*u.day, ork.force.deffe, demoa.prop.eclipse.events)
demoa.prop.eclipse.suntrans = demoa.prop.eclipse.genev['sun transition']
demoa.prop.eclipse.um12h = ork.prop.propagate(demoa.prop.eclipse.genev, 12*u.hour)
demoa.prop.eclipse.um12h05m = ork.prop.propagate(demoa.prop.eclipse.genev, astro.tq('12hr 5min'))
demoa.prop.eclipse.pvu12h05m = astro.splitsq(demoa.prop.eclipse.um12h05m[0])
prop5m3h = np.linspace(5.0*u.minute, 180.0*u.minute, 36) # Step every 5 minutes for an hour
p5m3h = ork.prop.propagate2(demoa.prop.eclipse.genev, prop5m3h, True)  # Propagate to each step, and include the initial state in the
# <Quantity [[-0.01778085, -0.01778085],   datetime.datetime(2025, 1, 1, 0, 0)
#            [-0.34137358, -0.34137358],   datetime.datetime(2025, 1, 1, 0, 5)
#            [-0.6424991 , -0.6424991 ],   datetime.datetime(2025, 1, 1, 0, 10)
#            [-0.87488452, -0.87488452],   datetime.datetime(2025, 1, 1, 0, 15)
#            [-0.88653962, -0.88653962],   datetime.datetime(2025, 1, 1, 0, 20)
#            [-0.66015284, -0.66015284],   datetime.datetime(2025, 1, 1, 0, 25)
#            [-0.35895314, -0.35895314],   datetime.datetime(2025, 1, 1, 0, 30)
#            [-0.03922183, -0.03922183],   datetime.datetime(2025, 1, 1, 0, 35)
#            [ 0.28506905,  0.28506905],   datetime.datetime(2025, 1, 1, 0, 40) FULL SUN
#            [ 0.60834213,  0.60834213],   datetime.datetime(2025, 1, 1, 0, 45)
#            [ 0.92541815,  0.92541815],   datetime.datetime(2025, 1, 1, 0, 50)
#            [ 1.2253585 ,  1.2253585 ],   datetime.datetime(2025, 1, 1, 0, 55)
#            [ 1.46916708,  1.46916708],   datetime.datetime(2025, 1, 1, 1, 0)
#            [ 1.51652395,  1.51652395],   datetime.datetime(2025, 1, 1, 1, 5)
#            [ 1.30774431,  1.30774431],   datetime.datetime(2025, 1, 1, 1, 10)
#            [ 1.00577482,  1.00577482],   datetime.datetime(2025, 1, 1, 1, 15)
#            [ 0.67814934,  0.67814934],   datetime.datetime(2025, 1, 1, 1, 20)
#            [ 0.34329324,  0.34329324],   datetime.datetime(2025, 1, 1, 1, 25)
#            [ 0.00987637,  0.00987637],   datetime.datetime(2025, 1, 1, 1, 30)
#            [-0.31478495, -0.31478495],   datetime.datetime(2025, 1, 1, 1, 35) UMBRA
#            [-0.61863019, -0.61863019],   datetime.datetime(2025, 1, 1, 1, 40)
#            [-0.86068561, -0.86068561],   datetime.datetime(2025, 1, 1, 1, 45)
#            [-0.89694722, -0.89694722],   datetime.datetime(2025, 1, 1, 1, 50)
#            [-0.68321815, -0.68321815],   datetime.datetime(2025, 1, 1, 1, 55)
#            [-0.38509842, -0.38509842],   datetime.datetime(2025, 1, 1, 2, 0)
#            [-0.06622372, -0.06622372],   datetime.datetime(2025, 1, 1, 2, 5)
#            [ 0.25789154,  0.25789154],   datetime.datetime(2025, 1, 1, 2, 10) FULL SUN
#            [ 0.58137472,  0.58137472],   datetime.datetime(2025, 1, 1, 2, 15)
#            [ 0.89919435,  0.89919435],   datetime.datetime(2025, 1, 1, 2, 20)
#            [ 1.2012902 ,  1.2012902 ],   datetime.datetime(2025, 1, 1, 2, 25)
#            [ 1.45307034,  1.45307034],   datetime.datetime(2025, 1, 1, 2, 30)
#            [ 1.52401611,  1.52401611],   datetime.datetime(2025, 1, 1, 2, 35)
#            [ 1.33028103,  1.33028103],   datetime.datetime(2025, 1, 1, 2, 40)
#            [ 1.03230662,  1.03230662],   datetime.datetime(2025, 1, 1, 2, 45)
#            [ 0.70592575,  0.70592575],   datetime.datetime(2025, 1, 1, 2, 50)
#            [ 0.37133659,  0.37133659],   datetime.datetime(2025, 1, 1, 2, 55)
#            [ 0.03758472,  0.03758472]]>  datetime.datetime(2025, 1, 1, 3, 0)]>









































################ Propagation Kepler element

demob = Munch()

demob.kep = element.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           posvel.dttm('2023-09-14T08:30:00'))
demob.pvt = ork.posvel.pvt(*demob.kep) # Convert Kepler elements to PVT

demob.gen = ork.prop.generate(demob.kep, 1*u.day)
demob.eph = ork.prop.propagate(demob.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
# posvel.pvt(keph, 35*u.min)
# fails with KeyError: 'No matches found for key 2023-09-14 09:05:00'
# but posvel.pvt(keph,['2023-09-14 09:05:00']) works

demob.kep0 = demob.gen.kepler() # Initial state as Kepler
demob.pvt0 = demob.gen.pvt()    # Initial state as PVT

demob.kep20m_pvt = posvel.pvt(demob.eph[4])      # PVT at 20 minutes from ephemeris
demob.kep20m_opvt = ork.posvel.tspvc(demob.eph[4]) # PVT at 20 minutes from ephemeris via Orekit
demob.kep20m_kep = demob.kep20m_opvt.kepler()     # Kepler elements at 20 minutes

################ Time series selection and manipulation

demoa.tssel = Munch()

demoa.tssel.pvt15m = posvel.pvt(demoa.prop.ephem, 3) # PVT for 15min by index
demoa.tssel.pvt35m = posvel.pvt(demoa.prop.ephem,'2025-01-01 00:35:00') # PVT for 35min by time
demoa.tssel.pvt45m = posvel.pvt(demoa.prop.ephem, 45*u.min) # PVT for 45min by relative time
demoa.tssel.pvt1h = posvel.pvt(demoa.prop.ephem) # PVT at the end of the ephemeris
#demoa.tssel.pvtshift = posvel.pvt(demoa.prop.pvt12h,1*u.day) # Shift the same posvel to 1 day later
demoa.tssel.kep1h = ork.element.kepler(demoa.tssel.pvt1h)  # Convert PVT to Kepler elements

################ Propagation with perturbations

demoa.prop.fe4x4 = ork.force.setgravity(4,4)
demoa.prop.gen4x4 = ork.prop.generate(demoa.init.pvt, 1*u.day, demoa.prop.fe4x4)
demoa.prop.ephem4x4 = ork.prop.propagate(demoa.prop.gen4x4, prop5m1h, True)
demoa.prop.ephem4x4_posdiff = posvel.magdiff(demoa.prop.ephem4x4['position'], demoa.prop.ephem['position'])

# It would be nice to have the ability to assemble table with any columns, appropriately renamed, generalize hcat

# Atmospheric drag
demoa.prop.fe4x4hpB01 = ork.force.dragforce(demoa.prop.fe4x4, 'hp')
demoa.prop.gen4x4hpB01 = ork.prop.generate(demoa.init.pvt, 1*u.day, demoa.prop.fe4x4hpB01)
demoa.prop.ephem4x4hpB01 = ork.prop.propagate(demoa.prop.gen4x4hpB01, prop5m1h, True)
demoa.prop.ephem4x4hpB01_posdiff = posvel.magdiff(demoa.prop.ephem4x4hpB01['position'], demoa.prop.ephem['position'])

# Lifetime - need a very low orbit to avoid a long integration, but don't go below 100km altitude, HP will fail
democ = Munch()
democ.kep = element.kepler({"sma":6600.0, "ecc":0.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25}, newyear)
democ.pvt = ork.posvel.pvt(*democ.kep) # Convert Kepler elements to PVT
democ.gen = ork.prop.generate(democ.pvt, 10*u.day, demoa.prop.fe4x4hpB01) # Ask for 10 days, but it only lasts about 6
democ.tspan = democ.gen.timerange() # <Quantity 518050.232645 s> ; time until altitude threshold is hit
democ.tspan_dhms = astro.tc(democ.tspan) # '5d 23hr 54min 10.233s'
democ.nhours = np.floor(democ.tspan.to(u.hour)) # Step by an hour for 143 hours, the maximum integer hour
democ.step1hmax = np.linspace(1.0*u.hour, democ.nhours, np.int64(democ.nhours))
democ.eph1hmax = ork.prop.propagate(democ.gen, democ.step1hmax, True) # Ephemeris every hour until it decays
democ.altperapo = ork.element.tselements(democ.eph1hmax, ["altper","altapo"]) # Altitudes of perigee and apogee every hour

############### Earth locations

import astropy.coordinates as coord
import ork.geog as oge
import geog
import geonames

eloc = Munch()

eloc.mcd = coord.EarthLocation.of_site('McDonald Observatory')
eloc.mcdsv = geog.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcd = astro.hcat(demoa.prop.ephem, eloc.mcdsv) # Ephemeris table with additional column for McDonald site vector
eloc.mcdsvork = oge.eciobs(eloc.mcd, demoa.prop.ephem.time, 'mcdonald sitevec')
eloc.cmcdork = astro.hcat(demoa.prop.ephem, eloc.mcdsvork) # Ephemeris table with additional column for McDonald site vector

# Difference between AstroPy and Orekit
eloc.mcdsv_apy = posvel.makepos(geog.eciobs(eloc.mcd, newyear))
eloc.mcdsv_ork = posvel.makepos(oge.eciobs(eloc.mcd, newyear))
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

eobs.ptobork = oge.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ptobapy = geog.eciobs(eloc.carbarn, eobs.ob, 'carbarn obs')
eobs.ob_ork_apy_dist = np.linalg.norm(eobs.ptobork.cartesian.xyz - eobs.ptobapy.cartesian.xyz).si
