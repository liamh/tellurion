""" Demonstration of state, element, propagation

See
demoa.keys()
demoa.init.keys()
demoa.prop.keys()
"""

import numpy as np
import astropy.units as u
from munch import Munch

import tellurion.core as tell
import tellurion.ork as tork

################ General

newyear = tell.dttm('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour

################ State and elements

demoa = Munch()
demoa.init = Munch()
demoa.init.check = Munch()

demoa.init.pos = [5740.13268349, 3314.06715   ,    0.]
demoa.init.vel = [-2.75082684,  4.76457184,  5.50165367]
demoa.init.pvt = tell.pvt((demoa.init.pos, demoa.init.vel, newyear))  # A tuple (Quantity, Time)
# Cartesian or Kepler transformation
demoa.init.altper = tork.elementval(demoa.init.pvt, 'altper')  # Altitude of perigee for the initial state
demoa.init.kep = tork.kepler(demoa.init.pvt)  # Convert PVT to Kepler elements
demoa.init.seekep = tell.splitsq(demoa.init.kep[0])  # Easier to read Kepler elements
demoa.init.cart = tork.cartesian(demoa.init.kep)  # Convert back to Cartesian, same as demoa.init.pvt

################ Propagation two-body

demoa.prop = Munch()
demoa.prop.check = Munch()

# The generator
demoa.prop.gen = tork.generate(demoa.init.pvt, 1*u.day) # Use generator for any propagation up to 1 day

# The example pvt as a CartesianOrbit
demoa.prop.ephem = tork.propagate(demoa.prop.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
demoa.prop.pvt12h = tork.propagate(demoa.prop.gen, 12*u.hour) # Propagate to a single time, as a PVT
demoa.prop.altperapo = tork.tselements(demoa.prop.ephem, ["altper","altapo"]) # TimeTable of altitudes of perigee and apogee
demoa.prop.pvapa = tell.hcat(demoa.prop.ephem, demoa.prop.altperapo) # Ephemeris table with additional columns for perigee and apogee altitude

# Eclipsing
demoa.prop.eclipse = Munch()
demoa.prop.eclipse.events = {'altitude': 125.0*u.km, 'eclipse': True, 'visibility': []}
demoa.prop.eclipse.genev = tork.generate(demoa.init.pvt, 1*u.day, tork.deffe, demoa.prop.eclipse.events)
demoa.prop.eclipse.suntrans = demoa.prop.eclipse.genev['sun transition']
demoa.prop.eclipse.ephem = tork.propagate(demoa.prop.eclipse.genev, np.linspace(5.0*u.minute, 5*60.0*u.minute, 60))
demoa.prop.eclipse.um12h = tork.propagate(demoa.prop.eclipse.genev, 12*u.hour)
demoa.prop.eclipse.um12h05m = tork.propagate(demoa.prop.eclipse.genev, tell.tq('12hr 5min'))

################ Propagation Kepler element

demob = Munch()

demob.kep = tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           tell.dttm('2023-09-14T08:30:00'))
demob.pvt = tork.cartesian(demob.kep) # Convert Kepler elements to PVT
demob.gen = tork.generate(demob.pvt, 1*u.day)
demob.eph = tork.propagate(demob.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
# tell.pvt(keph, 35*u.min)
# fails with KeyError: 'No matches found for key 2023-09-14 09:05:00'
# but tell.pvt(keph,['2023-09-14 09:05:00']) works

demob.kep20m_pvt = tell.pvt(demob.eph[4])      # PVT at 20 minutes from ephemeris

################ Time series selection and manipulation

demoa.tssel = Munch()

demoa.tssel.pvt15m = tell.pvt(demoa.prop.ephem, 3) # PVT for 15min by index
demoa.tssel.pvt35m = tell.pvt(demoa.prop.ephem,'2025-01-01 00:35:00') # PVT for 35min by time
demoa.tssel.pvt45m = tell.pvt(demoa.prop.ephem, 45*u.min) # PVT for 45min by relative time
demoa.tssel.pvt1h = tell.pvt(demoa.prop.ephem) # PVT at the end of the ephemeris
demoa.tssel.pvtshift = tell.pvt(demoa.prop.pvt12h,1*u.day) # Shift the same posvel to 1 day later
demoa.tssel.kep1h = tork.kepler(demoa.tssel.pvt1h)  # Convert PVT to Kepler elements

################ Propagation with perturbations

demoa.prop.fe4x4 = tork.setgravity(4,4)
demoa.prop.gen4x4 = tork.generate(demoa.init.pvt, 1*u.day, demoa.prop.fe4x4)
demoa.prop.ephem4x4 = tork.propagate(demoa.prop.gen4x4, prop5m1h, True)
demoa.prop.ephem4x4_posdiff = tell.magdiff(demoa.prop.ephem4x4['position'], demoa.prop.ephem['position'])

# It would be nice to have the ability to assemble table with any columns, appropriately renamed, generalize hcat

# Atmospheric drag
demoa.prop.fe4x4hpB01 = tork.dragforce(demoa.prop.fe4x4, 'hp')
demoa.prop.gen4x4hpB01 = tork.generate(demoa.init.pvt, 1*u.day, demoa.prop.fe4x4hpB01)
demoa.prop.ephem4x4hpB01 = tork.propagate(demoa.prop.gen4x4hpB01, prop5m1h, True)
demoa.prop.ephem4x4hpB01_posdiff = tell.magdiff(demoa.prop.ephem4x4hpB01['position'], demoa.prop.ephem['position'])

# Lifetime - need a very low orbit to avoid a long integration, but don't go below 100km altitude, HP will fail
democ = Munch()
democ.kep = tell.kepler({"sma":6600.0, "ecc":0.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25}, newyear)
democ.pvt = tork.cartesian(democ.kep) # Convert Kepler elements to PVT
democ.gen = tork.generate(democ.pvt, 10*u.day, demoa.prop.fe4x4hpB01) # Ask for 10 days, but it only lasts about 6

# PORTING IN PROGRESS - need to restore timerange()
# democ.tspan = democ.gen.timerange() # <Quantity 518050.232645 s> ; time until altitude threshold is hit
# democ.tspan_dhms = tell.tc(democ.tspan) # '5d 23hr 54min 10.233s'
# democ.nhours = np.floor(democ.tspan.to(u.hour)) # Step by an hour for 143 hours, the maximum integer hour
# democ.step1hmax = np.linspace(1.0*u.hour, democ.nhours, np.int64(democ.nhours))
# democ.eph1hmax = tork.propagate(democ.gen, democ.step1hmax, True) # Ephemeris every hour until it decays
# democ.altperapo = tork.tselements(democ.eph1hmax, ["altper","altapo"]) # Altitudes of perigee and apogee every hour
