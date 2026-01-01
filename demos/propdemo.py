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

newyear = tell.abstime('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour

################ State and elements

demoa = Munch()
demoa.init = Munch()
demoa.init.check = Munch()

demoa.init.pv = [5740.13268349, 3314.06715, 0., -2.75082684, 4.76457184, 5.50165367]
demoa.init.pvt = tell.pvtcart(demoa.init.pv, newyear)
# Cartesian or Kepler transformation
demoa.init.altper = tork.elementval(demoa.init.pvt, 'altper')  # Altitude of perigee for the initial state
demoa.init.kep = tork.kepler(demoa.init.pvt)  # Convert PVT to Kepler elements
demoa.init.seekep = demoa.init.kep[0].to_dict()  # Easier to read Kepler elements
demoa.init.cart = tork.pvt(demoa.init.kep)  # Convert back to Cartesian, same as demoa.init.pvt

################ Propagation two-body

demoa.propn = Munch()
demoa.propn.check = Munch()

### Numerical
# The generator
demoa.propn.gen = tork.prepare(demoa.init.pvt, 1*u.day) # Use generator for any propagation up to 1 day
# Ephemeris
demoa.propn.ephem = tork.propagate(demoa.propn.gen, prop5m1h, True)  # Propagate to each step, and include the initial state
demoa.propn.pvt = tork.propagate(demoa.propn.gen, prop5m1h, True, output='pvt') # As PVT
demoa.propn.radec = demoa.propn.pvt.spherical # ra, dec, distance
demoa.propn.pvt12h = tork.propagate(demoa.propn.gen, 12*u.hour) # Propagate to a single time, as a PVT
demoa.propn.altperapo = tork.tselements(demoa.propn.ephem, ["altper","altapo"]) # Time series of altitudes of perigee and apogee
demoa.propn.pvapa = tell.hcat(demoa.propn.ephem, demoa.propn.altperapo) # Ephemeris table with additional columns for perigee and apogee altitude

# STM calculation
demoa.props = Munch()
demoa.props.wstm = {'altitude': 125.0*u.km, 'eclipse': [], 'visibility': [], 'stm': True}
demoa.props.gen = tork.prepare(demoa.init.pvt, 1*u.day, demoa.props.wstm) # Use generator for any propagation up to 1 day
demoa.props.finstm = demoa.props.gen['final']['stm'] # STM of final state wrt initial state

# Eclipsing
demoa.propn.eclipse = Munch()
demoa.propn.eclipse.events = {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []}
demoa.propn.eclipse.genevpvt = tork.prepare(demoa.init.pvt, 8*u.hour, demoa.propn.eclipse.events, output='pvt')
demoa.propn.eclipse.suntrans = demoa.propn.eclipse.genevpvt['sun transition'].ephemeris()
demoa.propn.eclipse.ephempvt = tork.propagate(demoa.propn.eclipse.genevpvt, np.linspace(0.25*u.hour, 8*u.hour, 32), output='pvt')
demoa.propn.eclipse.ephem = demoa.propn.eclipse.ephempvt.ephemeris()
demoa.propn.eclipse.merged = \
    demoa.propn.eclipse.genevpvt['sun transition'].merge(demoa.propn.eclipse.ephempvt).ephemeris()
demoa.propn.eclipse.um4h = tork.propagate(demoa.propn.eclipse.genevpvt, 4*u.hour)
demoa.propn.eclipse.um4h15m = tork.propagate(demoa.propn.eclipse.genevpvt, tell.tq('4hr 15min'))

### Analytical
demoa.propa = Munch()
demoa.propa.gen = tork.prepare(demoa.init.pvt, 1*u.day, forceenv=tork.kepleranalytic())
demoa.propa.ephem = tork.propagate(demoa.propa.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
# Difference between analytical and numerical
demoa.propa.andiff = tell.magdiff(demoa.propa.ephem['position'], demoa.propn.ephem['position'])

################ Propagation Kepler element

demob = Munch()

demob.kep = tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           tell.abstime('2023-09-14T08:30:00'))
demob.pvt = tork.pvt(demob.kep) # Convert Kepler elements to PVT
demob.gen = tork.prepare(demob.pvt, 1*u.day)
demob.eph = tork.propagate(demob.gen, prop5m1h, True)  # Propagate to each step, and include the initial state in the ephemeris table
# tell.pvt(demob.eph, 35*u.min)
# fails with KeyError: 'No matches found for key 2023-09-14 09:05:00'
# but tell.pvt(demob.eph,['2023-09-14 09:05:00']) works

demob.kep20m_pvt = demob.eph[4]      # PVT at 20 minutes from ephemeris

################ Time series selection and manipulation

demoa.tssel = Munch()

demoa.tssel.pvt15m = demoa.propn.ephem[3].pvt() # PVT for 15min by index
demoa.tssel.pvt35m = demoa.propn.ephem.loc['2025-01-01 00:35:00']
demoa.tssel.pvt45m = demoa.propn.ephem.loc[tell.abstime(45*u.min, newyear)].pvt() # PVT for 45min by relative time
demoa.tssel.pvt1h = demoa.propn.ephem[-1].pvt() # PVT at the end of the ephemeris
demoa.tssel.kep1h = tork.kepler(demoa.tssel.pvt1h)  # Convert PVT to Kepler elements

################ Propagation with perturbations

demoa.propp = Munch()
demoa.propp.fe4x4 = tork.setgravity(4,4)
demoa.propp.gen4x4 = tork.prepare(demoa.init.pvt, 1*u.day, forceenv=demoa.propp.fe4x4)
demoa.propp.ephem4x4 = tork.propagate(demoa.propp.gen4x4, prop5m1h, True)
demoa.propp.ephem4x4_posdiff = tell.magdiff(demoa.propp.ephem4x4['position'], demoa.propn.ephem['position'])

# It would be nice to have the ability to assemble table with any columns, appropriately renamed, generalize hcat

# Atmospheric drag
demoa.propp.fe4x4hpB01 = tork.dragforce(demoa.propp.fe4x4, 'hp')
demoa.propp.gen4x4hpB01 = tork.prepare(demoa.init.pvt, 1*u.day, forceenv=demoa.propp.fe4x4hpB01)
demoa.propp.ephem4x4hpB01 = tork.propagate(demoa.propp.gen4x4hpB01, prop5m1h, True)
demoa.propp.ephem4x4hpB01_posdiff = tell.magdiff(demoa.propp.ephem4x4hpB01['position'], demoa.propn.ephem['position'])

# Atmospheric drag with Jacobian calculation
demoa.props.fe4x4hpB01 = tork.dragforce(demoa.propp.fe4x4, 'hp', [1.0, True], [1.0, True])
demoa.props.gen4x4hpB01 = tork.prepare(demoa.init.pvt, 1*u.day, demoa.props.wstm, forceenv=demoa.props.fe4x4hpB01)
demoa.props.finstm4x4hpB01 = demoa.props.gen4x4hpB01['final']['stm'] # STM of final state wrt initial state
# Compare demoa.props.finstm with demoa.props.finstm4x4hpB01, they are very different

# Effects of atmospheric drag in NTW and RSW (LVLH) relative coordinates
demoa.propp.ss4x4hpB01 = tork.propagate(demoa.propp.gen4x4hpB01, prop5m1h, True, output='ss')
demoa.propp.ss4x4 = tork.propagate(demoa.propp.gen4x4, prop5m1h, True, output='ss')
demoa.propp.ntw_4x4hpB01_to_4x4 = tork.ntw(demoa.propp.ss4x4hpB01, demoa.propp.ss4x4, tell.siunits)
demoa.propp.rsw_4x4hpB01_to_4x4 = tork.lvlh(demoa.propp.ss4x4hpB01, demoa.propp.ss4x4, tell.siunits)

# Lifetime - need a very low orbit to avoid a long integration, but don't go below 100km altitude, HP will fail
democ = Munch()
democ.kep = tell.kepler({"sma":6600.0, "ecc":0.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25}, newyear)
democ.pvt = tork.pvt(democ.kep) # Convert Kepler elements to PVT
democ.gen = tork.prepare(democ.pvt, 10*u.day, forceenv=demoa.propp.fe4x4hpB01) # Ask for 10 days, but it only lasts about 6
democ.tspan = tork.timerange(democ.gen) # <Quantity 518050.233819 s> ; time until altitude threshold is hit
democ.tspan_dhms = tell.tc(democ.tspan) # '5d 23hr 54min 10.233s'
democ.nhours = np.floor(democ.tspan.to(u.hour)) # Step by an hour for 143 hours, the maximum integer hour
democ.step1hmax = np.linspace(1.0*u.hour, democ.nhours, np.int64(democ.nhours))
democ.eph1hmax = tork.propagate(democ.gen, democ.step1hmax, True) # Ephemeris every hour until it decays
democ.altperapo = tork.tselements(democ.eph1hmax, ["altper","altapo"]) # Altitudes of perigee and apogee every hour
