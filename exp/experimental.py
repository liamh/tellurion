import numpy as np
import astropy.units as u
import posvel
import orkinit
import pvork
import prop
import force
import astro
import apdttm
import element
from astropy.timeseries import TimeSeries

################ Experimental

p0 = [5740.13268349, 3314.06715   ,    0.]
v0 = [-2.75082684,  4.76457184,  5.50165367]
newyear = apdttm.dttm('2025-01-01T00:00:00')
pvt0 = posvel.pvt((p0,v0,newyear))
cgen = prop.generate(pvt0, 86400.0) # Use cgen for any propagation up to 1 day

# The example pvt as a CartesianOrbit

proptimes = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
ceph = prop.propagate(cgen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table
ceph_has_pvt0 = posvel.pvt(ceph, 0) == pvt0
pvt15m = posvel.pvt(ceph, 3) # PVT for 15min by index
pvt35m = posvel.pvt(ceph,'2025-01-01 00:35:00') # PVT for 35min by time
pvt45m = posvel.pvt(ceph, 45*u.min) # PVT for 45min by relative time
pvt1h = posvel.pvt(ceph) # PVT at the end of the ephemeris
pvt12h = prop.propagate(cgen, 12*u.hour) # Propagate to a single time, as a PVT

pvtshift = posvel.pvt(pvt12h,1*u.day)

from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.orbits import Orbit, CartesianOrbit, OrbitType

ckep0 = cgen.initialState.orbit.keplerianorbit()

orb1h = pvork.orkpvt(*pvt1h).cartesianorbit()
kep1h = pvork.orkpvt(*pvt1h).keplerianorbit()
# pvt1h, pvork.pvtork(kep1h) are equal but can't be compared, u.allclose does not work on the pv part
# https://github.com/astropy/astropy/issues/17602


kep0 = element.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "raan":217.4, "ma":7.25}, \
                      apdttm.dttm('2023-09-14T08:30:00'))
kgen = prop.generate(kep0, 86400.0)
keph = prop.propagate(kgen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table


# from experimental import *
