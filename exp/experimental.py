import numpy as np
import astropy.units as u
import posvel
import orkinit
import pvork
import prop
import force
import astro
import apdttm
from astropy.timeseries import TimeSeries

################ Experimental

p0 = [5740.13268349, 3314.06715   ,    0.]
v0 = [-2.75082684,  4.76457184,  5.50165367]
newyear = apdttm.dttm('2025-01-01T00:00:00')
pvt0 = posvel.pvt((p0,v0,newyear))
gen = prop.generate(pvt0, 86400.0) # Use gen for any propagation up to 1 day

# The example pvt as a CartesianOrbit

proptimes = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
ts = prop.propagate(gen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table
ts_has_pvt0 = posvel.pvt(ts, 0) == pvt0
pvt15m = posvel.pvt(ts, 3) # PVT for 15min by index
pvt35m = posvel.pvt(ts,'2025-01-01 00:35:00') # PVT for 35min by time
pvt45m = posvel.pvt(ts, 45*u.min) # PVT for 45min by relative time
pvt1h = posvel.pvt(ts) # PVT at the end of the ephemeris
pvt12h = prop.propagate(gen, 12*u.hour) # Propagate to a single time, as a PVT

pvtshift = posvel.pvt(pvt12h,1*u.day)


# from experimental import *
