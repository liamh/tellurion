import numpy as np
import astropy.units as u
import posvel
import orkinit
import pvork
import prop
import force
import astro
import dttm
from astropy.timeseries import TimeSeries

################ Experimental

p0 = [5740.13268349, 3314.06715   ,    0.]
v0 = [-2.75082684,  4.76457184,  5.50165367]
newyear = dttm.to_dttm('2025-01-01T00:00:00')
pvt0 = posvel.makepvt((p0,v0,newyear))
gen = prop.generate(pvt0, 86400.0) # Use gen for any propagation up to 1 day

# The example pvt as a CartesianOrbit

proptimes = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
ts = prop.propagate(gen, proptimes, True)  # Propagate to each step, and include the initial state in the ephemeris table
ts_has_pvt0 = posvel.makepvt(ts[0]) == pvt0
pvt15m = posvel.makepvt(ts[3]) # PVT for 15min by index
pvt35m = posvel.makepvt(ts.loc[dttm.to_dttm('2025-01-01 00:35:00')]) # PVT for 35min by time
pvt12h = prop.propagate(gen, 12*u.hour) # Propagate to a single time, as a PVT

# from experimental import *
