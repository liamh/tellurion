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

myp = [5740.13268349, 3314.06715   ,    0.]
myv = [-2.75082684,  4.76457184,  5.50165367]
mypv = posvel.makepv(myp,myv)
newyear = dttm.todttm('2025-01-01T00:00:00')

# The example pvt as a CartesianOrbit

proptimes = [5*(1+tm) for tm in range(10)]*u.min
ts = prop.ephem((mypv, newyear), proptimes)
(pv3, t3) = posvel.makepvt(ts[3])

# from experimental import *
