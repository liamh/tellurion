from elements import *
from propagate import *
import math

# Semimajor axis and eccentricity from apogee and perigee altitudes
def zazptoae(zapo,zper):
    sma = (zapo+zper)/2.0+earthrad
    return([sma, abs((zapo-zper)/(2*sma))])

# Scalar multiplication of velocity vector
def rescalevel(pvt,mult):
    return(cartorb(pvt.position, pvt.velocity.scalarMultiply(mult), pvt.date))

# Find first-orbit perigee and apogee altitudes changing the velocity
# from the two-body velocity that meet the specified perigee and
# apogee
def perapo1orb(zapo, zper, otherels, epoch, velscale):
    ret = {'zapo': zapo, 'zper': zper}
    ret['ae'] = zazptoae(ret['zapo'], ret['zper'])
    args = [epoch] + ret['ae'] + otherels + [meananom]
    ret['kepoes'] = kepler_oes(*args)
    ret['scaled'] = rescalevel(orbpvt(ret['kepoes']), velscale)
    ret['period'] = period(ret['scaled'])
    ret['proptime'] = float(math.ceil(ret['period']+10.0))
    ret['prop1orb'] = prop(ret['scaled'], ret['proptime'])
    ret['alt1orb'] = [posmag(ephlookup(ret['prop1orb'], float(t)))-earthrad
                              for t in range(0,math.ceil(ret['proptime']),10)]
    ret['perapo1orb'] = [min(ret['alt1orb']), max(ret['alt1orb'])]
    return(ret)
