from elements import *
from propagate import *
import math
from scipy import signal
import numpy as np

# Semimajor axis and eccentricity from apogee and perigee altitudes
def zazptoae(zapo,zper):
    sma = (zapo+zper)/2.0+okct['earthrad']
    return([sma, abs((zapo-zper)/(2*sma))])

# Scalar multiplication of velocity vector
def rescale(pvt,posscale,velscale):
    return(cartorb(pvt.position.scalarMultiply(posscale), pvt.velocity.scalarMultiply(velscale), pvt.date))

# Find first-orbit perigee and apogee altitudes changing the velocity
# from the two-body velocity that meet the specified perigee and
# apogee
def peraponorb(zapo, zper, otherels, epoch, numorb=1, scale=[1.0,1.0]):
    ret = {'zapo': zapo, 'zper': zper}
    ret['ae'] = zazptoae(ret['zapo'], ret['zper'])
    args = [epoch] + ret['ae'] + otherels + [okct['meananom']]
    ret['kepoes'] = kepler_oes(*args)
    ret['scaled'] = rescale(orbpvt(ret['kepoes']), *scale)
    ret['period'] = period(ret['scaled'])
    ret['reqproptime'] = float(math.ceil(numorb*ret['period']+10.0))
    ret['prop'] = prop(ret['scaled'], ret['reqproptime'], scB010)
    ret['mindate'] = ret['prop'].minDate
    ret['maxdate'] = ret['prop'].maxDate
    ret['actproptime'] = ret['maxdate'].offsetFrom(ret['mindate'],okct['utc'])
    ret['shortfall'] = ret['reqproptime'] - ret['actproptime']
    ret['alt'] = [posmag(ephlookup(ret['prop'], float(t)))-okct['earthrad']
                  for t in range(0,math.floor(ret['actproptime']),1)]
    ret['perapo'] = [min(ret['alt']), max(ret['alt'])]
    ret['decay'] = decay(ret)
    return(ret)

# Find the apogee and perigee altitudes on each orbit
def decay(data):
    alts = data['alt']
    vinds = list(signal.argrelextrema(np.array(alts), np.less)[0])
    pinds = list(signal.argrelextrema(np.array(alts), np.greater)[0])
    return ([alts[i] for i in pinds], [alts[i] for i in vinds])
