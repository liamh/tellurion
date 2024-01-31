from elements import *
from propagate import *
import math
from scipy import signal
from scipy.optimize import fsolve
from org.hipparchus.geometry.euclidean.threed import Vector3D

import numpy as np

# Find the initial position and velocity by rescaling the values
# obtained from the two-body solution that matches the desired perigee
# and apogee altitudes
#
# Use peraposolve() to find these factors

# fixargs={'zapo': 325.e3, 'zper': 250.0e3, 'otherels':i45otherels, 'epoch': epoch2022}
# root = peraposolve([1.0, 1.0], fixargs)

# Semimajor axis and eccentricity from apogee and perigee altitudes
def zazptoae(zapo,zper):
    sma = (zapo+zper)/2.0+okct['earthrad']
    return([sma, abs((zapo-zper)/(2*sma))])

# Scalar multiplication of velocity vector
def rescale(pvt,posscale,velscale):
    #print('posscale: ', posscale,'  velscale: ', velscale)
    return(cartorb(pvt.position.scalarMultiply(float(posscale)), pvt.velocity.scalarMultiply(float(velscale)), pvt.date))

# Find perigee and apogee altitudes over the specified number of orbits
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

def decay(data):
    alts = data['alt']
    vinds = list(signal.argrelextrema(np.array(alts), np.less)[0])
    pinds = list(signal.argrelextrema(np.array(alts), np.greater)[0])
    return ([alts[i] for i in pinds], [alts[i] for i in vinds])

def perapodiff(scale, args):
    oneorb = peraponorb(args['zapo'], args['zper'], args['otherels'], args['epoch'], 1, scale)
    return oneorb['perapo']-np.array([args['zper'], args['zapo']])

def peraposolve(scale, args):
    root = fsolve(perapodiff, scale, args).tolist()
    diffs = perapodiff(root, args)
    #print ('first perigee error (m): ', diffs[0], ' first apogee error (m): ', diffs[1])
    return (root)

maxrevs = 1000

# Lifetime using full Cartesian integration
def ltfci(zper, zapo, inc, startyr, startmo, gravdeg, gravord):
    grav = GravityFieldFactory.getNormalizedProvider(gravdeg, gravord)
    okct['gravity']=grav
    otherels = [inc, 30.0, 0.0, 0.0]
    epoch = datm(startyr, startmo, 1)
    fixargs={'zapo': zapo, 'zper': zper, 'otherels':otherels, 'epoch': epoch}
    ret = {'gravdeg': gravdeg,
           'gravord': gravord,
           'epoch': epoch,
           'pvscale': peraposolve([1.0, 1.0], fixargs)
           }
    ret['proporb'] = peraponorb(zapo, zper, otherels, epoch, maxrevs, ret['pvscale'])
    if (ret['proporb']['shortfall']) > 1.0e5:
        ret['lifetime_days'] = ret['proporb']['actproptime']/day
    else:
        ret['lifetime_days'] = np.inf
    return(ret)

# Difference in lifetime between 0x0+drag and 10x10+drag force models
# Returned ltdiff: days difference, ltgrav: lifetime days with gravity, pctdif: percent difference
#  ltgravdiff(250.0e3, 325.0e3, 45.0, 2019, 12)
def ltgravdiff(zper, zapo, inc, startyr, startmo):
    ltgrav = ltfci(zper, zapo, inc, startyr, startmo, 10, 10)
    ltkep = ltfci(zper, zapo, inc, startyr, startmo, 0, 0)
    ltdiff = ltkep['lifetime_days']-ltgrav['lifetime_days']
    return({'zapo': zapo, 'zper': zper, 'inc': inc, 'year': startyr, 'month': startmo, # inputs
            'ltdiff': ltdiff, 'ltgrav': ltgrav['lifetime_days'],
            'pctdiff': 100.0*ltdiff/ltgrav['lifetime_days']}) # outputs
