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
    print('posscale: ', posscale,'  velscale: ', velscale)
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
    print ('first perigee error (m): ', diffs[0], ' first apogee error (m): ', diffs[1])
    return (root)
