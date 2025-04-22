import numpy as np
import astropy.units as u
from munch import Munch

import tellurion.core as tell
import tellurion.ork as tork

################################
####  Comparison of states  ####
################################

def kepequal(a, b):
    for nm in tell.kepeltma_names:
        np.testing.assert_allclose(a[0][nm],b[0][nm])
    return np.testing.assert_equal(a[1], b[1])

##########################
####   Kepler         ####
##########################

newyear = tell.abstime('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour

demoa = Munch()
demoa.init = Munch()
demoa.init.check = Munch()

demoa.init.pos = [5740.13268349, 3314.06715   ,    0.]
demoa.init.vel = [-2.75082684,  4.76457184,  5.50165367]
demoa.init.pvt = tell.pvt((demoa.init.pos, demoa.init.vel, newyear))  # A tuple (Quantity, Time)
# Cartesian or Kepler transformation
demoa.init.kep = tork.kepler(demoa.init.pvt)  # Convert PVT to Kepler elements
demoa.init.seekep = tell.splitsq(demoa.init.kep[0])  # Easier to read Kepler elements
demoa.init.cart = tork.cartesian(demoa.init.kep)  # Convert back to Cartesian, same as demoa.init.pvt

demob = Munch()
demob.kep = tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           tell.abstime('2023-09-14T08:30:00'))
demob.pvt = tork.cartesian(demob.kep) # Convert Kepler elements to PVT
demob.rekep = tork.kepler(demob.pvt)

#############################
####  SGP4 mean elements ####
#############################

sent3a = Munch()
# Created from
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# tell.mestrt(isssent['SENTINEL 3A'])
sent3a.mest = \
    tell.makemest(((7180.805, 0.0001011, 98.6294, \
                    170.7192, 89.9179, 270.2117, \
                    14.26734118, 3.28e-06, 0.0, 100.93, 801.944, 803.396, \
                    0.0019568549280000003), \
                   ('sma', 'ecc', 'inc', 'raan', 'argper', 'ma', 'memo', \
                    'memod', 'memodd', 'period', 'peralt', 'apoalt', 'B'), \
                   ('km', '', 'deg', 'deg', 'deg', 'deg', 'revolution / d', \
                    'revolution / d2', 'revolution / d3', 'min', 'km', \
                    'km', 'm2 / kg')), \
                  '2025-04-12 04:45:30.072', \
                  ('1 41335U 16011A   25102.19826472  .00000328  00000-0  15358-3 0  9992', \
                   '2 41335  98.6294 170.7192 0001011  89.9179 270.2117 14.26734118476609'), \
                  'SGP4', \
                  {'name': 'SENTINEL 3A', 'type': 'PAYLOAD', 'catid': 41335, \
                   'intldes': '2016-011A'})
sent3a.gen = tork.prepare(sent3a.mest, 1*u.day, {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []})

##########################
####      Tests       ####
##########################

def test_kepcart():
    np.testing.assert_allclose(tell.pvtsijd(demoa.init.pvt), \
                               np.array([ 5.74013268e+06,  3.31406715e+06,  0.00000000e+00, \
                                          -2.75082684e+03, 4.76457184e+03,  5.50165367e+03, \
                                          2.46067650e+06]))

def test_cartkep():
    return kepequal(demob.rekep, demob.kep)

def test_meanels():
    np.testing.assert_allclose(sent3a.gen['pvt0'].t.to_value('jd'), 2460777.69826472)
    np.testing.assert_allclose(sent3a.gen['pvt0'].pv.si.value[0], \
                               np.array([-7083008.74600191,  1198601.43504089, 17351.44198235]))
    np.testing.assert_allclose(sent3a.gen['pvt0'].pv.si.value[1], \
                               np.array([ 212.8109757 , 1100.43163449, 7365.81460009]))
