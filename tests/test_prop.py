from munch import Munch
import numpy as np
import astropy.units as u
import tellurion as tell

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

demoa.init.pv = np.array([5740.13268349, 3314.06715, 0., -2.75082684, 4.76457184, 5.50165367])
demoa.init.pvt = tell.pvtcart(demoa.init.pv, newyear)
# Cartesian or Kepler transformation
demoa.init.kep = demoa.init.pvt.kepler()  # Convert PVT to Kepler elements
demoa.init.cart = tell.pvt(demoa.init.kep)  # Convert back to Cartesian, same as demoa.init.pvt

demob = Munch()
demob.kep = tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           tell.abstime('2023-09-14T08:30:00'))
demob.pvt = tell.pvt(demob.kep) # Convert Kepler elements to PVT
demob.rekep = demob.pvt.kepler()

#############################
####  SGP4 mean elements ####
#############################

sent3a = Munch()
# Created from
# isssent = tell.spacetrack_latest(stclient, [25544, 41335])
# tell.mestrt(isssent['SENTINEL 3A'])
# sent3a.mest = \
#     tell.makemest({'sma': (np.float64(7180.799), None, 'km'),
#                     'ecc': (np.float64(8.89e-05), None, ''),
#                     'inc': (np.float64(98.6296), None, 'deg'),
#                     'raan': (np.float64(65.5838), None, 'deg'),
#                     'argper': (np.float64(97.7597), None, 'deg'),
#                     'ma': (np.float64(262.3685), None, 'deg'),
#                     'memo': (np.float64(14.26736057), None, 'revolution / d'),
#                     'memod': (np.float64(1.43e-06), None, 'revolution / d2'),
#                     'memodd': (np.float64(0.0), None, 'revolution / d3'),
#                     'period': (np.float64(100.93), None, 'min'),
#                     'peralt': (np.float64(802.025), None, 'km'),
#                     'apoalt': (np.float64(803.302), None, 'km'),
#                     'B': (np.float64(0.0009831673392), None, 'm2 / kg')},
#                    '2025-12-26 18:24:37.422',
#                    ('1 41335U 16011A   25360.76709979  .00000143  00000-0  77162-4 0  9990',
#                     '2 41335  98.6296  65.5838 0000889  97.7597 262.3685 14.26736057513475'),
#                    'SGP4',
#                    {'name': 'SENTINEL 3A',
#                     'type': 'PAYLOAD',
#                     'catid': 41335,
#                     'intldes': '2016-011A'})
# sent3a.gen = tell.prepare(sent3a.mest, 1*u.day, {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []})

##########################
####      Tests       ####
##########################

def test_kepcart():
    np.testing.assert_allclose(demoa.init.pvt.to_array(), \
                               np.array([ 5.74013268e+06,  3.31406715e+06,  0.00000000e+00, \

                                          -2.75082684e+03, 4.76457184e+03,  5.50165367e+03, \
                                          6.067600e+04]))

def test_cartkep():
    return kepequal(demob.rekep, demob.kep)

# Temporarily removed this test until a general serialization capability is added
# def test_meanels():
#     np.testing.assert_allclose(sent3a.gen['pvt0'].time.to_value('jd'), 2461036.26709979)
#     np.testing.assert_allclose(sent3a.gen['pvt0'].pv.si.value[0], \
#                                np.array([3007467.96362189, 6523970.98339558,   -7829.3325051 ]))
#     np.testing.assert_allclose(sent3a.gen['pvt0'].pv.si.value[1], \
#                                np.array([1029.78212853, -475.05478146, 7363.69699183]))
