# Create element sets and convert to both Cartesian and spherical
# (right ascension, declination, distance) position-velocity specification

# Run with pytest -q test_element.py

import numpy as np
import astropy.units as u
import tellurion.core as tell
import tellurion.ork as tork

def keppvt(dictels, time):
    """From the dictionary of element values, compute the Kepler
    elements (sma, ecc, inc, argper, raan, and ma or ta), the PVT
    transformation of those elements, and the spherical coordinates."""
    d = {'els': dictels}
    d['kep'] = tell.kepler(tork.allplane(d['els']), time)
    d['pvt'] = tork.pvt(d['kep'])
    d['pvt'].spherical
    return d

# A LEO circular orbit like the ISS
leo1 = keppvt({"altper": 350*u.km, "altapo": 350*u.km, \
                                  "inc":55.0*u.deg, "argper": 120.0*u.deg, \
                                  "raan": 20.0*u.deg, "ma": 30.0*u.deg},
              tell.abstime('2026-01-01 05:55:00'))

# A LEO near-circular orbit like the HST
leo2 = keppvt({"altper": 525.0*u.km, "altapo": 555.0*u.km, \
               "inc":28.5*u.deg, "argper": 40.0*u.deg, \
               "raan": 40.0*u.deg, "ma": -30.0*u.deg},
              tell.abstime('2026-01-01 09:30:00'))

# A GEO orbit
geo1 = keppvt({"memo":1.0*u.rev/u.sday, \
               "ecc":0.0, \
               "inc":0.0*u.deg, "argper": 120.0*u.deg, \
               "raan": 0.0*u.deg, "ma": 0.0*u.deg}, \
            tell.abstime('2026-01-01 20:30:00'))

# A GEO transfer orbit
ell1 = keppvt({"altper": 350*u.km, "altapo": tork.sma(1.0,True), \
                  "inc":0.0*u.deg, "argper": 120.0*u.deg, \
                  "raan": 0.0*u.deg, "ma": 90.0*u.deg},
                 tell.abstime('2026-01-01 05:55:00'))

# An elliptical orbit
ell2 = keppvt({"altper":160*u.km, "altapo":20250*u.km, \
               "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg}, \
              tell.abstime('2026-01-01 08:05:00'))

# Like Vanguard 1
vang1 = keppvt({"altper": 600.0*u.km, "altapo": 12000.0*u.km,
                "inc":36.0*u.deg, "argper": 140.0*u.deg, \
                "raan": 0.0*u.deg, "ma": 100.0*u.deg},
               tell.abstime('2026-01-01 14:45:00'))

# A semisynchronous orbit like GPS
gps1 = keppvt({"sma": tork.sma(2.0), "ecc": 0.0*u.dimensionless_unscaled,
                "inc":55.0*u.deg, "argper": 0.0*u.deg, \
                "raan": 120.0*u.deg, "ma": 77.0*u.deg},
               tell.abstime('2026-01-01 12:20:00'))

def test_allplane():
    smaecc = tork.allplane({"altper":160*u.km, "altapo":20250*u.km, \
              "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    alts = tork.allplane({"sma":8000*u.km, "ecc":0.1*u.dimensionless_unscaled, \
                          "inc":45*u.deg, "argper": 120.0*u.deg, "raan": 80.0*u.deg, \
                          "ma": 0.0*u.deg})


    np.testing.assert_allclose(smaecc['sma'], 16583.13646*u.km)
    np.testing.assert_allclose(smaecc['ecc'], 0.60573583*u.dimensionless_unscaled)
    np.testing.assert_allclose(alts['altper'], 821.86354*u.km)
    np.testing.assert_allclose(alts['altapo'], 2421.86354*u.km)
    np.testing.assert_allclose(geo1['kep'].els['sma'], 42164.1696233*u.km)
    # Test conversion to semimajor axis
    np.testing.assert_allclose(tork.sma(1e4*u.s), 10032.11910363*u.km)
    np.testing.assert_allclose(tork.sma(-20*(u.km/u.s)**2), 9965.0110375*u.km)
    np.testing.assert_allclose(ell1['kep'].els['sma'], 24446.15304165*u.km)

def test_posvel():
    np.testing.assert_allclose(leo1['pvt'].cartesian['position'].si.value, \
                               np.array([-6135286.90982537,  -179677.30883764,  2755683.36773215]))
    np.testing.assert_allclose(leo1['pvt'].cartesian['velocity'].si.value, \
                               np.array([-2308.74628271, -4909.03309813, -5460.30174535]))
    np.testing.assert_allclose(np.array(leo1['pvt'].spherical.si.value.tolist()), \
                               np.array([ 3.17087017e+00, 4.21989266e-01, 6.72813646e+06, \
                                          7.88434305e-04, -8.89601706e-04, -2.66453526e-12]))
    np.testing.assert_allclose(leo2['pvt'].cartesian['position'].si.value, \
                               np.array([4542282.92387955, 5170057.21334251,  565092.22923442]))
    np.testing.assert_allclose(leo2['pvt'].cartesian['velocity'].si.value, \
                               np.array([-5236.82832614,  4199.24378643,  3574.26419318]))
    np.testing.assert_allclose(np.array(leo2['pvt'].spherical.si.value.tolist()), \
                               np.array([8.49945142e-01, 8.19279154e-02, 6.90515423e+06,\
                                         9.74389282e-04, 5.19462927e-04, -8.25996492e+00]))
    np.testing.assert_allclose(geo1['pvt'].cartesian['position'].si.value, \
                               np.array([-21082084.8116475, 36515242.02324963, 0.]))
    np.testing.assert_allclose(geo1['pvt'].cartesian['velocity'].si.value, \
                               np.array([-2662.73375321, -1537.33004919, -0.]))
    np.testing.assert_allclose(np.array(geo1['pvt'].spherical.si.value.tolist()), \
                               np.array([2.0943951023931953, 0.0, 42164169.623295024, \
                                         7.292115855377073e-05, 0.0, 0.0]), \
                               atol=1.0e-7, rtol=0.0)
    np.testing.assert_allclose(ell1['pvt'].cartesian['position'].si.value, \
                               np.array([3698345.45792049, -34232447.31322606, -0.]))
    np.testing.assert_allclose(ell1['pvt'].cartesian['velocity'].si.value, \
                               np.array([2148.20304004, -1494.36501206, 0.]))
    np.testing.assert_allclose(np.array(ell1['pvt'].spherical.si.value.tolist()), \
                               np.array([4.82000783e+00, 0.00000000e+00, 3.44316454e+07, \
                                         5.73676739e-05, 0.00000000e+00, 1.71646077e+03]))
    np.testing.assert_allclose(ell2['pvt'].cartesian['position'].si.value, \
                               np.array([6538136.46,       0.  ,       0.  ]))
    np.testing.assert_allclose(ell2['pvt'].cartesian['velocity'].si.value, \
                               np.array([-0., 8695.15748257, 4721.08531441]))
    np.testing.assert_allclose(np.array(ell2['pvt'].spherical.si.value.tolist()), \
                               np.array([0.00000000e+00, 0.00000000e+00, 6.53813646e+06, \
                                         1.32991374e-03, 7.22084243e-04, 0.00000000e+00]))
    np.testing.assert_allclose(vang1['pvt'].cartesian['position'].si.value, \
                               np.array([  3313526.7287987 , -12405644.53079198,  -9013228.33893749]))
    np.testing.assert_allclose(vang1['pvt'].cartesian['velocity'].si.value, \
                               np.array([4321.62378008, -676.43293734, -491.45729632]))
    np.testing.assert_allclose(np.array(vang1['pvt'].spherical.si.value.tolist()), \
                               np.array([4.973394318345154, -0.6120236334818134, 15688140.76609011, \
                                         0.00031156788681119144, 3.91331269159631e-05, 1730.0341520321017]))
    np.testing.assert_allclose(gps1['pvt'].cartesian['position'].si.value, \
                               np.array([-15843456.17405487,  -2247776.5833332 ,  21200462.73723146]))
    np.testing.assert_allclose(gps1['pvt'].cartesian['velocity'].si.value, \
                               np.array([ 1454.40855181, -3518.76365779,   713.82704166]))
    np.testing.assert_allclose(np.array(gps1['pvt'].spherical.si.value.tolist()), \
                               np.array([3.28252622e+00, 9.24230194e-01, 2.65617624e+07, 2.30480399e-04, \
                                      4.46083005e-05, 1.11022302e-13]))
