# Run with pytest -q test_element.py

import numpy as np
import astropy.units as u
import tellurion.core as tell
import tellurion.ork as tork

def test_allplane():
    smaecc = tork.allplane({"altper":160*u.km, "altapo":20250*u.km, \
              "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    alts = tork.allplane({"sma":8000*u.km, "ecc":0.1*u.dimensionless_unscaled, \
                          "inc":45*u.deg, "argper": 120.0*u.deg, "raan": 80.0*u.deg, \
                          "ma": 0.0*u.deg})
    geo = tork.allplane({"memo":1.0*u.rev/u.sday, \
                         "ecc":0.0*u.dimensionless_unscaled, \
                         "inc":0.0*u.deg, "argper": 120.0*u.deg, \
                         "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    gto = tork.allplane({"altper": 350*u.km, "altapo": tork.sma(1.0,True), \
                         "ecc":0.0*u.dimensionless_unscaled, \
                         "inc":0.0*u.deg, "argper": 120.0*u.deg, \
                         "raan": 0.0*u.deg, "ma": 0.0*u.deg})
    ell = tork.allplane({"altper":160*u.km, "altapo":20250*u.km, \
                         "inc":28.5*u.deg, "argper": 0.0*u.deg, "raan": 0.0*u.deg, "ma": 0.0*u.deg})


    np.testing.assert_allclose(smaecc['sma'], 16583.13646*u.km)
    np.testing.assert_allclose(smaecc['ecc'], 0.60573583*u.dimensionless_unscaled)
    np.testing.assert_allclose(alts['altper'], 821.86354*u.km)
    np.testing.assert_allclose(alts['altapo'], 2421.86354*u.km)
    np.testing.assert_allclose(geo['sma'], 42164.1696233*u.km)
    # Test conversion to semimajor axis
    np.testing.assert_allclose(tork.sma(1e4*u.s), 10032.11910363*u.km)
    np.testing.assert_allclose(tork.sma(-20*(u.km/u.s)**2), 9965.0110375*u.km)
    np.testing.assert_allclose(gto['sma'], 24446.15304165*u.km)
