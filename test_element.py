# Run with pytest -q test_element.py

import numpy as np
import astropy.units as u
from munch import Munch
import tellurion.core as tell
import tellurion.ork as tork

def test_allplane():
    byalts = tell.kepler({"altper":160*u.km, "altapo":20250*u.km, \
                          "inc":28.5*u.deg, "argper": 0.0*u.deg, \
                          "raan": 0.0*u.deg, "ma": 0.0*u.deg}, \
                          tell.dttm('2022-02-15T08:30:00'))
    smaecc = tork.allplane(byalts)
    bysmaecc = tell.kepler({"sma":8000, "ecc":0.1, \
                          "inc":45, "argper": 120.0, "raan": 80.0, "ma": 0.0}, \
                          tell.dttm('2025-02-01T12:30:00'))
    alts = tork.allplane(bysmaecc)
    geo = tork.allplane(tell.kepler({"memo":1.0*u.rev/u.sday, \
                                     "ecc":0.0*u.dimensionless_unscaled, \
                                     "inc":0.0*u.deg, "argper": 120.0*u.deg, \
                                     "raan": 0.0*u.deg, "ma": 0.0*u.deg}))

    np.testing.assert_allclose(smaecc[0]['sma'], 16583.13646*u.km)
    np.testing.assert_allclose(smaecc[0]['ecc'], 0.60573583*u.dimensionless_unscaled)
    np.testing.assert_allclose(alts[0]['altper'], 821.86354*u.km)
    np.testing.assert_allclose(alts[0]['altapo'], 2421.86354*u.km)
    np.testing.assert_allclose(geo['sma'], 42164.1696233*u.km)
