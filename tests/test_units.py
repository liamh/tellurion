# Run with pytest -q test_astro.py

import datetime
import numpy as np
import astropy.units as u
from munch import Munch

import tellurion as tell

def test_angnorm():
    np.testing.assert_allclose(tell.normalizeangle(225*u.deg, u.rev/2), -135*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(45*u.deg, u.rev/2), 45*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(-45*u.deg, 1*u.rev), 315*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(1135*u.deg, u.rev/2), 55*u.deg)

# Test Kepler - this should give an error (inclination out of range)
# tell.kepler({"sma":8000.0, "ecc":0.1, "inc":-42.0, "argper":66.0, "raan":217.4, "ma":7.25}, tell.abstime('2023-09-14T08:30:00'))
# Test Kepler - this should normalize the raan
# tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25}, tell.abstime('2023-09-14T08:30:00'))
