# Run with pytest -q test_astro.py

import numpy as np
import astropy.units as u
from munch import Munch

import tellurion.core as tell

# Create a structured quantity
#sqnorm = tell.makesq([12345.0, 245.0], ('alength','anangle'), ('km', 'deg'))

def test_angnorm():
    np.testing.assert_allclose(tell.normalizeangle(225*u.deg, u.rev/2), -135*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(45*u.deg, u.rev/2), 45*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(-45*u.deg, 1*u.rev), 315*u.deg)
    np.testing.assert_allclose(tell.normalizeangle(1135*u.deg, u.rev/2), 55*u.deg)
