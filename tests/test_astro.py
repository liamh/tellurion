# Run with pytest -q test_astro.py

import datetime
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

def test_abstime():
    newyear = tell.abstime('2025-01-01T00:00:00')
    prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
    np.testing.assert_allclose(newyear.to_value('jd'), 2460676.5)
    np.testing.assert_allclose(tell.abstime(datetime.datetime(2025, 4, 17, 22, 54, 8, 684006)).to_value('jd'), \
                               2460783.454267176)
    np.testing.assert_allclose(tell.abstime(['2025-01-01T00:00:00', '2025-01-02T00:00:00', \
                                             '2025-01-03T00:00:00']).to_value('jd'), \
                               np.array([2460676.5, 2460677.5, 2460678.5]))
    np.testing.assert_allclose(tell.abstime(prop5m1h, newyear).to_value('jd'), \
                               np.array([2460676.50347222, 2460676.50694444, 2460676.51041667, \
                                         2460676.51388889, 2460676.51736111, 2460676.52083333, \
                                         2460676.52430556, 2460676.52777778, 2460676.53125, \
                                         2460676.53472222, 2460676.53819444, 2460676.54166667]))
    np.testing.assert_allclose(tell.abstime(5*u.hour, newyear).to_value('jd'), 2460676.7083333335)
    np.testing.assert_allclose(tell.abstime('12d 17hr 23min 33.1s', newyear).to_value('jd'), \
                               2460689.2246886576)
    np.testing.assert_allclose(tell.abstime([5*u.hour, '12d 17hr 23min 33.1s'], newyear).to_value('jd'), \
                               np.array([2460676.7083333335, 2460689.2246886576]))


# Test Kepler - this should give an error (inclination out of range)
# tell.kepler({"sma":8000.0, "ecc":0.1, "inc":-42.0, "argper":66.0, "raan":217.4, "ma":7.25}, tell.abstime('2023-09-14T08:30:00'))
# Test Kepler - this should normalize the raan
# tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25}, tell.abstime('2023-09-14T08:30:00'))
