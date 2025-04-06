import numpy as np
import astropy.units as u
from munch import Munch

import tellurion.core as tell
import tellurion.ork as tork

################################
####  Comparison of states  ####
################################

def pvtequal(a, b):
    np.testing.assert_allclose(a.pv['position'],b.pv['position'], rtol=1e-5, atol=1e-12)
    np.testing.assert_allclose(a.pv['velocity'],b.pv['velocity'], rtol=1e-5, atol=1e-12)
    np.testing.assert_equal(a[1], b[1])

def kepequal(a, b):
    for nm in tell.kepeltma_names:
        np.testing.assert_allclose(a[0][nm],b[0][nm])
    np.testing.assert_equal(a[1], b[1])

##########################
####   Definitions    ####
##########################

newyear = tell.dttm('2025-01-01T00:00:00')
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

def test_kepcart():
    pvtequal(demoa.init.cart, demoa.init.pvt)

demob = Munch()
demob.kep = tell.kepler({"sma":8000.0, "ecc":0.1, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25},
                           tell.dttm('2023-09-14T08:30:00'))
demob.pvt = tork.cartesian(*demob.kep) # Convert Kepler elements to PVT
demob.rekep = tork.kepler(demob.pvt)

def test_cartkep():
    kepequal(demob.rekep, demob.kep)
