import unittest

class TestElementTransform(unittest.TestCase):
    pos = [5740.13268349, 3314.06715   ,    0.]
    vel = [-2.75082684,  4.76457184,  5.50165367]
    pvt = tell.pvt((pos, vel, tell.dttm('2025-01-01T00:00:00')))

    def test_cart_to_kep(self):
        kep = tork.kepler(demoa.init.pvt)
