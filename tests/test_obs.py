import numpy as np
import astropy.units as u
import astropy.coordinates as coord
import tellurion as tell

newyear = tell.abstime('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour
kickapoo = tell.earthloc(lon='98°45′49.82″W', lat='33°33′08.50″N')

def test_siderealtime():
    np.testing.assert_allclose(tell.siderealtime(newyear), coord.Longitude(100.89974637, u.deg))
    np.testing.assert_allclose(tell.siderealtime(newyear + prop5m1h), \
                               coord.Longitude([102.15316876, 103.40659114, 104.66001353, 105.91343592, \
                                                107.16685831, 108.4202807 , 109.67370309, 110.92712548, \
                                                112.18054787, 113.43397026, 114.68739264, 115.94081503], u.deg))
    np.testing.assert_allclose(tell.siderealtime(newyear + prop5m1h, kickapoo), \
                               coord.Longitude([200.91700764, 202.17043003, 203.42385242, 204.67727481, \
                                                205.9306972 , 207.18411959, 208.43754198, 209.69096437, \
                                                210.94438676, 212.19780914, 213.45123153, 214.70465392], u.deg))
