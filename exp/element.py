import numpy as np
import astropy.units as u
import astro

def kepler(oes, dttm=None, units=(astro.prefunits['length'], astro.prefunits['angle'])):
    '''Make a Kepler orbital element set with either mean or true anomaly as the time element.'''
    if 'ma' in oes:
        keptype = [('sma', 'f8'), ('ecc', 'f8'), ('inc', 'f8'), ('argper', 'f8'), ('raan', 'f8'), ('ma', 'f8')]
    else:
        keptype = [('sma', 'f8'), ('ecc', 'f8'), ('inc', 'f8'), ('argper', 'f8'), ('raan', 'f8'), ('ta', 'f8')]
    npa = np.array((oes[keptype[0][0]], oes[keptype[1][0]], oes[keptype[2][0]], \
                    oes[keptype[3][0]], oes[keptype[4][0]], oes[keptype[5][0]]), dtype = keptype)
    kep = u.Quantity(npa, u.StructuredUnit((units[0], u.dimensionless_unscaled, units[1], \
                                            units[1], units[1], units[1])))
    if dttm==None:
        return kep
    else:
        return (kep, dttm)
