import numpy as np
import astropy.units as u
import astro

# kep1 = kepler({"ecc":0.1, "sma":8000.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25})
def kepler(oes, dttm=None, units=(astro.prefunits['length'], astro.prefunits['angle'])):
    '''Make a Kepler orbital element set with either mean or true anomaly as the time element.'''
    if 'ma' in oes:
        keppt = {"sma":'length', "ecc":'dimensionless', "inc":'angle', "argper":'angle', "raan":'angle', "ma":'angle'}
    else:
        keppt = {"sma":'length', "ecc":'dimensionless', "inc":'angle', "argper":'angle', "raan":'angle', "ta":'angle'}
    ordoes = {k:oes[k] for k in keppt.keys()}
    kepsq = astro.makesq(ordoes, phystype=keppt)
    if dttm==None:
        return kepsq
    else:
        return (kepsq, dttm)
