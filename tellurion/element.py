import numpy as np
import astropy.units as u
import astropy.constants # astropy.constants.R_earth
from astropy.time import Time
from astropy.timeseries import TimeSeries
from . import astro

# kep1 = kepler({"ecc":0.1, "sma":8000.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25})
# kep2 = kepler({"zper":250.0, "zapo":350.0, "inc":22.0, "argper":66.0, "raan":68.0, "ma":7.25})
def kepler(oes, dttm=None, units=(astro.prefunits['length'], astro.prefunits['angle'])):
    '''Make a Kepler orbital element set with either mean or true anomaly as the time element.'''
    if 'ma' in oes:
        timeelt = {"ma":'angle'}
    else:
        timeelt = {"ta":'angle'}

    if 'sma' in oes and 'ecc' in oes:
        plane = {"ecc":'dimensionless', "sma":'length'}
    elif 'rper' in oes and 'rapo' in oes:
        plane = {"rper":'length', "rapo":'length'}
    elif 'zper' in oes and 'zapo' in oes:
        plane = {"zper":'length', "zapo":'length'}

    keppt = {"inc":'angle', "argper":'angle', "raan":'angle'} | plane | timeelt
    ordoes = {k:oes[k] for k in keppt.keys()}
    kepsq = astro.makesq(ordoes, phystype=keppt)
    if dttm==None:
        return kepsq
    else:
        return (kepsq, dttm)

def iskepels(obj):
    '''The argument is a Keper element set or (kepels, epoch)'''
    if type(obj) is tuple and len(obj)==2 and type(obj[1])==Time:
        return iskepels(obj[0])
    else:
        knames = {'argper', 'ecc', 'inc', 'ma', 'raan', 'sma'}
        return type(obj) is u.Quantity and not(knames - set(obj.dtype.names))
