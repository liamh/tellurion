import collections
import numpy as np
import astropy.units as u
import astropy.constants # astropy.constants.R_earth
from astropy.time import Time
from astropy.timeseries import TimeSeries
from . import astro
from . import nquant
from . import posvel

###############################
#### State function values ####
###############################

def statefnval (orbstate, quantname, sfdict, constants=None):
    """
    Compute value(s) given by function(s) of the orbital state
    Arguments
      orbstate:   Representation of orbital state in any form
      quantname:  Name(s) of quantity; may be a list, e.g. ["sma", "ecc"]
      constants:  Quantities independent of the orbital state
    """
    if isinstance(quantname, list):
        return [statefnval(orbstate, itm, sfdict, constants) for itm in quantname]
    else:
        lookup = sfdict[quantname]
        getter = lookup["getter"]
        if '__code__' in dir(getter) and len(getter.__code__.co_varnames) > 1:
            orkval = getter(orbstate, constants)
        else:
            orkval = getter(orbstate)
        orkunit = lookup["orkunit"]
        return u.Quantity(orkval, orkunit).to(astro.prefunits[lookup["phystype"]])

def sfdict(sfvbl):
    '''Make a state function dictionary of the state function variables'''
    keys = ["name", "description", "phystype", "orkunit", "getter"]
    return dict(zip([ev[0] for ev in sfvbl], [dict(zip(keys,ev)) for ev in sfvbl]))

##############################
####  Element sets        ####
##############################

ElementSetT = collections.namedtuple('ElementSetT', 'els t')

kepeltma_names = ["ecc", "sma", "inc", "argper", "raan", "ma"]
kepeltta_names = ["ecc", "sma", "inc", "argper", "raan", "ta"]
timeelements=['ta', 'ma']

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
    elif 'memo' in oes and 'ecc' in oes:
        plane = {"ecc":'dimensionless', "memo":'angular speed'}
    elif 'radper' in oes and 'radapo' in oes:
        plane = {"radper":'length', "radapo":'length'}
    elif 'altper' in oes and 'altapo' in oes:
        plane = {"altper":'length', "altapo":'length'}

    keppt = {"inc":'angle', "argper":'angle', "raan":'angle'} | plane | timeelt
    ordoes = {k:oes[k] for k in keppt.keys()}
    kepsq = nquant.structquant(ordoes, phystype=keppt)
    if not astro.isupperhalfplane(kepsq['inc']):
        raise ValueError('Inclination must be between 0 and 180 degrees, inclusive')
    kepsqn = astro.normalizeangle(kepsq, u.rev/2, timeelements)
    if dttm==None:
        return kepsqn
    else:
        return ElementSetT(kepsqn, dttm)

def iskepels(obj, est=True):
    '''The argument is a Keper element set or (kepels, epoch); if
    `est` is `True`, then it is a properly constructed `ElementSetT`,
    and if `False`, it is element values only, without an epoch time.'''
    if type(obj) is ElementSetT and est:
        return iskepels(obj.els, False) and posvel.isdttm(obj.t)
    else:
        return type(obj) is u.Quantity \
            and (not(set(kepeltma_names) - set(obj.dtype.names)) \
                 or not(set(kepeltta_names) - set(obj.dtype.names)))
