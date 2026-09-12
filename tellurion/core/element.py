"""Orbital elements user interface."""

import astropy.units as u
import numpy as np
from astropy.time import Time

from tellurion.astro import quantity_utils as quant, units as tunits

# ----------------------------
# State function values
# ----------------------------


def statefnval(orbstate, quantname, sfdict, constants=None):
    """
    Compute value(s) given by function(s) of the orbital state.

    Arguments
      orbstate:   Representation of orbital state in any form
      quantname:  Name(s) of quantity; may be a list, e.g. ["sma", "ecc"]
      constants:  Quantities independent of the orbital state

    """
    if isinstance(quantname, list):
        return [statefnval(orbstate, itm, sfdict, constants)
                for itm in quantname]
    else:
        lookup = sfdict[quantname]
        getter = lookup["getter"]
        if "__code__" in dir(getter) and len(getter.__code__.co_varnames) > 1:
            orkval = getter(orbstate, constants)
        else:
            orkval = getter(orbstate)
        orkunit = lookup["orkunit"]
        return u.Quantity(orkval, orkunit).to(tunits.prefunits[lookup["phystype"]])


def sfdict(sfvbl):
    """Make a state function dictionary of the state function variables."""
    keys = ["name", "description", "phystype", "orkunit", "getter"]
    return dict(zip([ev[0] for ev in sfvbl],
                    [dict(zip(keys, ev)) for ev in sfvbl]))


# ----------------------------
#  Element sets
# ----------------------------


kepeltma_names = ["ecc", "sma", "inc", "argper", "raan", "ma"]
kepeltta_names = ["ecc", "sma", "inc", "argper", "raan", "ta"]
timeelements = ["ta", "ma"]


class ElementSetT:
    """
    An orbital element set paired with an epoch time.

    Parameters
    ----------
    elements : `~astropy.units.Quantity`
        Structured Quantity array of orbital elements (e.g. sma, ecc, inc, …).
    time : `~astropy.time.Time`
        The epoch associated with the element set.
    """

    def __init__(self, elements, time):
        if not isinstance(time, Time):
            raise TypeError(f"time must be an astropy Time, got {type(time)}")
        self.elements = elements
        self.time = time

    def __repr__(self):
        return f"ElementSetT(elements={self.elements!r}, time={self.time!r})"

    def __eq__(self, other):
        if not isinstance(other, ElementSetT):
            return NotImplemented
        return (np.all(self.elements == other.elements) and self.time == other.time)

    def __iter__(self):
        yield self.elements
        yield self.time

    def __getitem__(self, idx):
        return (self.elements, self.time)[idx]

    def pvt(self):
        """
        Convert this element set to a Cartesian state.

        Returns
        -------
        tellurion.core.posvel.PositionVelocityT
            Cartesian state at the same epoch as this element set.

        Notes
        -----
        This is the instance-method form of :func:`tellurion.pvt` for callers
        who already have an :class:`ElementSetT`.
        Unlike :class:`~tellurion.core.posvel.PositionVelocityT` conversion
        methods, this method does not accept ``mean_time_element`` because the
        element set already encodes whether it uses mean or true angular
        elements.
        This is an Orekit-backed conversion. Calling it may trigger lazy
        initialization of the JVM and Orekit data on first use.
        """
        from importlib import import_module

        ork_element = import_module("tellurion.ork.element")
        return ork_element.pvt(self)

# kep1 = kepler({"ecc":0.1, "sma":8000.0, "inc":42.0, "argper":66.0, "raan":217.4, "ma":7.25})
# kep2 = kepler({"zper":250.0, "zapo":350.0, "inc":22.0, "argper":66.0, "raan":68.0, "ma":7.25})
def kepler(oes, dttm=None, unitlookup=tunits.prefunits):
    """Make a Kepler orbital element set with either mean or true
    anomaly as the time element. Units not specified default to those
    given in unitlookup."""

    keppt = {"inc":"angle", "argper":"angle", "raan":"angle", "ma":"angle", "ta":"angle", \
             "ecc":"dimensionless", "sma":"length", "memo":"angular speed", \
             "radper":"length", "radapo":"length", "altper":"length", "altapo":"length"}
    isscalar = not(hasattr(dttm, "isscalar")) or dttm.isscalar
    kepsq = quant.make_quantity(oes, keppt, isscalar, unitlookup)
    # Possibly check kepsq['inc'] is in upper halfplane with
    # .is_within_bounds('0d', '180d') on Angle instances
    #    raise ValueError('Inclination must be between 0 and 180 degrees, inclusive')
    kepsqn = tunits.normalizeangle(kepsq, u.rev/2, timeelements)
    if dttm is None:
        return kepsqn
    else:
        return ElementSetT(elements=kepsqn, time=dttm)

def iskepels(obj, est=True):
    if type(obj) is ElementSetT and est:
        return iskepels(obj.elements, False) and type(obj.time) is Time
    else:
        return type(obj) is u.Quantity \
            and obj.dtype.names is not None \
            and (not(set(kepeltma_names) - set(obj.dtype.names)) \
                 or not(set(kepeltta_names) - set(obj.dtype.names)))

# --- Equinoctial elements ---
# a, ex=e·cos(ω+Ω), ey=e·sin(ω+Ω), hx=tan(i/2)·cos(Ω), hy=tan(i/2)·sin(Ω), λ
equeltma_names = ["sma", "ex", "ey", "hx", "hy", "ml"]   # mean longitude
equeltta_names = ["sma", "ex", "ey", "hx", "hy", "tl"]   # true longitude
equtimeelements = ["ml", "tl"]

# --- Circular elements ---
# a, ex=e·cos(αω), ey=e·sin(αω), i, Ω, u (latitude argument)
circeltma_names = ["sma", "cex", "cey", "inc", "raan", "mla"]  # mean latitude arg
circeltta_names = ["sma", "cex", "cey", "inc", "raan", "tla"]  # true latitude arg
circtimeelements = ["mla", "tla"]

_equpt = {
    "sma":  "length",
    "ex":   "dimensionless", "ey":   "dimensionless",
    "hx":   "dimensionless", "hy":   "dimensionless",
    "ml":   "angle",         "tl":   "angle",
}

_circpt = {
    "sma":  "length",
    "cex":  "dimensionless", "cey":  "dimensionless",
    "inc":  "angle",
    "raan": "angle",
    "mla":  "angle",         "tla":  "angle",
}


def equinoctial(oes, dttm=None, unitlookup=tunits.prefunits):
    """Make an equinoctial orbital element set.

    Parameters
    ----------
    oes : dict or structured Quantity
        Keys: sma, ex, ey, hx, hy, and one of ml (mean longitude) or tl
        (true longitude).
    dttm : `~astropy.time.Time`, optional
        Epoch; if supplied, returns an ``ElementSetT``.
    """
    isscalar = not hasattr(dttm, "isscalar") or dttm.isscalar
    eqsq = quant.make_quantity(oes, _equpt, isscalar, unitlookup)
    eqsqn = tunits.normalizeangle(eqsq, u.rev / 2, equtimeelements)
    return ElementSetT(elements=eqsqn, time=dttm) if dttm is not None else eqsqn

def isequels(obj, est=True):
    if type(obj) is ElementSetT and est:
        return isequels(obj.elements, False) and type(obj.time) is Time
    return (type(obj) is u.Quantity
            and obj.dtype.names is not None
            and (not set(equeltma_names) - set(obj.dtype.names)
                 or not set(equeltta_names) - set(obj.dtype.names)))

def circular(oes, dttm=None, unitlookup=tunits.prefunits):
    """Make a circular orbital element set.

    Parameters
    ----------
    oes : dict or structured Quantity
        Keys: sma, cex, cey, inc, raan, and one of mla (mean latitude
        argument) or tla (true latitude argument).
    dttm : `~astropy.time.Time`, optional
        Epoch; if supplied, returns an ``ElementSetT``.
    """
    isscalar = not hasattr(dttm, "isscalar") or dttm.isscalar
    csq = quant.make_quantity(oes, _circpt, isscalar, unitlookup)
    csqn = tunits.normalizeangle(csq, u.rev / 2, circtimeelements)
    return ElementSetT(elements=csqn, time=dttm) if dttm is not None else csqn

def iscircels(obj, est=True):
    if type(obj) is ElementSetT and est:
        return iscircels(obj.elements, False) and type(obj.time) is Time
    return (type(obj) is u.Quantity
            and obj.dtype.names is not None
            and (not set(circeltma_names) - set(obj.dtype.names)
                 or not set(circeltta_names) - set(obj.dtype.names)))
