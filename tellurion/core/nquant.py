"""From a set of named quantities or structured quantity, make a dictionary (`namedquant`) or structured quantity (`structquant`)"""

import numpy as np
import astropy.units as u
from tellurion.core import astro

# obsdict = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed'}
# nq = namedquant({'azim': 110.0, 'elev': 81.0, 'range': 662.1}, phystype=obsdict)
# sq = structquant(nq)
# nq2 = namedquant(sq)
# nq and nq2 are the same
#  {'azim': <Quantity 110. deg>,
#   'elev': <Quantity 81. deg>,
#   'range': <Quantity 662.1 km>}
# sq is the structured quantity
#  <Quantity (110., 81., 662.1) (deg, deg, km)>
# sq.value.dtype.names
#  ('azim', 'elev', 'range')

def namedquant(values, names=None, units=None, phystype=None, unitlookup=astro.prefunits):
    """Make a dictionary of quantities by name from numbers or quantities (including structured quantities).

    Parameters
    ----------
    values : dict of numbers, dict of `u.Quantity`, list or array of numbers, list or array of `u.Quantity`
      If a dict is given, the keys are the names of the structure elements, and `names` is ignored
      If numbers are given for values, `units` or `phystype` must be specified
      If `u.Quanitity` are given for values, they are converted to `units` or `phystype` and `unitlookup`
      Each value can be scalars or vectors
    names : the names of the fields of the structure;
      size must match number of values (ignored if `values` is a dict), optional
    units : the units the structure; size must match number of values, optional
    phystype : the physical type (dimension) of the each argument, units are looked
      up in `unitlookup`, may specified instead of `units`, optional
      Must be a dict if `values` is a dict
    unitlookup : a dict of physical type strings and preferred units, defaults to `prefunits`, optional

    Returns
    -------
    A dictionary with names as keys and `u.Quantity` as values

    Examples
    --------
    >>> namedquant([[1,2,3],[4,5,6]], ('pos','vel'), phystype=('length', 'speed'))
    {'pos': <Quantity [1., 2., 3.] km>, 'vel': <Quantity [4., 5., 6.] km / s>}

    >>> namedquant([12345.0, 45.0], ('sma','inc'), phystype=('length', 'angle'))
    {'sma': <Quantity 12345. km>, 'inc': <Quantity 45. deg>}

    >>> namedquant([12345.0, 45.0], ('sma','inc'), ('km', 'deg'))
    {'sma': <Quantity 12345. km>, 'inc': <Quantity 45. deg>}

    >>> namedquant([12345.0*u.km, 45.0*u.deg], ('sma','inc'), units = ('meter', 'radian'))
    {'sma': <Quantity 12345000. m>, 'inc': <Quantity 0.78539816 rad>}

    >>> namedquant([12345000*u.m, 0.125*u.rev], ('sma','inc'), phystype = ('length', 'angle'))
    {'sma': <Quantity 12345. km>, 'inc': <Quantity 45. deg>}

    >>> namedquant({'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, units = ('meter', 'radian'))
    {'sma': <Quantity 12345000. m>, 'inc': <Quantity 0.78539816 rad>}

    """
    if type(values) is u.Quantity and type(values.unit) is u.StructuredUnit:
        return _splitsq(values)
    else:
        (converted_vals, units, names) = _namedquant_tuple(values, names, units, phystype, unitlookup)
        vals = zip(names, [val*u.Unit(unit) for val, unit in zip(converted_vals, units)])
        return {name: value for name, value in vals}

def structquant(values, names=None, units=None, phystype=None, unitlookup=astro.prefunits):
    """Make a structured quantity from numbers or quantities

    Parameters
    ----------
    values : dict of numbers, dict of `u.Quantity`, list or array of numbers, list or array of `u.Quantity`
      If a dict is given, the keys are the names of the structure elements, and `names` is ignored
      If numbers are given for values, `units` or `phystype` must be specified
      If `u.Quanitity` are given for values, they are converted to `units` or `phystype` and `unitlookup`
      Each value can be scalars or vectors
    names : the names of the fields of the structure;
      size must match number of values (ignored if `values` is a dict), optional
    units : the units the structure; size must match number of values, optional
    phystype : the physical type (dimension) of the each argument, units are looked
      up in `unitlookup`, may specified instead of `units`, optional
      Must be a dict if `values` is a dict
    unitlookup : a dict of physical type strings and preferred units, defaults to `prefunits`, optional

    Returns
    -------
    A structured `u.Quantity`

    Examples
    --------
    >>> structquant([[1,2,3],[4,5,6]], ('pos','vel'), phystype=('length', 'speed'))
    <Quantity ([1., 2., 3.], [4., 5., 6.]) (km, km / s)>

    >>> structquant([12345.0, 45.0], ('sma','inc'), phystype=('length', 'angle'))
    <Quantity (12345., 45.) (km, deg)>

    >>> structquant([12345.0, 45.0], ('sma','inc'), ('km', 'deg'))
    <Quantity (12345., 45.) (km, deg)>

    >>> structquant([12345.0*u.km, 45.0*u.deg], ('sma','inc'), units = ('meter', 'radian'))
    <Quantity (12345000., 0.78539816) (m, rad)>

    >>> structquant([12345000*u.m, 0.125*u.rev], ('sma','inc'), phystype = ('length', 'angle'))
    <Quantity (12345., 45.) (km, deg)>

    >>> structquant({'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, units = ('meter', 'radian'))
    <Quantity (12345000., 0.78539816) (m, rad)>

    >>> obsdict = {'azim': 'angle', 'elev': 'angle', 'range': 'length', 'rangerate': 'speed'}
    >>> structquant({'azim': 110.0, 'elev': 81.0, 'range': 662.1}, phystype=obsdict)
    <Quantity (110., 81., 662.1) (deg, deg, km)>

    """
    (converted_vals, units, names) = _namedquant_tuple(values, names, units, phystype, unitlookup)

    def scvec(size):
        if size==1:
            return f"f8"
        else:
            return f"({size},)f8"
    def pqlen(item):
        if type(item) is u.Quantity:
            if type(item.value) is np.ndarray:
                return len(item)
            else:
                return 1
        else:
            if type(item) is list or type(item) is tuple or type(item) is np.ndarray:
                return len(item)
            else:
                return 1

    # Build the structured array
    sizes = [pqlen(v) for v in converted_vals]
    dtype = [(n, scvec(s)) for (n, s) in zip(names, sizes)]
    npa = np.array(tuple(converted_vals), dtype=dtype)
    return u.Quantity(npa, u.StructuredUnit(units))

def _namedquant_tuple(values, names=None, units=None, phystype=None, unitlookup=astro.prefunits):
    def conv(val, unit):
        if type(val) is u.Quantity:
            return val.to(unit).value
        else:
            return val

    # Extract names and values from dict if needed
    if type(values) is dict:
        names = tuple(values.keys())
        vals = tuple(values.values())
    else:
        vals = values

    # Determine units for each value
    if units is None:
        if phystype is None:
            # Extract units from any quantities that have them
            units = []
            for v in vals:
                if type(v) is u.Quantity:
                    units.append(v.unit)
                else:
                    raise ValueError("If values are not Quantities, units or phystype must be specified")
            units = tuple(units)
        else:
            # Use phystype to look up units
            if type(values) is dict:
                units = tuple([unitlookup[phystype[nm]] for nm in names])
            else:
                units = tuple([unitlookup[pt] for pt in phystype])
    else:
        # Units explicitly provided
        if type(values) is dict and type(units) is dict:
            units = tuple(units[nm] for nm in names)
        # else units is already a tuple/list

    # Convert all values to the target units, handling both Quantities and raw values
    converted_vals = []
    for v, unit in zip(vals, units):
        converted_vals.append(conv(v, unit))

    return (converted_vals, units, names)

def _splitsq(stqu, quant=True, readably=False):
    '''Make a dict of names and quantities, or values, names, and units from the structured quantity. To recreate a structured quantity, set `quant` to False; this results in input for structquant. For example, `structquant(*_splitsq(sq, False))` copies `sq`.'''
    if quant:
        vnu = _splitsq(stqu, False)
        qs = [v*u for (v, u) in zip(vnu[0], vnu[2])]
        return {un:val for (val, un) in zip(qs, vnu[1])}
    elif readably:
        return (stqu.value.tolist(), stqu.dtype.names, tuple(un.to_string() for un in stqu.unit.values()))
    else:
        return (stqu.value.tolist(), stqu.dtype.names, stqu.unit.values())
