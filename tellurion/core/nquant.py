"""Make a dictionary (`namedquant`) of named quantities or structured quantity (`structquant`)"""

import string
import numpy as np
import astropy.units as u

# Extract value as np.array from quantities
u.Quantity.to_array = lambda self: np.array(self.value.tolist()) if self.isscalar \
    else self.value.view(np.float64).reshape(self.value.shape + (-1,))

# Make the structured quantity a singleton vector if it is a scalar
u.Quantity.tovector = lambda self: u.Quantity([self]) if self.isscalar else self

# Concatenate rows of the same structures
u.Quantity.vstack = lambda self, second: (self.vstack(second[0]).vstack(second[1:]) if len(second)>1 \
                                          else self.vstack(second[0])) if type(second) is list \
                                          else np.hstack((self.tovector(), second.tovector()))

# Now
# values.ndim = 1 => scalar value per name
# values.ndim = 2 => vector value per name
# Need option
# values.ndim = 2 => vector of scalar value per name, i.e., len(pv) exists
# Need
# values.ndim = 3 => vector of vectors as in pv example https://docs.astropy.org/en/stable/units/structured_units.html, take as input
# np.array([[[1., 0., 0.], [ 0.,  0.125,  0.]], [[0., 1., 0.], [-0.125,  0.,  0.]]])

def namedquant(values, names=string.ascii_letters, units=None, unitlookup={}):
    """Make a dictionary of quantities by name from numbers or quantities (including structured quantities).

    If units are provided either in the values or via the units
    argument, and unitlookup is not None, then the quantities will be
    converted to the units given in unitlookup. If the values all have
    units and units or unitlookup are provided, units will be
    converted.

    Parameters
    ----------
    values : A dict, list, or np.array of values; each value may be a number or u.Quantity
    names : the names for each part
    units : the units the structure; size must match number of values, optional
    unitlookup : a dict of physical type strings and preferred units, defaults to `prefunits`, optional

    Returns
    -------
    A dictionary with names as keys and `u.Quantity` as values

    Examples
    --------
    tell.namedquant([[1,2,3],[4,5,6]], ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.namedquant(np.array([[1,2,3],[4,5,6]]), ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.namedquant([12345.0, 45.0], ['sma','inc'], ['km', 'deg'])
    tell.namedquant([12345.0*u.km, 45.0*u.deg], ['sma','inc'], ['meter', 'radian'])
    tell.namedquant({'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, units = {'sma': 'meter', 'inc': 'radian'})
    tell.namedquant([{'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, {'sma' : 23456.0*u.km, 'inc' : -45.0*u.deg}], units = {'sma': 'meter', 'inc': 'radian'})
    """
    # if units and unitlookup:
    #     # Make the quantities and then convert
    #     return namedquant(namedquant(values, units=units, phystype=None, unitlookup=None), \
    #                       units=None, phystype=phystype, unitlookup=unitlookup)
    # el
    if type(values) is u.Quantity and type(values.unit) is u.StructuredUnit:
        if values.isscalar:
            return _splitsq(values)
        else:
            return [_splitsq(v) for v in values]
    else:
        (rvals, runits, rnames, vect) = _namedquant_tuple(values, names, units, unitlookup)
        if vect:
            return [{name: value for name, value in zip(rnames, [val*u.Unit(unit) for val, unit in zip(rv, runits)])} \
                    for rv in rvals]
        else:
            vals = zip(rnames, [val*u.Unit(unit) for val, unit in zip(rvals, runits)])
            return {name: value for name, value in vals}

def structquant(values, names=string.ascii_letters, units=None, unitlookup={}):
    """Make a structured quantity from numbers or quantities

    If units are provided either in the values or via the units
    argument, and unitlookup is not None, then the quantities will be
    converted to the units given in unitlookup. If the values all have
    units and units or unitlookup are provided, units will be
    converted.

    Parameters
    ----------
    values : A dict, list, or np.array of values; each value may be a number or u.Quantity
    names : the names for each part
    units : the units the structure; size must match number of values, optional
    unitlookup : a dict of physical type strings and preferred units, defaults to `prefunits`, optional

    Returns
    -------
    A structured `u.Quantity`

    Examples
    --------
    tell.structquant([[1,2,3],[4,5,6]], ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.structquant([[[1,2,3],[4,5,6]], [[-1,-2,-3],[-4,-5,-6]]], ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.structquant([[[1,2,3],6,[4,5,6]], [[10,20,30],60,[40,50,60]]], ['pos','sum','vel'], ['length', 'length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.structquant(np.array([[1,2,3],[4,5,6]]), ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})
    tell.structquant([12345.0, 45.0], ['sma','inc'], ['km', 'deg'])
    tell.structquant([12345.0*u.km, 45.0*u.deg], ['sma','inc'], ['meter', 'radian'])
    tell.structquant({'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, units = {'sma': 'meter', 'inc': 'radian'})
    tell.structquant([{'sma' : 12345.0*u.km, 'inc' : 45.0*u.deg}, {'sma' : 23456.0*u.km, 'inc' : -45.0*u.deg}], units = {'sma': 'meter', 'inc': 'radian'})
    """
    # tell.structquant(np.array([[[1,2,3],[4,5,6]], [[10,20,30],[40,50,60]]]), ['pos','vel'], ['length', 'speed'], {'length': 'km', 'speed': 'km/s'})



    # if units and unitlookup:
    #     return structquant(namedquant(values, units=units, phystype=phystype, unitlookup=unitlookup), \
    #                        unitlookup=unitlookup)

    (rvals, runits, rnames, vect) = _namedquant_tuple(values, names, units, unitlookup)

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
    if vect:
        sizes = [pqlen(v) for v in rvals[0]]
    else:
        sizes = [pqlen(v) for v in rvals]
    dtype = [(n, scvec(s)) for (n, s) in zip(rnames, sizes)]
    if vect:
        npa = np.array([tuple(v) for v in rvals], dtype=dtype)
    else:
        npa = np.array(tuple(rvals), dtype=dtype)
    return u.Quantity(npa, u.StructuredUnit(runits))

def depth(L): return isinstance(L, list) and max(map(depth, L))+1

def _namedquant_tuple(values, names=string.ascii_letters, units=None, unitlookup={}):
    """Return a tuple (values, units, names) where values are numbers or u.Quantity,
    units are names of units corresponding to the values, and `names`
    are names of the quantities.

    If any of the `values` are u.Quantity, then the ones that are
    numbers are made into u.Quantity with unit
    u.dimensionless_unscaled. The values will be converted to new
    units based on `units` and `unitlookup`.

    If all the `values` are numbers, then `units` and `unitlookup`
    determine what units they get.

    Parameters
    ----------
    values : dictionary, list, or array. May contain numbers, or u.Quantity (unless array).
      If a dict is given, the keys are the names of the structure elements, and `names` is ignored.
      If numbers are given for values, `units` or `unitlookup` must be specified
      If `u.Quantity` are given for values, they are converted using
      `units` or `unitlookup` if either (or both) are specified.
      >>>> Each value can be scalars or vectors
    names : the names of the fields or keys, required if `values` is a list
    units : the units or physical types; must be a dictionary with the
      same keys or list with the same length as values; optional
    unitlookup : a dict of physical type strings and preferred units, optional
      If units contains the name of a physical type when
      In order for a physical type to be recognized in unitlookup when
      it needs to be inferred from `values`, it must be the first named
      physical type returned by `._physicaltype_list()`,
      e.g. `energy density` instead of 'pressure':
        u.get_physical_type('pressure')._physical_type_list
        ['energy density', 'pressure', 'stress']

    Returns
    -------
    Tuple of (values, units, names)

    """

    # Extract names and values from dict if needed
    if type(values) is dict:
        names = tuple(values.keys())
        vals = tuple(values.values())
        vect = False
    elif type(values) is list and all(isinstance(item, dict) for item in values) \
         and all(set(d.keys()) == set(values[0].keys()) for d in values):
        # Make a non-scalar structured/named quantity from a list of dicts
        names = tuple(values[0].keys())
        vals = [tuple(w.values()) for w in values]
        vect = True
    else: # values is not a dict or list of dicts
        # [[1,2],[3,4],[5,6]] # could be vector of scalars or scalar of vectors
        # It will always be interpreted as a scalar of vectors
        vals = values
        vect = depth(values)==3
    if vect:
        vals0 = vals[0]
    else:
        vals0 = vals

    # If there are any u.Quantity in values, all values without units are assumed to be u.dimensionless.
    if any([hasattr(v, 'unit') for v in vals0]):
        if vect:
            vals = [[v*u.dimensionless_unscaled for v in w] for w in vals]
        else:
            vals = [v*u.dimensionless_unscaled for v in vals]
        units_in_values = [v.unit.to_string() for v in vals0]
        phtys_in_values = [u.get_physical_type(v.unit)._physical_type_list[0] for v in vals0]
    else:
        units_in_values = None

    def lookuppt(un, default=None):
        return unitlookup.get(un) or un or default

    # If names of units correspond to a physical type in `unitlookup`, substitute those names
    newunits = units_in_values
    if units_in_values: # Convert the u.Quantity given
        if type(units) is dict:
            newunits = tuple([lookuppt(units.get(nm), val.unit.to_string()) for nm, val in zip(names, vals0)])
        elif type(units) is list:
            newunits = tuple([lookuppt(un, val.unit.to_string()) for un, val in zip(units, vals0)])
        else:
            newunits = tuple([unitlookup.get(pt) or un for un, pt in zip(units_in_values, phtys_in_values)])
        if vect:
            vals = [tuple([val.to(u.Unit(un)).value for val, un in zip(w, newunits)]) for w in vals]
        else:
            vals = tuple([val.to(u.Unit(un)).value for val, un in zip(vals, newunits)])
    else:  # Assign units for numbers
        if type(units) is dict:
            newunits = tuple([lookuppt(units.get(nm)) for nm in names])
        elif type(units) is list:
            newunits = tuple([lookuppt(un) for un in units])
        else:
            newunits = tuple(['1']*len(vals0))
    return (vals, newunits, names[0:len(vals0)], vect)

# generate_variations_nquant({'alength':1000.0*u.km, 'amass':200.0*u.kg, 'atime':3.0*u.hour}, {'alength':'m', 'amass':'g', 'atime':'s'})
def generate_variations_nquant(dict_with_units, new_units, \
                               unitlookup={'length':'meter', 'mass':'gram', 'time': 'second'}, \
                               fn=_namedquant_tuple):
    """Generate and apply variations of arguments for testing"""
    names = list(dict_with_units.keys())
    dict_no_units = dict((k, v.value) for k, v in dict_with_units.items())
    dict_units = dict((k, v.unit.to_string()) for k, v in dict_with_units.items())
    dict_phystype = dict((k, list(u.get_physical_type(v)._physical_type)[0]) \
                         for k, v in dict_with_units.items())
    dict_ul = dict((list(u.get_physical_type(v)._physical_type)[0], un) \
                         for v, un in zip(dict_with_units.values(), list(dict_units.values())))
    list_no_units = list(dict_no_units.values())
    array_no_units = np.array(list_no_units)
    list_units = list(dict_units.values())
    list_new_units = list(new_units.values())
    list_phystype = list(dict_phystype.values())
    list_with_units = list(dict_with_units.values())
    argsets = [
        # values are numbers; units if given are assigned
        {'values': dict_no_units}, # as dict, no units
        {'values': dict_no_units, 'units': dict_units}, # as dict, units in argument `units`
        {'values': dict_no_units, 'units': dict_phystype, 'unitlookup': dict_ul}, # as dict, physical type in argument `units`
        {'values': list_no_units, 'names': names}, # as list, no units
        {'values': list_no_units, 'names': names, 'units': list_units},  # as list, units in argument `units`
        {'values': list_no_units, 'names': names, 'units': list_units, 'unitlookup': dict_ul},  # as list, physical type in argument `units`
        {'values': array_no_units, 'names': names}, # as array, no units
        {'values': array_no_units, 'names': names, 'units': list_units},  # as array, units in argument `units`
        {'values': array_no_units, 'names': names, 'units': list_units, 'unitlookup': dict_ul},  # as array, physical type in argument `units`

        # values are u.Quantity; units are as given in values or converted
        {'values': dict_with_units}, # as dict, no conversion
        {'values': dict_with_units, 'units': new_units}, # as dict, convert to units in `new_units`
        {'values': dict_with_units, 'unitlookup': unitlookup}, # as dict, convert units, by inferred physical type
        {'values': dict_with_units, 'units': dict_phystype, 'unitlookup': unitlookup},  # as dict, convert units by explicit physical type
        {'values': dict_with_units, 'units': dict_units},  # as dict, "convert" to original units, explicit
        {'values': dict_with_units, 'unitlookup': dict_ul}, # as dict, "convert" to original units, inferred
        {'values': dict_with_units, 'units': dict_phystype, 'unitlookup': dict_ul}, # as dict, "convert" to original units by physical type
        {'values': list_with_units, 'names': names}, # as list, no conversion
        {'values': list_with_units, 'names': names, 'units': list_new_units}, # as list, convert to units in `new_units`
        {'values': list_with_units, 'names': names, 'unitlookup': unitlookup}, # as dict, convert units
        {'values': list_with_units, 'names': names, 'units': list_phystype, 'unitlookup': unitlookup}, # as dict, convert explicitly stated units
        {'values': list_with_units, 'names': names, 'units': list_units},  # as list, "convert" to original units
        {'values': list_with_units, 'names': names, 'units': list_units, 'unitlookup': dict_ul},  # as list, physical type in argument `units`
        {'values': list_with_units, 'names': names, 'units': list_phystype, 'unitlookup': dict_ul}, # as list, "convert" to original units by physical type
        ]
    return (argsets, [fn(**a) for a in argsets])


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
