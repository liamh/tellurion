"""Make a dictionary (`namedquant`) of named quantities or structured quantity (`structquant`)"""

import numbers, string
import collections.abc
import numpy as np
from numpy.lib import recfunctions as rfn
import astropy.units as u

############################################################
## Make structured quantities from arrays and other sq
############################################################

################### Used only in cartesianpv()
def sq(values, names_size_units, unitlookup={}):
    """Make a struct quantity from arrays. Make a scalar struct if
    `values` is a 1d np.array, or a vector struct if it is a 2d
    np.array.  `names_size_units` is a sequence of 3-tuples (name size
    units).

    sq(np.array([[1,2,3,6,4,5,6], [10,20,30,60,40,50,60]]), \
                                    [('pos', 3, 'length'), ('sum', 1, 'length'), ('vel', 3, 'speed')], \
                                    {'length': 'km', 'speed': 'km/s'})

    """
    def lookuppt(un, default=None):
        return unitlookup.get(un) or un or default
    dtype = [(nsu[0], scvec(nsu[1])) for nsu in names_size_units]
    unit = tuple([lookuppt(nsu[2]) for nsu in names_size_units])
    return compose_sq(values, dtype, unit)

# Make the structured quantity a singleton vector if it is a scalar
u.Quantity.tovector = lambda self: u.Quantity([self]) if self.isscalar else self

# Extract value as np.array from quantities
u.Quantity.to_array = lambda self: rfn.structured_to_unstructured(self.value) if isinstance(self.unit, u.StructuredUnit) else self.value

################### Used only in spacetrack to serialize/deserialize tles

# Decompose SQ into an array, dtype, and unit string
# If the object is not a structured quantity, then `name` will be used for the name and if further it is a vector (1d array), isscalar is used to determine whether that should be interpreted as a scalar structure (row vector) or vector structure (column vector)
u.Quantity.decompose = lambda self, name=None, isscalar=False: (self.to_array(), self.dtype.descr, self.unit.to_string()) \
    if self.dtype.names else (self.value, (name, '<f8', (1,)) if isscalar else (name, '<f8'), self.unit.to_string())

def compose_sq(values, dtype, unit_string):
    """Make a structured quantity from an np.array, dtype, and unit string"""
    return rfn.unstructured_to_structured(np.array(values), dtype=dtype)*u.Unit(unit_string)

################# Used to see what the SQ consists of
# Make a dictionary by structure components
u.Quantity.to_dict = lambda self: {nm: self[nm] for nm in self.dtype.names}


################ Used only by spacetrack
def dict_decompose(d):
    return {k: v.decompose() for k, v in d.items()}

def dict_compose(d):
    return {k: v[0]*u.Unit(v[2]) for k, v in d.items()}

################ Used only by sphericalpv,
##### replace with vstack() below, but that assumes their already pq and don't need to be converted
def _sq_nvsu(names, values, shapes, units):
    """Make a structured quantity from names, values, shapes, and units"""
    def nfsz(name, shape):
        if isinstance(shape, numbers.Number):
            size = shape
        elif len(shape) == 2:
            return (name, '<f8', (shape[1],))
        else:
            return (name, 'f8')
    dtype = [nfsz(nm,sh) for (nm,sh) in zip(names, shapes)]
    try:
        arr = np.array(list(zip(*values)), dtype=dtype)
    except:
        arr = np.array(tuple(values), dtype=dtype)
    sq = u.Quantity(arr, u.StructuredUnit(units))
    return sq

def sq_from_dict(d):
    """Make a structured quantity from a dictionary of quantities (u.Quantity)"""
    nvsu = (d.keys(), \
           [val.value for val in d.values()], \
           [val.shape for val in d.values()], \
           tuple([val.unit for val in d.values()]))
    return _sq_nvsu(*nvsu)

# apd = {'altper': demoa.propn.altperapo['altper'], 'altapo': demoa.propn.altperapo['altapo'], 'position': demoa.propn.pvt.cartesian['position']}
# apq = sq_from_dict(apd)


####################################################################
##### Build structured quantities
####################################################################

def _make_structured_quantity(q, name, scalar=False):
    """
    Convert a non-structured Quantity into a structured Quantity with one field.

    Parameters
    ----------
    q : astropy.units.Quantity
        Non-structured Quantity to convert (scalar or array)
    name : str
        Field name for the structured array
    scalar : bool, optional
        If True, treat the entire input as a single scalar element in the
        structured array (resulting in isscalar=True for the structured quantity).
        If False (default), each element becomes a separate row in the
        structured array (resulting in isscalar=False).

    Returns
    -------
    structured_qty : astropy.units.Quantity
        Structured quantity with StructuredUnit containing one field

    Examples
    --------
    # Scalar input
    >>> make_structured_quantity(2*u.m, 'distance')
    # Creates single-element structured array

    # Vector input, scalar=False (default)
    >>> make_structured_quantity([1, 2, 3]*u.m, 'distance', scalar=False)
    # Creates 3-element structured array, each element is a scalar
    # isscalar = False

    # Vector input, scalar=True
    >>> make_structured_quantity([1, 2, 3]*u.m, 'my3vec', scalar=True)
    # Creates 1-element structured array, the element is a 3-vector
    # isscalar = True
    """
    if scalar or q.isscalar:
        # Treat entire input as a single structured element
        # This handles both scalar inputs and vectors that should be kept as single elements
        dtype = [(name, q.dtype, q.shape)]
        struct_array = np.empty(1, dtype=dtype)
        struct_array[name][0] = q.value

        struct_unit = u.StructuredUnit((q.unit,), names=(name,))
        result = u.Quantity(struct_array, unit=struct_unit)

        # Return the scalar element (not the 1-element array)
        return result[0]
    else:
        # Each element of q becomes a separate row in the structured array
        q_array = np.atleast_1d(q)
        dtype = [(name, q.dtype)]
        struct_array = np.empty(len(q_array), dtype=dtype)
        struct_array[name] = q_array.value

        struct_unit = u.StructuredUnit((q.unit,), names=(name,))
        return u.Quantity(struct_array, unit=struct_unit)

u.Quantity.structure = lambda self, name: _make_structured_quantity(self, name)

def hstack(sqs):
    """Concatenate the quantities with different structure components
    and the same number of rows; sqs is a tuple or list of structured
    quantities.

    """
    isscalars = [sq.isscalar for sq in sqs]
    if all(isscalars):
        compat = True
    elif not(any(isscalars)):
        lens = [len(sq) for sq in sqs]
        compat = all(x==units[0] for x in units)
    else:
        compat = False
    if compat:
        return rfn.merge_arrays(sqs, flatten=True)
    else:
        raise ValueError("All lengths must be the same")

def vstack(sqs):
    """Concatenate rows the structured quantities with identical structure."""
    units = [sq.unit for sq in sqs]
    if all(x==units[0] for x in units):
        return rfn.stack_arrays([sq.value for sq in sqs])*units[0]
    else:
        raise ValueError("All units must be the same")

############################################################
## Make scalar sq from lists and dicts: structquant()
############################################################

def scvec(size):
    if size==1:
        return f"f8"
    else:
        return f"({size},)f8"

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
    (rvals, runits, rnames, vect) = _namedquant_tuple(values, names, units, unitlookup)
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
      Each value can be scalars or vectors
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
    Tuple of (values, units, names, vect)

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
            vals0 = vals[0]
        else:
            vals = [v*u.dimensionless_unscaled for v in vals]
            vals0 = vals
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
