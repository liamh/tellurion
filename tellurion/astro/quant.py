"""Make and convert quantities"""

import numpy as np
from numpy.lib import recfunctions as rfn
import astropy.units as u
import astropy.coordinates as coord

####################################################################
##### Build quantities
####################################################################

def make_quantity(value, unit=None, isscalar = False, unitlookup={}):
    """Make an unstructured quantity given a value, a unit (which may
    be a physical type), and optionally a dictionary `unitlookup`
    which maps physical types to the names units.

    value: number, array-like (argument to np.array), quantity, or dict
    unit: unit-like (argument to u.Unit), or a dict with the same keys as `value`

    unitlookup: dict, keys are strings that name a physical type,
    and every physical type used in `unit` should be present in the
    dict

    """
    def qv(value, sngun):
        if type(value) in (u.Quantity, coord.Angle, coord.Longitude, coord.Latitude, coord.Distance):
            q = u.Quantity(value)
        elif type(sngun) is str:
            unt = unitlookup.get(sngun) or sngun
            if unt:
                q = value*u.Unit(unt)
            else:
                q = value*u.Unit('1')
        else:
            raise TypeError("Could not assign a unit to the value")
        return q

    if type(value) is dict:
        if type(unit) is dict:
            quants = [_make_structured_quantity(qv(value[key], unit[key]), key, isscalar) for key in value]
        else:
            quants = [_make_structured_quantity(qv(value[key], None), key, isscalar) for key in value]
        return hstack(tuple(quants))
    else:
        return qv(value, unit)

def changeunits(qsq, unitlookup={}):
    """Change the units for the quantity or structured quantity to the
    system of units. The physical types of all units that occur in qsq
    must have a matching key in the `unitlookup` dictionary."""

    def getanypt(un):
        ptl = u.get_physical_type(un)._physical_type
        return next((unitlookup.get(k) for k in ptl if unitlookup.get(k) is not None), None)

    if type(qsq.unit) is u.StructuredUnit:
        pt = [getanypt(un) for un in qsq.unit.values()]
        if all(pt):
            tounits = u.StructuredUnit(tuple(pt))
        else:
            tounits = None
    else:
        tounits = getanypt(qsq.unit)

    if tounits:
        return qsq.to(tounits)
    else:
        raise ValueError("Unit physical type not found in `unitlookup`")

####################################################################
##### Build structured quantities
####################################################################

# Make the scalar structured quantity a singleton vector
u.Quantity.tovector = lambda self: u.Quantity([self]) if self.isscalar else self

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

    # 2D array input
    >>> make_structured_quantity(np.ones((13, 3))*u.m, 'coords')
    # Creates 13-element structured array, each element is a 3-vector
    # isscalar = False
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

        # Determine the shape of each element
        # For shape (13, 3), we want 13 rows, each with shape (3,)
        if q_array.ndim > 1:
            element_shape = q_array.shape[1:]
            n_elements = q_array.shape[0]
            dtype = [(name, q.dtype, element_shape)]
        else:
            n_elements = len(q_array)
            dtype = [(name, q.dtype)]

        struct_array = np.empty(n_elements, dtype=dtype)
        struct_array[name] = q_array.value

        struct_unit = u.StructuredUnit((q.unit,), names=(name,))
        return u.Quantity(struct_array, unit=struct_unit)

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
        compat = all(x==lens[0] for x in lens)
    else:
        compat = False
    if compat:
        ma = rfn.merge_arrays(sqs, flatten=True)
        if all(isscalars):
            return ma[0]
        else:
            return ma
    else:
        raise ValueError("All lengths must be the same")

def vstack(sqs):
    """Concatenate rows the structured quantities with identical structure."""
    units = [sq.unit for sq in sqs]
    if all(x==units[0] for x in units):
        return rfn.stack_arrays([sq.tovector().value for sq in sqs])*units[0]
    else:
        raise ValueError("All units must be the same")

####################################################################
##### Convert structured and unstructured quantities
####################################################################

# Make a dictionary by structure components, used to see what the SQ
# contents because it's not clear in the default print form
u.Quantity.to_dict = lambda self: {nm: self[nm] for nm in self.dtype.names} \
    if isinstance(self.unit, u.StructuredUnit) else {'': self.value}

# Extract value as np.array from quantities
u.Quantity.to_array = lambda self: rfn.structured_to_unstructured(self.value) \
    if isinstance(self.unit, u.StructuredUnit) else self.value
