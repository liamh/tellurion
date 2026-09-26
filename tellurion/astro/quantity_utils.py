"""Utilities for creating and manipulating structured Quantity objects.

This module provides utilities for working with both regular and structured
astropy Quantity objects, including functions to create structured quantities,
convert between unit systems, and stack quantities.

"""

from typing import Any, Dict, Optional, Sequence, Union

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from numpy.lib import recfunctions as rfn

__all__ = [
    'make_quantity',
    'change_units',
    'hstack',
    'vstack',
    'quantity_to_dict',
    'quantity_to_array',
]


def sifloat(q):
    return float(q.si.value)

####################################################################
# Build quantities
####################################################################

def make_quantity(
    value: Union[float, np.ndarray, u.Quantity, Dict[str, Any]],
    unit: Optional[Union[u.Unit, str, Dict[str, Union[u.Unit, str]]]] = None,
    is_scalar: bool = False,
    unit_lookup: Optional[Dict[str, Union[u.Unit, str]]] = None
) -> u.Quantity:
    """Create a Quantity from a value and unit specification.

    This function creates either a regular or structured Quantity depending
    on the input types. It supports physical type lookups through a unit
    dictionary, allowing you to specify units by their physical type name
    (e.g., "length", "velocity") rather than specific units.

    Parameters
    ----------
    value : number, array-like, Quantity, or dict
        The numerical value(s) to create a Quantity from. If a dict, each
        key becomes a field name in a structured Quantity.
    unit : Unit, str, or dict, optional
        The unit(s) for the value(s). Can be:
        - A Unit object or string parseable by `u.Unit()`
        - A physical type name (must be in `unit_lookup`)
        - A dict with same keys as `value` (if `value` is a dict)
        If None, assumes dimensionless (unit='1').
    is_scalar : bool, optional
        Only applicable when `value` is a dict. If True, treat the entire
        input as a single element in the structured quantity (resulting in
        isscalar=True). If False (default), each element becomes a separate
        row. Raises ValueError if True when value is not a dict.
        Default is False.
    unit_lookup : dict, optional
        Mapping from physical type names (str) to unit specifications.
        Used when `unit` contains physical type names instead of actual units.
        Keys should be physical type strings, values should be Unit-like objects.

    Returns
    -------
    quantity : Quantity
        A regular Quantity if `value` is not a dict, or a structured Quantity
        with StructuredUnit if `value` is a dict.

    Raises
    ------
    TypeError
        If unable to assign units to the value.
    ValueError
        If `value` and `unit` are both dicts with mismatched keys.

    Examples
    --------
    Create a simple quantity:

    >>> import tellurion as tell
    >>> tell.make_quantity(5.0, 'meter')
    <Quantity 5. m>

    Create a quantity using physical type lookup:

    >>> unit_lookup = {'length': 'm', 'time': 's'}
    >>> tell.make_quantity(10.0, 'length', unit_lookup=unit_lookup)
    <Quantity 10. m>

    Create a structured quantity from a dict:

    >>> values = {'x': 1.0, 'y': 2.0}
    >>> units = {'x': 'm', 'y': 'm'}
    >>> tell.make_quantity(values, units)
    <Quantity (1., 2.) m>

    Create a structured quantity with different units:

    >>> values = {'position': 10.0, 'velocity': 5.0}
    >>> units = {'position': 'm', 'velocity': 'm/s'}
    >>> tell.make_quantity(values, units)
    <Quantity (10., 5.) (m, m / s)>

    Wrap an array as a single vector element:

    >>> values = {'position': [1, 2, 3]}
    >>> tell.make_quantity(values, 'm', is_scalar=True)
    <Quantity ([1, 2, 3],) m>

    Wrap an existing Quantity into a structured format:

    >>> existing_q = 5.0 * u.m
    >>> tell.make_quantity({'distance': existing_q.value}, existing_q.unit)
    <Quantity (5.,) m>

    Notes
    -----
    When `value` is a dict and `unit` is also a dict, they must have matching keys.
    The resulting structured Quantity will have fields corresponding to the keys.

    """
    if unit_lookup is None:
        unit_lookup = {}

    def _quantity_from_value(
            val: Any, unit_spec: Optional[Union[u.Unit, str]]
    ) -> u.Quantity:
        """Helper to create a Quantity from a value and unit specification."""
        # If already a Quantity-like object, convert to Quantity
        if isinstance(val, (u.Quantity, coord.Angle, coord.Longitude,
                           coord.Latitude, coord.Distance)):
            return u.Quantity(val)

        # If unit_spec is a string, try to look it up or parse it
        if isinstance(unit_spec, str):
            # Try to get from lookup first (for physical types)
            resolved_unit = unit_lookup.get(unit_spec, unit_spec)
            if resolved_unit:
                return val * u.Unit(resolved_unit)
            else:
                # Dimensionless
                return val * u.dimensionless_unscaled

        # If unit_spec is already a Unit
        if isinstance(unit_spec, u.UnitBase):
            return val * unit_spec

        # If unit_spec is None
        if unit_spec is None:
            # Check if val has units already
            if hasattr(val, 'unit'):
                return u.Quantity(val)
            # Otherwise dimensionless
            return val * u.dimensionless_unscaled

        raise TypeError(
            f"Cannot assign unit to value. unit_spec type {type(unit_spec)} "
            f"not recognized."
        )

    # Handle dict input - create structured quantity
    if isinstance(value, dict):
        if isinstance(unit, dict):
            # Both value and unit are dicts
            # Check that all value keys have corresponding units
            value_keys = set(value.keys())
            unit_keys = set(unit.keys())

            if not value_keys.issubset(unit_keys):
                missing_keys = value_keys - unit_keys
                raise ValueError(
                    f"All keys in value dict must have corresponding units. "
                    f"Missing units for keys: {missing_keys}"
                )

            # Only use the units that correspond to values
            quants = [
                _make_structured_quantity(_quantity_from_value(value[key], unit[key]), key, is_scalar)
                for key in value
            ]
        else:
            # value is dict but unit is not - apply same unit to all
            quants = [
                _make_structured_quantity(_quantity_from_value(value[key], unit), key, is_scalar)
                for key in value
            ]
        return hstack(quants)
    elif is_scalar:
        raise ValueError(
            "is_scalar=True is only valid when value is a dict. "
            "For regular quantities, scalar/array behavior is determined by the value itself."
        )

    # Handle non-dict input - create regular quantity
    return _quantity_from_value(value, unit)

def change_units(
    quantity: u.Quantity,
    unit_lookup: Dict[str, Union[u.Unit, str]]
) -> u.Quantity:
    """Convert a Quantity to units specified by physical type lookup.

    This function changes the units of a regular or structured Quantity
    according to a mapping from physical types to desired units. For each
    unit in the Quantity, it determines the physical type and looks up the
    corresponding target unit in `unit_lookup`.

    Parameters
    ----------
    quantity : Quantity
        The Quantity (regular or structured) to convert.
    unit_lookup : dict
        Mapping from physical type names to target units. Keys should be
        physical type strings (e.g., 'length', 'time'), values should be
        Unit-like objects or strings parseable by `u.Unit()`.

    Returns
    -------
    converted : Quantity
        The Quantity converted to the units specified in `unit_lookup`.

    Raises
    ------
    ValueError
        If a physical type in the Quantity's units is not found in `unit_lookup`.

    Examples
    --------
    Convert a simple quantity:

    >>> import tellurion as tell
    >>> q = 1000.0 * u.meter
    >>> unit_lookup = {'length': 'km'}
    >>> tell.change_units(q, unit_lookup)
    <Quantity 1. km>

    Convert a structured quantity:

    >>> # Create structured quantity with cm and hours
    >>> q = tell.make_quantity({'height': 100, 'duration': 2}, {'height': 'cm', 'duration': 'hour'})
    >>> unit_lookup = {'length': 'm', 'time': 's'}
    >>> tell.change_units(q, unit_lookup)
    <Quantity (1., 7200.) (m, s)>

    Notes
    -----
    This function uses `astropy.units.get_physical_type()` to determine the
    physical type of each unit. All physical types present in the Quantity
    must have corresponding entries in `unit_lookup`.

    """
    def _get_target_unit(unit: u.Unit) -> Optional[u.Unit]:
        """Get the target unit for a given unit based on its physical type."""
        # Get all physical types for this unit
        physical_types = u.get_physical_type(unit)._physical_type

        # Try to find a match in unit_lookup
        for phys_type in physical_types:
            target = unit_lookup.get(phys_type)
            if target is not None:
                return u.Unit(target) if isinstance(target, str) else target

        return None

    # Handle structured units
    if isinstance(quantity.unit, u.StructuredUnit):
        target_units = [_get_target_unit(unit) for unit in quantity.unit.values()]

        if all(tu is not None for tu in target_units):
            target_structured_unit = u.StructuredUnit(tuple(target_units))
            return quantity.to(target_structured_unit)
        else:
            # Find which physical type is missing
            missing_types = []
            for unit in quantity.unit.values():
                if _get_target_unit(unit) is None:
                    phys_types = u.get_physical_type(unit)._physical_type
                    missing_types.extend(phys_types)
            raise ValueError(
                f"Physical types {set(missing_types)} not found in unit_lookup. "
                f"Available types: {set(unit_lookup.keys())}"
            )

    # Handle regular units
    target_unit = _get_target_unit(quantity.unit)
    if target_unit is not None:
        return quantity.to(target_unit)
    else:
        phys_types = u.get_physical_type(quantity.unit)._physical_type
        raise ValueError(
            f"Physical types {phys_types} not found in unit_lookup. "
            f"Available types: {set(unit_lookup.keys())}"
        )


####################################################################
# Build structured quantities (internal)
####################################################################

def _make_structured_quantity(
    quantity: u.Quantity,
    name: str,
    is_scalar: bool = False
) -> u.Quantity:
    """Convert a regular Quantity into a structured Quantity with one field.

    Internal function used by make_quantity(). Not part of the public API.

    Parameters
    ----------
    quantity : Quantity
        Regular (non-structured) Quantity to convert. Can be scalar or array.
    name : str
        Field name for the structured array.
    is_scalar : bool, optional
        If True, treat the entire input as a single scalar element in the
        structured array (resulting in a scalar structured Quantity).
        If False (default), each element becomes a separate row in the
        structured array (resulting in an array structured Quantity).
        Default is False.

    Returns
    -------
    structured_quantity : Quantity
        Structured Quantity with StructuredUnit containing one field.
        - If `is_scalar=True` or input is scalar: returns a scalar structured Quantity
        - If `is_scalar=False` and input is array: returns an array structured Quantity

    Notes
    -----
    This is an internal implementation detail. Users should use make_quantity()
    with a dict for value to create structured quantities.

    """
    # Handle scalar inputs or explicit scalar treatment
    if is_scalar or quantity.isscalar:
        # Create a 1-element structured array
        dtype = [(name, quantity.dtype, quantity.shape)]
        struct_array = np.empty(1, dtype=dtype)
        struct_array[name][0] = quantity.value

        struct_unit = u.StructuredUnit((quantity.unit,), names=(name,))
        result = u.Quantity(struct_array, unit=struct_unit)

        # Return the scalar element (not the 1-element array)
        return result[0]

    # Handle array inputs - each element becomes a separate row
    q_array = np.atleast_1d(quantity)

    # Determine the shape of each element
    # For shape (N, M, ...), we want N rows, each with shape (M, ...)
    if q_array.ndim > 1:
        element_shape = q_array.shape[1:]
        n_elements = q_array.shape[0]
        dtype = [(name, quantity.dtype, element_shape)]
    else:
        n_elements = len(q_array)
        dtype = [(name, quantity.dtype)]

    # Create the structured array
    struct_array = np.empty(n_elements, dtype=dtype)
    struct_array[name] = q_array.value

    struct_unit = u.StructuredUnit((quantity.unit,), names=(name,))
    return u.Quantity(struct_array, unit=struct_unit)


####################################################################
# Stack structured quantities
####################################################################

def hstack(quantities: Sequence[u.Quantity]) -> u.Quantity:
    """Horizontally stack structured quantities (concatenate fields).

    Combines multiple structured Quantities by adding their fields together.
    All input quantities must have the same number of rows (or all be scalar).
    This is analogous to adding columns to a table.

    Parameters
    ----------
    quantities : sequence of Quantity
        Structured Quantities to concatenate. Must all have the same length
        (if arrays) or all be scalar.

    Returns
    -------
    stacked : Quantity
        Structured Quantity containing all fields from input quantities.
        If all inputs are scalar, returns a scalar. Otherwise returns an array.

    Raises
    ------
    ValueError
        If quantities have inconsistent lengths (mixing scalars with arrays
        of different lengths).

    Examples
    --------
    Stack two scalar structured quantities:

    >>> import tellurion as tell
    >>> q1 = tell.make_quantity({'x': 5}, 'm')
    >>> q2 = tell.make_quantity({'y': 10}, 'm')
    >>> tell.hstack([q1, q2])
    <Quantity (5., 10.) m>

    Stack two array structured quantities:

    >>> q1 = tell.make_quantity({'x': [1, 2, 3]}, 'm')
    >>> q2 = tell.make_quantity({'y': [4, 5, 6]}, 'm')
    >>> tell.hstack([q1, q2])
    <Quantity [(1., 4.), (2., 5.), (3., 6.)] m>

    Stack quantities with different units:

    >>> q1 = tell.make_quantity({'position': [1, 2]}, 'm')
    >>> q2 = tell.make_quantity({'velocity': [3, 4]}, 'm/s')
    >>> tell.hstack([q1, q2])
    <Quantity [(1., 3.), (2., 4.)] (m, m / s)>

    See Also
    --------
    vstack : Vertically stack structured quantities (concatenate rows)
    astropy.table.hstack : Similar function for Table objects

    Notes
    -----
    This function uses numpy's `merge_arrays` under the hood to combine
    the structured arrays while preserving all field information and units.

    """
    # Check if all are scalar or all have the same length
    is_scalar_list = [q.isscalar for q in quantities]

    if all(is_scalar_list):
        # All scalar - straightforward merge
        compatible = True
    elif not any(is_scalar_list):
        # All arrays - check lengths match
        lengths = [len(q) for q in quantities]
        compatible = all(length == lengths[0] for length in lengths)
    else:
        # Mix of scalar and array - not compatible
        compatible = False

    if not compatible:
        lengths_info = [
            f"{i}: {'scalar' if q.isscalar else f'length {len(q)}'}"
            for i, q in enumerate(quantities)
        ]
        raise ValueError(
            f"All structured quantities must have the same length or all be scalar. "
            f"Got: {', '.join(lengths_info)}"
        )

    # Merge the arrays
    merged_array = rfn.merge_arrays(quantities, flatten=True)

    # Return scalar element if all inputs were scalar
    if all(is_scalar_list):
        return merged_array[0]
    else:
        return merged_array


def vstack(quantities: Sequence[u.Quantity]) -> u.Quantity:
    """Vertically stack structured quantities (concatenate rows).

    Combines multiple structured Quantities by stacking their rows.
    All input quantities must have identical structure (same field names
    and units). This is analogous to adding rows to a table.

    Parameters
    ----------
    quantities : sequence of Quantity
        Structured Quantities to stack. Must all have the same StructuredUnit
        (same field names and units for each field).

    Returns
    -------
    stacked : Quantity
        Structured Quantity containing all rows from input quantities,
        stacked vertically.

    Raises
    ------
    ValueError
        If quantities have different StructuredUnits (different fields or units).

    Examples
    --------
    Stack two structured quantities with the same structure:

    >>> import tellurion as tell
    >>> q1 = tell.make_quantity({'x': [1, 2]}, 'm')
    >>> q2 = tell.make_quantity({'x': [3, 4]}, 'm')
    >>> tell.vstack([q1, q2])
    <Quantity [(1.,), (2.,), (3.,), (4.,)] m>

    Stack scalar quantities:

    >>> q1 = tell.make_quantity({'x': 5}, 'm')
    >>> q2 = tell.make_quantity({'x': 10}, 'm')
    >>> tell.vstack([q1, q2])
    <Quantity [(5.,), (10.,)] m>

    Stack multi-field structured quantities:

    >>> q1 = tell.make_quantity({'x': [1, 2], 'y': [3, 4]}, 'm')
    >>> q2 = tell.make_quantity({'x': [5, 6], 'y': [7, 8]}, 'm')
    >>> tell.vstack([q1, q2])
    <Quantity [(1., 3.), (2., 4.), (5., 7.), (6., 8.)] (m, m)>

    See Also
    --------
    hstack : Horizontally stack structured quantities (concatenate fields)
    astropy.table.vstack : Similar function for Table objects

    Notes
    -----
    This function uses numpy's `stack_arrays` to combine the rows while
    maintaining the structured array format. Scalar inputs are automatically
    converted to 1-element arrays before stacking.

    """
    # Check that all units are the same
    units = [q.unit for q in quantities]
    if not all(unit == units[0] for unit in units):
        unit_strs = [str(unit) for unit in units]
        raise ValueError(
            f"All structured quantities must have the same unit structure. "
            f"Got units: {', '.join(unit_strs)}"
        )

    # Convert scalars to 1-element vectors for stacking
    quantities_as_vectors = [
        u.Quantity([q]) if q.isscalar else q
        for q in quantities
    ]

    # Stack the arrays
    stacked_values = rfn.stack_arrays(
        [q.value for q in quantities_as_vectors],
        usemask=False
    )

    return u.Quantity(stacked_values, unit=units[0])


####################################################################
# Convert structured and unstructured quantities
####################################################################

def quantity_to_dict(quantity: u.Quantity) -> Dict[str, Any]:
    """Convert a structured Quantity to a dictionary.

    Extracts the fields of a structured Quantity into a dictionary where
    keys are field names and values are the corresponding Quantities.
    For regular (non-structured) Quantities, returns a dict with empty
    string key.

    .. note::

        Importing :mod:`tellurion` monkey-patches this functionality directly onto
        :class:`astropy.units.Quantity` as :meth:`~astropy.units.Quantity.to_dict`.

    Parameters
    ----------
    quantity : Quantity
        Quantity to convert (structured or regular).

    Returns
    -------
    dict
        Dictionary representation of the Quantity.
        - For structured Quantities: {field_name: sub_quantity, ...}
        - For regular Quantities: {'': array_value}

    Examples
    --------
    Convert a structured quantity:

    >>> import tellurion as tell
    >>> q = tell.make_quantity({'x': 5, 'y': 10}, 'm')
    >>> tell.quantity_to_dict(q)
    {'x': <Quantity 5. m>, 'y': <Quantity 10. m>}
    >>> q.to_dict()
    {'x': <Quantity 5. m>, 'y': <Quantity 10. m>}

    Convert a regular quantity:

    >>> q = 5 * u.m
    >>> tell.quantity_to_dict(q)
    {'': 5.0}

    """
    if isinstance(quantity.unit, u.StructuredUnit):
        return {name: quantity[name] for name in quantity.dtype.names}
    else:
        return {'': quantity.value}

u.Quantity.to_dict = quantity_to_dict

def quantity_to_array(quantity: u.Quantity) -> np.ndarray:
    """Extract numerical values from a Quantity as a numpy array.

    Converts a Quantity to a pure numpy array, discarding unit information.
    For structured Quantities, converts to an unstructured array using
    numpy's structured_to_unstructured function.

    .. note::

        Importing :mod:`tellurion` monkey-patches this functionality directly onto
        :class:`astropy.units.Quantity` as :meth:`~astropy.units.Quantity.to_array`.

    Parameters
    ----------
    quantity : Quantity
        Quantity to convert (structured or regular).

    Returns
    -------
    array : ndarray
        Numpy array containing the numerical values.
        - For regular Quantities: returns the .value attribute
        - For structured Quantities: returns flattened unstructured array

    Examples
    --------
    Convert a regular quantity:

    >>> import tellurion as tell
    >>> q = [1, 2, 3] * u.m
    >>> tell.quantity_to_array(q)
    array([1., 2., 3.])

    Convert a structured quantity:

    >>> q = tell.make_quantity({'x': [1, 2], 'y': [3, 4]}, 'm')
    >>> tell.quantity_to_array(q)
    array([[1., 3.],
           [2., 4.]])

    Notes
    -----
    This function is useful when you need to pass quantity values to
    functions that don't understand astropy units, or for numerical
    computations where units have already been verified to be consistent.

    """
    if isinstance(quantity.unit, u.StructuredUnit):
        return rfn.structured_to_unstructured(quantity.value)
    else:
        return quantity.value

u.Quantity.to_array = quantity_to_array
