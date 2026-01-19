import pytest
import numpy as np
import astropy.units as u
import tellurion as tell
from tellurion.core.nquant import _namedquant_tuple

@pytest.mark.parametrize(
    "input_dict,expected_output",
    [
        # Test case 0
        (
            {'values': {'alength': np.float64(1000.0),
                'amass': np.float64(200.0),
                'atime': np.float64(3.0)}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('1', '1', '1'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 1
        (
            {'values': {'alength': np.float64(1000.0),
                'amass': np.float64(200.0),
                'atime': np.float64(3.0)},
             'units': {'alength': 'km', 'amass': 'kg', 'atime': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 2
        (
            {'values': {'alength': np.float64(1000.0),
                'amass': np.float64(200.0),
                'atime': np.float64(3.0)},
             'units': {'alength': 'length', 'amass': 'mass', 'atime': 'time'},
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 3
        (
            {'values': [np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             'names': ['alength', 'amass', 'atime']},
            ([np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             ('1', '1', '1'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 4
        (
            {'values': [np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             'names': ['alength', 'amass', 'atime'],
             'units': ['km', 'kg', 'h']},
            ([np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 5
        (
            {'values': [np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             'names': ['alength', 'amass', 'atime'],
             'units': ['km', 'kg', 'h'],
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ([np.float64(1000.0), np.float64(200.0), np.float64(3.0)],
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 6 - Quantity objects
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 7 - Quantity objects with unit conversion
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'units': {'alength': 'm', 'amass': 'g', 'atime': 's'}},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('m', 'g', 's'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 8 - Quantity objects with unitlookup
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'unitlookup': {'length': 'meter', 'mass': 'gram', 'time': 'second'}},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('meter', 'gram', 'second'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 9 - Quantity objects with units and unitlookup
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'units': {'alength': 'length', 'amass': 'mass', 'atime': 'time'},
             'unitlookup': {'length': 'meter', 'mass': 'gram', 'time': 'second'}},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('meter', 'gram', 'second'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 10 - Quantity objects with explicit units matching internal
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'units': {'alength': 'km', 'amass': 'kg', 'atime': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 11 - Quantity objects with unitlookup matching internal
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 12 - Quantity objects with units and unitlookup matching
        (
            {'values': {'alength': 1000.*u.km,
                'amass': 200.*u.kg,
                'atime': 3.*u.h},
             'units': {'alength': 'length', 'amass': 'mass', 'atime': 'time'},
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ('alength', 'amass', 'atime'))
        ),
        # Test case 13 - List of Quantity objects with names
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime']},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 14 - List of Quantity objects with unit conversion
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'units': ['m', 'g', 's']},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('m', 'g', 's'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 15 - List of Quantity objects with unitlookup
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'unitlookup': {'length': 'meter', 'mass': 'gram', 'time': 'second'}},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('meter', 'gram', 'second'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 16 - List of Quantity objects with units and unitlookup
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'units': ['length', 'mass', 'time'],
             'unitlookup': {'length': 'meter', 'mass': 'gram', 'time': 'second'}},
            ((np.float64(1000000.0), np.float64(200000.0), np.float64(10800.0)),
             ('meter', 'gram', 'second'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 17 - List of Quantity objects with matching units
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'units': ['km', 'kg', 'h']},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 18 - List of Quantity objects with matching units and unitlookup
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'units': ['km', 'kg', 'h'],
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
        # Test case 19 - List of Quantity objects with units as lookups
        (
            {'values': [1000.*u.km, 200.*u.kg, 3.*u.h],
             'names': ['alength', 'amass', 'atime'],
             'units': ['length', 'mass', 'time'],
             'unitlookup': {'length': 'km', 'mass': 'kg', 'time': 'h'}},
            ((np.float64(1000.0), np.float64(200.0), np.float64(3.0)),
             ('km', 'kg', 'h'),
             ['alength', 'amass', 'atime'])
        ),
    ],
    ids=[f"test_case_{i}" for i in range(20)]
)
def test_namedquant_tuple(input_dict, expected_output):
    """Test _namedquant_tuple function with various input configurations."""
    result = _namedquant_tuple(**input_dict)

    # Unpack expected output
    expected_values, expected_units, expected_names = expected_output
    result_values, result_units, result_names, result_vect = result

    # Compare values
    assert not result_vect
    if isinstance(expected_values, tuple):
        assert isinstance(result_values, tuple), "Result values should be a tuple"
        assert len(result_values) == len(expected_values), "Value lengths should match"
        for rv, ev in zip(result_values, expected_values):
            assert np.isclose(rv, ev), f"Values should match: {rv} != {ev}"
    else:  # list
        assert isinstance(result_values, list), "Result values should be a list"
        assert len(result_values) == len(expected_values), "Value lengths should match"
        for rv, ev in zip(result_values, expected_values):
            assert np.isclose(rv, ev), f"Values should match: {rv} != {ev}"

    # Compare units
    assert result_units == expected_units, f"Units should match: {result_units} != {expected_units}"

    # Compare names
    assert result_names == expected_names, f"Names should match: {result_names} != {expected_names}"
