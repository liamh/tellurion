"""
Tests for quantity_utils module.

These tests cover:
- Creating regular and structured quantities
- Unit conversion with physical type lookup
- Stacking operations (hstack and vstack)
- Conversion utilities (to_dict, to_array)
- Edge cases and error conditions
"""

import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
import astropy.units as u
from astropy.coordinates import Angle, Distance

# Import the functions to test
from tellurion.astro.quantity_utils import (
    make_quantity,
    change_units,
    hstack,
    vstack,
    quantity_to_dict,
    quantity_to_array,
)


class TestMakeQuantity:
    """Tests for make_quantity function."""

    def test_simple_scalar(self):
        """Test creating a simple scalar quantity."""
        q = make_quantity(5.0, 'meter')
        assert isinstance(q, u.Quantity)
        assert q.value == 5.0
        assert q.unit == u.meter
        assert q.isscalar

    def test_simple_array(self):
        """Test creating a simple array quantity."""
        q = make_quantity([1, 2, 3], 'm')
        assert isinstance(q, u.Quantity)
        assert_array_equal(q.value, [1, 2, 3])
        assert q.unit == u.m
        assert not q.isscalar

    def test_from_existing_quantity(self):
        """Test creating quantity from existing Quantity."""
        q_input = 10 * u.km
        q = make_quantity(q_input, 'm')
        # Should preserve the input quantity
        assert q.unit == u.km  # Unit from input takes precedence

    def test_from_angle(self):
        """Test creating quantity from Angle."""
        angle = Angle(90, unit='deg')
        q = make_quantity(angle, 'deg')
        assert isinstance(q, u.Quantity)
        assert q.unit == u.deg

    def test_from_distance(self):
        """Test creating quantity from Distance."""
        dist = Distance(10, unit='pc')
        q = make_quantity(dist, 'pc')
        assert isinstance(q, u.Quantity)
        assert q.unit == u.pc

    def test_with_unit_lookup(self):
        """Test using unit_lookup for physical type resolution."""
        unit_lookup = {'length': 'meter', 'time': 'second'}
        q = make_quantity(5.0, 'length', unit_lookup=unit_lookup)
        assert q.unit == u.meter
        assert q.value == 5.0

    def test_dimensionless(self):
        """Test creating dimensionless quantity."""
        q = make_quantity(5.0, None)
        assert q.unit == u.dimensionless_unscaled

    def test_dict_value_single_unit(self):
        """Test creating structured quantity from dict with single unit."""
        values = {'x': 1.0, 'y': 2.0}
        q = make_quantity(values, 'm')
        assert isinstance(q.unit, u.StructuredUnit)
        assert q['x'].value == 1.0
        assert q['y'].value == 2.0
        assert q['x'].unit == u.m
        assert q['y'].unit == u.m

    def test_dict_value_dict_unit(self):
        """Test creating structured quantity from dicts."""
        values = {'position': 10.0, 'velocity': 5.0}
        units = {'position': 'm', 'velocity': 'm/s'}
        q = make_quantity(values, units)
        assert isinstance(q.unit, u.StructuredUnit)
        assert q['position'].value == 10.0
        assert q['velocity'].value == 5.0
        assert q['position'].unit == u.m
        assert q['velocity'].unit == u.m / u.s

    def test_dict_mismatched_keys(self):
        """Test error when dict keys don't match."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'x': 'm', 'z': 's'}  # 'z' instead of 'y'
        with pytest.raises(ValueError, match="Keys in value dict.*must match"):
            make_quantity(values, units)

    def test_unit_lookup_with_dict(self):
        """Test unit_lookup with dict values."""
        values = {'pos': 100.0, 'vel': 10.0}
        units = {'pos': 'length', 'vel': 'velocity'}
        unit_lookup = {'length': 'cm', 'velocity': 'km/s'}
        q = make_quantity(values, units, unit_lookup=unit_lookup)
        assert q['pos'].unit == u.cm
        assert q['vel'].unit == u.km / u.s

    def test_is_scalar_flag(self):
        """Test is_scalar flag for structured quantities."""
        values = {'vec': [1, 2, 3]}
        q = make_quantity(values, 'm', is_scalar=True)
        assert q.isscalar
        assert_array_equal(q['vec'].value, [1, 2, 3])

    def test_wrap_existing_quantity(self):
        """Test wrapping an existing Quantity into structured format."""
        existing_q = 5.0 * u.m
        q = make_quantity({'distance': existing_q.value}, existing_q.unit)
        assert isinstance(q.unit, u.StructuredUnit)
        assert q['distance'].value == 5.0
        assert q['distance'].unit == u.m

    def test_dict_mismatched_keys(self):
        """Test error when value keys are missing from unit dict."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'x': 'm', 'z': 's'}  # Missing 'y'
        with pytest.raises(ValueError, match="Missing units for keys"):
            make_quantity(values, units)

    def test_dict_exact_match(self):
        """Test when value and unit dicts have exact same keys."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'x': 'm', 'y': 's'}
        q = make_quantity(values, units)
        assert set(q.dtype.names) == {'x', 'y'}
        assert q['x'].unit == u.m
        assert q['y'].unit == u.s

    def test_dict_unit_superset(self):
        """Test when unit dict has more keys than value dict."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'x': 'm', 'y': 's', 'z': 'kg', 'vx': 'm/s'}  # Extra keys
        q = make_quantity(values, units)
        # Should only create fields for x and y
        assert set(q.dtype.names) == {'x', 'y'}
        assert q['x'].value == 1.0
        assert q['y'].value == 2.0
        assert q['x'].unit == u.m
        assert q['y'].unit == u.s

    def test_dict_value_superset_fails(self):
        """Test that having more value keys than unit keys fails."""
        values = {'x': 1.0, 'y': 2.0, 'z': 3.0}
        units = {'x': 'm', 'y': 's'}  # Missing 'z'
        with pytest.raises(ValueError, match="Missing units for keys.*z"):
            make_quantity(values, units)

    def test_dict_no_overlap_fails(self):
        """Test that non-overlapping keys fail."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'a': 'm', 'b': 's'}  # Completely different keys
        with pytest.raises(ValueError, match="Missing units for keys"):
            make_quantity(values, units)

    def test_dict_partial_overlap_fails(self):
        """Test that partial overlap with missing keys fails."""
        values = {'x': 1.0, 'y': 2.0}
        units = {'y': 's', 'z': 'kg'}  # Missing 'x', has extra 'z'
        with pytest.raises(ValueError, match="Missing units for keys.*x"):
            make_quantity(values, units)

class TestChangeUnits:
    """Tests for change_units function."""

    def test_simple_conversion(self):
        """Test simple unit conversion."""
        q = 1000.0 * u.m
        unit_lookup = {'length': 'km'}
        result = change_units(q, unit_lookup)
        assert result.value == 1.0
        assert result.unit == u.km

    def test_no_conversion_needed(self):
        """Test when already in correct units."""
        q = 5.0 * u.km
        unit_lookup = {'length': 'km'}
        result = change_units(q, unit_lookup)
        assert result.value == 5.0
        assert result.unit == u.km

    def test_structured_quantity(self):
        """Test converting structured quantity."""
        q = make_quantity({'height': 100, 'duration': 2}, {'height': 'cm', 'duration': 'hour'})

        unit_lookup = {'length': 'm', 'time': 's'}
        result = change_units(q, unit_lookup)

        assert result['height'].value == 1.0
        assert result['height'].unit == u.m
        assert result['duration'].value == 7200.0
        assert result['duration'].unit == u.s

    def test_missing_physical_type(self):
        """Test error when physical type not in lookup."""
        q = 5.0 * u.m
        unit_lookup = {'time': 's'}  # Missing 'length'
        with pytest.raises(ValueError, match="Physical types.*not found"):
            change_units(q, unit_lookup)

    def test_array_quantity(self):
        """Test converting array quantity."""
        q = [1, 2, 3] * u.m
        unit_lookup = {'length': 'cm'}
        result = change_units(q, unit_lookup)
        assert_array_equal(result.value, [100, 200, 300])
        assert result.unit == u.cm


class TestHstack:
    """Tests for hstack function."""

    def test_scalar_quantities(self):
        """Test stacking scalar structured quantities."""
        q1 = make_quantity({'x': 5}, 'm')
        q2 = make_quantity({'y': 10}, 'm')
        result = hstack([q1, q2])

        assert result.isscalar
        assert result['x'].value == 5
        assert result['y'].value == 10
        assert 'x' in result.dtype.names
        assert 'y' in result.dtype.names

    def test_array_quantities(self):
        """Test stacking array structured quantities."""
        q1 = make_quantity({'x': [1, 2, 3]}, 'm')
        q2 = make_quantity({'y': [4, 5, 6]}, 'm')
        result = hstack([q1, q2])

        assert not result.isscalar
        assert len(result) == 3
        assert_array_equal(result['x'].value, [1, 2, 3])
        assert_array_equal(result['y'].value, [4, 5, 6])

    def test_different_units(self):
        """Test stacking quantities with different units."""
        q1 = make_quantity({'position': [1, 2]}, 'm')
        q2 = make_quantity({'velocity': [3, 4]}, 'm/s')
        result = hstack([q1, q2])

        assert result['position'].unit == u.m
        assert result['velocity'].unit == u.m / u.s

    def test_mismatched_lengths(self):
        """Test error when lengths don't match."""
        q1 = make_quantity({'x': [1, 2]}, 'm')
        q2 = make_quantity({'y': [3, 4, 5]}, 'm')
        with pytest.raises(ValueError, match="same length"):
            hstack([q1, q2])

    def test_mixed_scalar_array(self):
        """Test error when mixing scalar and array."""
        q1 = make_quantity({'x': 5}, 'm')
        q2 = make_quantity({'y': [1, 2]}, 'm')
        with pytest.raises(ValueError, match="same length"):
            hstack([q1, q2])

    def test_multiple_quantities(self):
        """Test stacking more than two quantities."""
        q1 = make_quantity({'x': [1, 2]}, 'm')
        q2 = make_quantity({'y': [3, 4]}, 'm')
        q3 = make_quantity({'z': [5, 6]}, 'm')
        result = hstack([q1, q2, q3])

        assert len(result) == 2
        assert set(result.dtype.names) == {'x', 'y', 'z'}


class TestVstack:
    """Tests for vstack function."""

    def test_array_quantities(self):
        """Test stacking array quantities."""
        q1 = make_quantity({'x': [1, 2]}, 'm')
        q2 = make_quantity({'x': [3, 4]}, 'm')
        result = vstack([q1, q2])

        assert len(result) == 4
        assert_array_equal(result['x'].value, [1, 2, 3, 4])
        assert result['x'].unit == u.m

    def test_scalar_quantities(self):
        """Test stacking scalar quantities."""
        q1 = make_quantity({'x': 5}, 'm')
        q2 = make_quantity({'x': 10}, 'm')
        result = vstack([q1, q2])

        assert len(result) == 2
        assert_array_equal(result['x'].value, [5, 10])

    def test_multi_field(self):
        """Test stacking multi-field structured quantities."""
        # Create first structured quantity
        q1 = make_quantity({'x': [1, 2], 'y': [3, 4]}, 'm')

        # Create second structured quantity
        q2 = make_quantity({'x': [5, 6], 'y': [7, 8]}, 'm')

        # Stack them
        result = vstack([q1, q2])

        assert len(result) == 4
        assert_array_equal(result['x'].value, [1, 2, 5, 6])
        assert_array_equal(result['y'].value, [3, 4, 7, 8])

    def test_mismatched_units(self):
        """Test error when units don't match."""
        q1 = make_quantity({'x': [1, 2]}, 'm')
        q2 = make_quantity({'x': [3, 4]}, 'km')
        with pytest.raises(ValueError, match="same unit structure"):
            vstack([q1, q2])

    def test_mismatched_structure(self):
        """Test error when structures don't match."""
        q1 = make_quantity({'x': [1, 2]}, 'm')
        q2 = make_quantity({'x': [3, 4], 'y': [5, 6]}, 'm')

        with pytest.raises(ValueError, match="same unit structure"):
            vstack([q1, q2])


class TestQuantityToDict:
    """Tests for quantity_to_dict function."""

    def test_structured_quantity(self):
        """Test converting structured quantity to dict."""
        q = make_quantity({'x': 5, 'y': 10}, 'm')

        result = quantity_to_dict(q)

        assert isinstance(result, dict)
        assert set(result.keys()) == {'x', 'y'}
        assert result['x'].value == 5
        assert result['y'].value == 10

    def test_regular_quantity(self):
        """Test converting regular quantity to dict."""
        q = 5 * u.m
        result = quantity_to_dict(q)

        assert isinstance(result, dict)
        assert '' in result
        assert result[''] == 5.0

    def test_array_structured_quantity(self):
        """Test with array structured quantity."""
        q = make_quantity({'x': [1, 2, 3], 'y': [4, 5, 6]}, 'm')

        result = quantity_to_dict(q)

        assert set(result.keys()) == {'x', 'y'}
        assert_array_equal(result['x'].value, [1, 2, 3])
        assert_array_equal(result['y'].value, [4, 5, 6])


class TestQuantityToArray:
    """Tests for quantity_to_array function."""

    def test_regular_scalar(self):
        """Test converting regular scalar quantity."""
        q = 5 * u.m
        result = quantity_to_array(q)
        assert result == 5.0

    def test_regular_array(self):
        """Test converting regular array quantity."""
        q = [1, 2, 3] * u.m
        result = quantity_to_array(q)
        assert_array_equal(result, [1, 2, 3])

    def test_structured_quantity(self):
        """Test converting structured quantity."""
        q = make_quantity({'x': [1, 2], 'y': [3, 4]}, 'm')

        result = quantity_to_array(q)

        expected = np.array([[1, 3], [2, 4]])
        assert_array_equal(result, expected)

    def test_2d_structured(self):
        """Test with 2D data in structured quantity."""
        q = make_quantity({'coords': np.array([[1, 2], [3, 4], [5, 6]])}, 'm')

        result = quantity_to_array(q)

        # Should flatten the structure
        assert result.shape == (3, 2)


class TestEdgeCases:
    """Tests for edge cases and special scenarios."""

    def test_empty_array(self):
        """Test with empty array."""
        q = make_quantity({'empty': []}, 'm')
        assert len(q) == 0

    def test_zero_values(self):
        """Test with zero values."""
        q = make_quantity(0.0, 'm')
        assert q.value == 0.0
        assert q.unit == u.m

    def test_negative_values(self):
        """Test with negative values."""
        q = make_quantity(-5.0, 'm')
        assert q.value == -5.0

    def test_complex_units(self):
        """Test with complex compound units."""
        q = make_quantity(5.0, 'm/(s**2)')
        assert q.unit == u.m / u.s**2

    def test_very_large_numbers(self):
        """Test with very large numbers."""
        q = make_quantity(1e100, 'm')
        assert q.value == 1e100

    def test_very_small_numbers(self):
        """Test with very small numbers."""
        q = make_quantity(1e-100, 'm')
        assert q.value == 1e-100


class TestIntegration:
    """Integration tests combining multiple operations."""

    def test_full_workflow(self):
        """Test a complete workflow."""
        # Create structured quantity from dict
        values = {'x': 100.0, 'y': 200.0, 'z': 300.0}
        units = {'x': 'cm', 'y': 'cm', 'z': 'cm'}
        q = make_quantity(values, units)

        # Convert units
        unit_lookup = {'length': 'm'}
        q_converted = change_units(q, unit_lookup)

        # Check conversion
        assert q_converted['x'].value == 1.0
        assert q_converted['y'].value == 2.0
        assert q_converted['z'].value == 3.0

        # Convert to dict
        q_dict = quantity_to_dict(q_converted)
        assert set(q_dict.keys()) == {'x', 'y', 'z'}

        # Convert to array
        q_array = quantity_to_array(q_converted)
        assert_array_equal(q_array, [1.0, 2.0, 3.0])

    def test_stacking_workflow(self):
        """Test combining hstack and vstack operations."""
        # Create first row
        row1 = make_quantity({'x': 1, 'y': 2}, 'm')

        # Create second row
        row2 = make_quantity({'x': 3, 'y': 4}, 'm')

        # Stack rows
        table = vstack([row1, row2])

        assert len(table) == 2
        assert_array_equal(table['x'].value, [1, 3])
        assert_array_equal(table['y'].value, [2, 4])

    def test_multidimensional_data(self):
        """Test with multidimensional data (e.g., position vectors)."""
        # Create position data: 10 positions, each 3D
        positions = np.random.rand(10, 3)
        velocities = np.random.rand(10, 3)

        q_pos = make_quantity({'position': positions}, 'm')
        q_vel = make_quantity({'velocity': velocities}, 'm/s')

        # Combine into phase space
        phase_space = hstack([q_pos, q_vel])

        assert len(phase_space) == 10
        assert phase_space['position'].shape == (10, 3)
        assert phase_space['velocity'].shape == (10, 3)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
