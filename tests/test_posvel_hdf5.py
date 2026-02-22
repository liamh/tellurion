"""
Comprehensive test suite for posvel_hdf5 module.

Tests HDF5 serialization and deserialization of PositionT and PositionVelocityT objects.
"""
import pytest
import numpy as np
import astropy.units as u
import astropy.time
from numpy.testing import assert_allclose, assert_array_equal
import tempfile
import os

# Import the modules to test
from tellurion.core import posvel
from tellurion.core import posvel_hdf5  # This registers the HDF5 serializers
from tellurion.astro import units

# Import HDF5 save/load functions
from fsc.hdf5_io import save, load

# Import astropy-hdf5io to register AstroPy serializers
import astropy_hdf5io


# Test data: satellite states (km, km/s, MJD)
# Same data as in test_posvel.py for consistency
SATELLITE_STATES_m = np.array([
    [5740132.68349, 3314067.15, 0., -2750.82684, 4764.5718400000005, 5501.65367, 60676.],
    [4581815.8086414365, 4512263.175300685, 1616826.633579345, -4891.448994373261, 3141.5916390539796, 5166.422663426453, 60676.00347222222],
    [2865499.627198302, 5161045.3831239175, 3036846.5973876943, -6432.386721110699, 1140.4152997483332, 4203.821979784255, 60676.006944444445],
    [801265.5186635376, 5183536.72810986, 4088441.727473798, -7187.958175073694, -990.1159832057442, 2736.513492596365, 60676.010416666664],
    [1359960.4135013185, 4580209.018283869, 4646557.569878746, -7074.061712797389, -2989.0709511248806, 948.4194787381766, 60676.01388888889],
    [3358332.5947880326, 3427320.1009210977, 4647312.570345052, -6115.173294019763, -4617.482818602055, -941.2707751617605, 60676.01736111111],
    [4956791.986645598, 1865868.1956239555, 4094285.249660138, -4436.283289105959, -5686.778784033234, -2706.75324733447, 60676.020833333336],
    [5968503.631257458, 83291.80532639672, 3056384.6340837656, -2243.2158897858735, -6078.384156990985, -4142.427147818573, 60676.024305555555],
    [6277207.599631215, -1709253.233206629, 1658347.0778792806, 204.18625084661036, -5753.624175544211, -5084.877823796807, 60676.02777777778],
    [5848976.975063754, -3301215.6206595753, 65551.89665418287, 2621.940233178246, -4754.696329165089, -5428.657923351877, 60676.03125],
    [4734934.460366363, -4506102.992803526, -1534932.4332061035, 4731.935799049163, -3198.1666192098874, -5135.661435806881, 60676.03472222222],
    [3065182.546311113, -5182038.024056856, -2955185.2981953565, 6290.263956710347, -1262.3432221895998, -4238.353275838781, 60676.038194444445],
    [1034645.953016234, -5247728.3845713185, -4027343.115521708, 7112.606684996648, 830.5532720344688, -2837.023108899671, 60676.041666666664]
])

SATELLITE_STATES = np.hstack((SATELLITE_STATES_m[:, 0:6] / 1000, SATELLITE_STATES_m[:, 6].reshape(13, 1)))


class TestPositionTHDF5:
    """Tests for PositionT HDF5 serialization."""

    def test_save_load_scalar_with_time(self):
        """Test saving and loading a scalar PositionT with time."""
        pos = SATELLITE_STATES[0, 0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')
        
        pt = posvel.PositionT(time=time, cartesian=pos)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pt, filename)
            
            # Load
            loaded_pt = load(filename)
            
            # Verify
            assert isinstance(loaded_pt, posvel.PositionT)
            assert loaded_pt.isscalar is True
            assert loaded_pt.time == time
            assert_allclose(loaded_pt.cartesian.to(u.km).value, pos.value, rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_scalar_without_time(self):
        """Test saving and loading a scalar PositionT without time."""
        pos = SATELLITE_STATES[0, 0:3] * u.km
        
        pt = posvel.PositionT(time=None, cartesian=pos)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pt, filename)
            
            # Load
            loaded_pt = load(filename)
            
            # Verify
            assert isinstance(loaded_pt, posvel.PositionT)
            assert loaded_pt.time is None
            assert_allclose(loaded_pt.cartesian.to(u.km).value, pos.value, rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_array(self):
        """Test saving and loading an array of PositionT."""
        positions = SATELLITE_STATES[:5, :3] * u.km
        times = astropy.time.Time(SATELLITE_STATES[:5, 6], format='mjd')
        
        pt = posvel.PositionT(time=times, cartesian=positions)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pt, filename)
            
            # Load
            loaded_pt = load(filename)
            
            # Verify
            assert isinstance(loaded_pt, posvel.PositionT)
            assert loaded_pt.isscalar is False
            assert len(loaded_pt) == 5
            assert_allclose(loaded_pt.time.mjd, times.mjd, rtol=1e-10)
            assert_allclose(loaded_pt.cartesian.to(u.km).value, positions.value, rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_with_aux(self):
        """Test saving and loading PositionT with auxiliary attributes."""
        pos = SATELLITE_STATES[0, 0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')
        aux = {'label': 'test_satellite', 'flag': 'active', 'mission_id': '12345'}
        
        pt = posvel.PositionT(time=time, cartesian=pos, aux=aux)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pt, filename)
            
            # Load
            loaded_pt = load(filename)
            
            # Verify
            assert isinstance(loaded_pt, posvel.PositionT)
            # Note: Simple string values will be preserved
            assert 'label' in loaded_pt.aux or 'flag' in loaded_pt.aux
        finally:
            os.unlink(filename)

    def test_save_load_preserves_units(self):
        """Test that units are preserved through save/load cycle."""
        # Use meters instead of km
        pos = SATELLITE_STATES[0, 0:3] * 1000 * u.m
        time = astropy.time.Time(60676.0, format='mjd')
        
        pt = posvel.PositionT(time=time, cartesian=pos)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pt, filename)
            
            # Load
            loaded_pt = load(filename)
            
            # Verify units are preserved (can convert back to original)
            assert_allclose(loaded_pt.cartesian.to(u.m).value, pos.value, rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_roundtrip_multiple_times(self):
        """Test that multiple save/load cycles don't degrade data."""
        positions = SATELLITE_STATES[:3, :3] * u.km
        times = astropy.time.Time(SATELLITE_STATES[:3, 6], format='mjd')
        
        pt = posvel.PositionT(time=times, cartesian=positions)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Multiple save/load cycles
            current_pt = pt
            for i in range(3):
                save(current_pt, filename)
                current_pt = load(filename)
            
            # Verify final result matches original
            assert_allclose(current_pt.cartesian.to(u.km).value, 
                          positions.value, rtol=1e-10)
            assert_allclose(current_pt.time.mjd, times.mjd, rtol=1e-10)
        finally:
            os.unlink(filename)


class TestPositionVelocityTHDF5:
    """Tests for PositionVelocityT HDF5 serialization."""

    def test_save_load_scalar_with_time(self):
        """Test saving and loading a scalar PositionVelocityT with time."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        
        pvt = posvel.pvtcart(state, time)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            assert loaded_pvt.isscalar is True
            assert loaded_pvt.time == time
            assert loaded_pvt.has_velocity is True
            
            # Verify position
            assert_allclose(loaded_pvt.position_vector.to(u.km).value,
                          state[:3], rtol=1e-10)
            
            # Verify velocity
            vel_loaded = loaded_pvt.cartesian['velocity'].to(u.km / u.s).value
            assert_allclose(vel_loaded, state[3:6], rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_with_time_in_array(self):
        """Test saving and loading with time in the array."""
        state = SATELLITE_STATES[0]  # includes time as last element
        
        pvt = posvel.pvtcart(state, None)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            assert_allclose(loaded_pvt.time.mjd, SATELLITE_STATES[0, 6], rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_array(self):
        """Test saving and loading an array of PositionVelocityT."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            assert loaded_pvt.isscalar is False
            assert len(loaded_pvt) == len(SATELLITE_STATES)
            
            # Verify times
            assert_allclose(loaded_pvt.time.mjd, SATELLITE_STATES[:, 6], rtol=1e-10)
            
            # Verify positions
            assert_allclose(loaded_pvt.position_vector.to(u.km).value,
                          SATELLITE_STATES[:, :3], rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_with_aux(self):
        """Test saving and loading PositionVelocityT with auxiliary attributes."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        aux = {'satellite': 'ISS', 'orbit_number': '12345'}
        
        pvt = posvel.pvtcart(state, time)
        pvt.aux = aux
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            # Aux attributes should be preserved (at least as strings)
            assert 'satellite' in loaded_pvt.aux or 'orbit_number' in loaded_pvt.aux
        finally:
            os.unlink(filename)

    def test_save_load_preserves_structured_quantity(self):
        """Test that structured quantity format is preserved."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify structured format
            assert hasattr(loaded_pvt.cartesian, 'dtype')
            assert loaded_pvt.cartesian.dtype.names is not None
            assert 'position' in loaded_pvt.cartesian.dtype.names
            assert 'velocity' in loaded_pvt.cartesian.dtype.names
        finally:
            os.unlink(filename)

    def test_roundtrip_accuracy(self):
        """Test that save/load preserves numerical accuracy."""
        states = SATELLITE_STATES[:5]
        pvt = posvel.pvtcart(states, None)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Convert both to arrays and compare
            original_array = pvt.to_array()
            loaded_array = loaded_pvt.to_array()
            
            assert_allclose(loaded_array, original_array, rtol=1e-12)
        finally:
            os.unlink(filename)

    def test_save_load_subset(self):
        """Test saving and loading a subset of states."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)
        
        # Extract subset
        pvt_subset = pvt[2:7]
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save subset
            save(pvt_subset, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert len(loaded_pvt) == 5
            assert_allclose(loaded_pvt.time.mjd, SATELLITE_STATES[2:7, 6], rtol=1e-10)
        finally:
            os.unlink(filename)


class TestNestedStructures:
    """Tests for saving PositionT/PositionVelocityT in nested structures."""

    def test_save_load_in_dict(self):
        """Test saving PositionVelocityT nested in a dictionary."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        
        data = {
            'orbit': pvt,
            'epoch': time,
            'mission_id': 'TEST-001'
        }
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(data, filename)
            
            # Load
            loaded_data = load(filename)
            
            # Verify
            assert isinstance(loaded_data, dict)
            assert 'orbit' in loaded_data
            assert isinstance(loaded_data['orbit'], posvel.PositionVelocityT)
            assert_allclose(loaded_data['orbit'].time.mjd, time.mjd, rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_in_list(self):
        """Test saving multiple PositionVelocityT objects in a list."""
        pvt1 = posvel.pvtcart(SATELLITE_STATES[0], None)
        pvt2 = posvel.pvtcart(SATELLITE_STATES[1], None)
        pvt3 = posvel.pvtcart(SATELLITE_STATES[2], None)
        
        data = [pvt1, pvt2, pvt3]
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(data, filename)
            
            # Load
            loaded_data = load(filename)
            
            # Verify
            assert isinstance(loaded_data, list)
            assert len(loaded_data) == 3
            for i, pvt in enumerate(loaded_data):
                assert isinstance(pvt, posvel.PositionVelocityT)
                assert_allclose(pvt.time.mjd, SATELLITE_STATES[i, 6], rtol=1e-10)
        finally:
            os.unlink(filename)

    def test_save_load_mixed_with_quantities(self):
        """Test saving PositionVelocityT mixed with other AstroPy types."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)
        
        data = {
            'state': pvt,
            'altitude': 408 * u.km,
            'period': 92.5 * u.min,
            'inclination': 51.6 * u.deg
        }
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(data, filename)
            
            # Load
            loaded_data = load(filename)
            
            # Verify
            assert isinstance(loaded_data['state'], posvel.PositionVelocityT)
            assert_allclose(loaded_data['altitude'].to(u.km).value, 408, rtol=1e-10)
            assert_allclose(loaded_data['period'].to(u.min).value, 92.5, rtol=1e-10)
        finally:
            os.unlink(filename)


class TestEdgeCases:
    """Tests for edge cases and special situations."""

    def test_empty_aux_dict(self):
        """Test saving and loading with empty aux dictionary."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        
        pvt = posvel.pvtcart(state, time)
        # Explicitly set empty aux
        pvt.aux = {}
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            assert loaded_pvt.aux == {}
        finally:
            os.unlink(filename)

    def test_single_element_array(self):
        """Test saving a 1-element array (not scalar)."""
        states = SATELLITE_STATES[:1]  # Single row but as array
        pvt = posvel.pvtcart(states, None)
        
        # Force it to be non-scalar
        pvt.isscalar = False
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Should still be length 1
            assert len(loaded_pvt) == 1
        finally:
            os.unlink(filename)

    def test_very_large_array(self):
        """Test saving and loading a large array of states."""
        # Create a larger dataset by repeating
        large_states = np.tile(SATELLITE_STATES, (10, 1))
        
        pvt = posvel.pvtcart(large_states, None)
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(pvt, filename)
            
            # Load
            loaded_pvt = load(filename)
            
            # Verify
            assert len(loaded_pvt) == len(large_states)
            assert_allclose(loaded_pvt.position_vector.to(u.km).value,
                          large_states[:, :3], rtol=1e-10)
        finally:
            os.unlink(filename)


class TestCompatibility:
    """Tests for compatibility and integration."""

    def test_position_to_position_velocity_compatibility(self):
        """Test that PositionT and PositionVelocityT can be saved/loaded independently."""
        pos = SATELLITE_STATES[0, :3] * u.km
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        
        pt = posvel.PositionT(time=time, cartesian=pos)
        pvt = posvel.pvtcart(state, time)
        
        with tempfile.NamedTemporaryFile(suffix='_pt.hdf5', delete=False) as f:
            pt_filename = f.name
        with tempfile.NamedTemporaryFile(suffix='_pvt.hdf5', delete=False) as f:
            pvt_filename = f.name
        
        try:
            # Save both
            save(pt, pt_filename)
            save(pvt, pvt_filename)
            
            # Load both
            loaded_pt = load(pt_filename)
            loaded_pvt = load(pvt_filename)
            
            # Verify types are preserved
            assert isinstance(loaded_pt, posvel.PositionT)
            assert isinstance(loaded_pvt, posvel.PositionVelocityT)
            assert not loaded_pt.has_velocity
            assert loaded_pvt.has_velocity
        finally:
            os.unlink(pt_filename)
            os.unlink(pvt_filename)

    def test_interoperability_with_native_astropy(self):
        """Test that native AstroPy types can be saved alongside posvel objects."""
        from astropy.coordinates import SkyCoord
        
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)
        
        # Create some native AstroPy objects
        coord = SkyCoord(ra=10 * u.degree, dec=40 * u.degree, distance=1000 * u.pc)
        time = astropy.time.Time('2023-01-01T00:00:00')
        
        data = {
            'satellite_state': pvt,
            'target_coordinate': coord,
            'observation_time': time
        }
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save
            save(data, filename)
            
            # Load
            loaded_data = load(filename)
            
            # Verify all types preserved
            assert isinstance(loaded_data['satellite_state'], posvel.PositionVelocityT)
            assert isinstance(loaded_data['target_coordinate'], SkyCoord)
            assert isinstance(loaded_data['observation_time'], astropy.time.Time)
        finally:
            os.unlink(filename)


class TestErrorHandling:
    """Tests for error conditions and validation."""

    def test_file_not_found(self):
        """Test that loading from non-existent file raises appropriate error."""
        with pytest.raises((FileNotFoundError, OSError)):
            load('nonexistent_file.hdf5')

    def test_corrupted_file_handling(self):
        """Test handling of corrupted HDF5 file."""
        # Create a text file with .hdf5 extension
        with tempfile.NamedTemporaryFile(mode='w', suffix='.hdf5', delete=False) as f:
            filename = f.name
            f.write("This is not an HDF5 file")
        
        try:
            with pytest.raises(Exception):  # Should raise some kind of error
                load(filename)
        finally:
            os.unlink(filename)


class TestIntegrationWorkflows:
    """Integration tests for complete workflows."""

    def test_full_orbit_propagation_workflow(self):
        """Test a complete workflow: create, save, load, modify, save again."""
        # Create initial orbit
        pvt = posvel.pvtcart(SATELLITE_STATES[:5], None)
        
        with tempfile.NamedTemporaryFile(suffix='_initial.hdf5', delete=False) as f:
            initial_file = f.name
        with tempfile.NamedTemporaryFile(suffix='_modified.hdf5', delete=False) as f:
            modified_file = f.name
        
        try:
            # Save initial
            save(pvt, initial_file)
            
            # Load and modify
            loaded_pvt = load(initial_file)
            loaded_pvt.aux['mission'] = 'ISS'
            
            # Save modified
            save(loaded_pvt, modified_file)
            
            # Load modified
            final_pvt = load(modified_file)
            
            # Verify
            assert len(final_pvt) == 5
            assert 'mission' in final_pvt.aux
        finally:
            os.unlink(initial_file)
            os.unlink(modified_file)

    def test_merge_from_files(self):
        """Test loading multiple files and merging."""
        pvt1 = posvel.pvtcart(SATELLITE_STATES[:5], None)
        pvt2 = posvel.pvtcart(SATELLITE_STATES[5:10], None)
        
        with tempfile.NamedTemporaryFile(suffix='_1.hdf5', delete=False) as f:
            file1 = f.name
        with tempfile.NamedTemporaryFile(suffix='_2.hdf5', delete=False) as f:
            file2 = f.name
        
        try:
            # Save both
            save(pvt1, file1)
            save(pvt2, file2)
            
            # Load both
            loaded1 = load(file1)
            loaded2 = load(file2)
            
            # Merge
            merged = loaded1.merge(loaded2)
            
            # Verify
            assert len(merged) == 10
        finally:
            os.unlink(file1)
            os.unlink(file2)

    def test_extract_subset_and_save(self):
        """Test extracting a subset and saving it separately."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)
        
        # Extract subset
        subset = pvt[3:8]
        
        with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
            filename = f.name
        
        try:
            # Save subset
            save(subset, filename)
            
            # Load
            loaded_subset = load(filename)
            
            # Verify it's the correct subset
            assert len(loaded_subset) == 5
            assert_allclose(loaded_subset.time.mjd, SATELLITE_STATES[3:8, 6], rtol=1e-10)
        finally:
            os.unlink(filename)


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])