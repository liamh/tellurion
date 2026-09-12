"""
Comprehensive test suite for posvel module.

Tests position and velocity handling with AstroPy-style conventions.
"""
import importlib

import pytest
import numpy as np
import astropy.units as u
import astropy.time
from astropy.timeseries import TimeSeries
from numpy.testing import assert_allclose, assert_array_equal
from types import SimpleNamespace

# Import the module to test
from tellurion.core import posvel
from tellurion.astro import units


# Test data: satellite states (km, km/s, MJD)
# Dump data in SI at full precision
#  np.set_printoptions(precision=16)
#  demoa.propn.ephem.pvt().to_array()
SATELLITE_STATES_m = np.array([[ 5740132.68349, 3314067.15, 0., -2750.82684, 4764.5718400000005, 5501.65367, 60676.],
       [4581815.8086414365, 4512263.175300685, 1616826.633579345, -4891.448994373261, 3141.5916390539796, 5166.422663426453, 60676.00347222222],
       [2865499.627198302, 5161045.3831239175, 3036846.5973876943, -6432.386721110699, 1140.4152997483332, 4203.821979784255, 60676.006944444445],
       [ 801265.5186635376, 5183536.72810986, 4088441.727473798, -7187.958175073694, -990.1159832057442, 2736.513492596365, 60676.010416666664],
       [1359960.4135013185, 4580209.018283869, 4646557.569878746, -7074.061712797389, -2989.0709511248806, 948.4194787381766, 60676.01388888889],
       [3358332.5947880326, 3427320.1009210977, 4647312.570345052, -6115.173294019763, -4617.482818602055, -941.2707751617605, 60676.01736111111],
       [4956791.986645598, 1865868.1956239555, 4094285.249660138, -4436.283289105959, -5686.778784033234, -2706.75324733447, 60676.020833333336],
       [5968503.631257458, 83291.80532639672, 3056384.6340837656, -2243.2158897858735, -6078.384156990985, -4142.427147818573, 60676.024305555555],
       [6277207.599631215, -1709253.233206629, 1658347.0778792806, 204.18625084661036, -5753.624175544211, -5084.877823796807, 60676.02777777778],
       [5848976.975063754, -3301215.6206595753, 65551.89665418287, 2621.940233178246, -4754.696329165089, -5428.657923351877, 60676.03125],
       [4734934.460366363, -4506102.992803526, -1534932.4332061035, 4731.935799049163, -3198.1666192098874, -5135.661435806881, 60676.03472222222],
       [3065182.546311113, -5182038.024056856, -2955185.2981953565, 6290.263956710347, -1262.3432221895998, -4238.353275838781, 60676.038194444445],
       [1034645.953016234, -5247728.3845713185, -4027343.115521708, 7112.606684996648, 830.5532720344688, -2837.023108899671, 60676.041666666664]])

SATELLITE_STATES = np.hstack((SATELLITE_STATES_m[:,0:6]/1000, SATELLITE_STATES_m[:,6].reshape(13,1)))

class TestPositionT:
    """Tests for PositionT class (position without velocity)."""

    def test_create_from_cartesian(self):
        """Test creating PositionT from Cartesian coordinates."""
        pos = SATELLITE_STATES[0,0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        pt = posvel.PositionT(time=time, cartesian=pos)

        assert pt.time == time
        assert pt.isscalar is True
        assert_allclose(pt.cartesian.to(u.km).value, pos.value, rtol=1e-10)

    def test_create_without_time(self):
        """Test creating PositionT without time information."""
        pos = SATELLITE_STATES[0,0:3] * u.km

        pt = posvel.PositionT(time=None, cartesian=pos)

        assert pt.time is None
        assert_allclose(pt.cartesian.to(u.km).value, pos.value, rtol=1e-10)

    def test_position_vector_property(self):
        """Test position_vector property extraction."""
        pos = SATELLITE_STATES[0,0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        pt = posvel.PositionT(time=time, cartesian=pos)

        assert_allclose(pt.position_vector.to(u.km).value, pos.value, rtol=1e-10)

    def test_has_velocity_false(self):
        """Test that PositionT reports no velocity."""
        pos = SATELLITE_STATES[0,0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        pt = posvel.PositionT(time=time, cartesian=pos)

        assert pt.has_velocity is False

    def test_copy(self):
        """Test copying PositionT."""
        pos = SATELLITE_STATES[0,0:3] * u.km
        time = astropy.time.Time(60676.0, format='mjd')
        aux = {'label': 'test'}

        pt1 = posvel.PositionT(time=time, cartesian=pos, aux=aux)
        pt2 = pt1.copy()

        # Modify original
        pt1.aux['label'] = 'modified'

        # Copy should be independent
        assert pt2.aux['label'] == 'test'
        assert_allclose(pt2.cartesian.to(u.km).value, pos.value, rtol=1e-10)

    def test_multiple_positions(self):
        """Test creating PositionT with multiple positions."""
        positions = SATELLITE_STATES[:, :3] * u.km
        times = astropy.time.Time(SATELLITE_STATES[:, 6], format='mjd')

        pt = posvel.PositionT(time=times, cartesian=positions)

        assert pt.isscalar is False
        assert len(pt) == len(SATELLITE_STATES)
        assert_allclose(pt.cartesian.to(u.km).value, positions.value, rtol=1e-10)

    def test_indexing(self):
        """Test indexing into multiple positions."""
        positions = SATELLITE_STATES[:5, :3] * u.km
        times = astropy.time.Time(SATELLITE_STATES[:5, 6], format='mjd')

        pt = posvel.PositionT(time=times, cartesian=positions)

        # Test single index
        pt_single = pt[0]
        assert pt_single.isscalar is True
        assert_allclose(pt_single.cartesian.to(u.km).value,
                       positions[0].value, rtol=1e-10)

        # Test slice
        pt_slice = pt[1:3]
        assert len(pt_slice) == 2
        assert_allclose(pt_slice.cartesian.to(u.km).value,
                       positions[1:3].value, rtol=1e-10)

    def test_concatenate(self):
        """Test concatenating PositionT objects."""
        pos1 = SATELLITE_STATES[0, :3] * u.km
        pos2 = SATELLITE_STATES[1, :3] * u.km
        time1 = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')
        time2 = astropy.time.Time(SATELLITE_STATES[1, 6], format='mjd')

        pt1 = posvel.PositionT(time=time1, cartesian=pos1)
        pt2 = posvel.PositionT(time=time2, cartesian=pos2)

        pt_combined = pt1.concatenate(pt2)

        assert len(pt_combined) == 2
        assert pt_combined.isscalar is False

    def test_to_array_with_time(self):
        """Test converting to array with time."""
        pos = np.array([5740132.6835, 3314067.15, 0.0]) * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        pt = posvel.PositionT(time=time, cartesian=pos)
        arr = pt.to_array()

        assert arr.shape == (1, 4)  # [x, y, z, time]
        assert_allclose(arr[0, :3], pos.to(u.m).value, rtol=1e-10)
        assert_allclose(arr[0, 3], 60676.0, rtol=1e-10)

    def test_to_array_without_time(self):
        """Test converting to array without time."""
        pos = np.array([5740132.6835, 3314067.15, 0.0]) * u.km

        pt = posvel.PositionT(time=None, cartesian=pos)
        arr = pt.to_array()

        assert arr.shape == (3,)  # [x, y, z]
        assert_allclose(arr, pos.to(u.m).value, rtol=1e-10)

    def test_ephemeris(self):
        """Test creating ephemeris TimeSeries."""
        positions = SATELLITE_STATES[:3, :3] * u.km
        times = astropy.time.Time(SATELLITE_STATES[:3, 6], format='mjd')

        pt = posvel.PositionT(time=times, cartesian=positions)
        ts = pt.ephemeris()

        assert isinstance(ts, TimeSeries)
        assert len(ts) == 3


class TestPositionVelocityT:
    """Tests for PositionVelocityT class (position and velocity)."""

    def test_create_from_cartesian(self):
        """Test creating PositionVelocityT from Cartesian coordinates."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt = posvel.pvtcart(state, time)

        assert isinstance(pvt, posvel.PositionVelocityT)
        assert pvt.time == time
        assert pvt.isscalar is True
        assert pvt.has_velocity is True

    def test_position_extraction(self):
        """Test extracting position-only as PositionT."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt = posvel.pvtcart(state, time)
        pt = pvt.position

        assert isinstance(pt, posvel.PositionT)
        assert pt.has_velocity is False

    def test_pv_legacy_property(self):
        """Test legacy 'pv' property alias."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt = posvel.pvtcart(state, time)

        # pv should be alias for cartesian
        assert np.all(pvt.pv == pvt.cartesian)

    def test_copy(self):
        """Test copying PositionVelocityT."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt1 = posvel.pvtcart(state, time)
        pvt2 = pvt1.copy()

        # Verify independence
        assert pvt2.time == pvt1.time
        assert np.all(pvt2.cartesian == pvt1.cartesian)

    def test_multiple_states(self):
        """Test creating PositionVelocityT with multiple states."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)  # time in last column

        assert pvt.isscalar is False
        assert len(pvt) == len(SATELLITE_STATES)

    def test_indexing(self):
        """Test indexing into multiple states."""
        pvt = posvel.pvtcart(SATELLITE_STATES[:5], None)

        # Test single index
        pvt_single = pvt[0]
        assert pvt_single.isscalar is True

        # Test slice
        pvt_slice = pvt[1:3]
        assert len(pvt_slice) == 2

    def test_concatenate(self):
        """Test concatenating PositionVelocityT objects."""
        state1 = SATELLITE_STATES[0]
        state2 = SATELLITE_STATES[1]

        pvt1 = posvel.pvtcart(state1, None)
        pvt2 = posvel.pvtcart(state2, None)

        pvt_combined = pvt1.concatenate(pvt2)

        assert len(pvt_combined) == 2
        assert pvt_combined.isscalar is False

    def test_merge_and_timeorder(self):
        """Test merging and time ordering."""
        # Create states out of time order
        state1 = SATELLITE_STATES[2]
        state2 = SATELLITE_STATES[0]
        state3 = SATELLITE_STATES[1]

        pvt1 = posvel.pvtcart(state1, None)
        pvt2 = posvel.pvtcart(state2, None)
        pvt3 = posvel.pvtcart(state3, None)

        # Merge and order
        pvt_merged = pvt1.merge(pvt2).merge(pvt3)

        assert len(pvt_merged) == 3
        # Verify time ordering
        times = pvt_merged.time.mjd
        assert np.all(times[:-1] <= times[1:])

    def test_to_array(self):
        """Test converting to array."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)

        arr = pvt.to_array()

        assert arr.shape == (1, 7)  # [x, y, z, vx, vy, vz, time]
        # Position in meters
        assert_allclose(arr[0, :3], state[:3] * 1000, rtol=1e-10)
        # Velocity in m/s
        assert_allclose(arr[0, 3:6], state[3:6] * 1000, rtol=1e-10)
        # Time
        assert_allclose(arr[0, 6], state[6], rtol=1e-10)

    def test_to_array_multiple(self):
        """Test converting multiple states to array."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)
        arr = pvt.to_array()

        assert arr.shape == (len(SATELLITE_STATES), 7)
        # Verify first state
        assert_allclose(arr[0, :3], SATELLITE_STATES[0, :3] * 1000, rtol=1e-10)
        assert_allclose(arr[0, 3:6], SATELLITE_STATES[0, 3:6] * 1000, rtol=1e-10)

    def test_ephemeris(self):
        """Test creating ephemeris TimeSeries."""
        pvt = posvel.pvtcart(SATELLITE_STATES[:5], None)
        ts = pvt.ephemeris()

        assert isinstance(ts, TimeSeries)
        assert len(ts) == 5
        # Check that position and velocity columns exist
        assert 'position' in ts.colnames or any('position' in col for col in ts.colnames)

    def test_ephemeris_with_elapsed(self):
        """Test ephemeris with elapsed time column."""
        pvt = posvel.pvtcart(SATELLITE_STATES[:5], None)
        ts = pvt.ephemeris(elapsed=True)

        assert 'elapsed' in ts.colnames

    def test_pvt_method(self):
        """Test pvt() method returns self."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)

        assert pvt.pvt() is pvt

    def test_kepler_method_delegates_to_ork_converter(self, monkeypatch):
        """Test kepler() stays as the core wrapper and delegates lazily."""
        pvt = posvel.pvtcart(SATELLITE_STATES[0], None)
        sentinel = object()
        calls = []

        def fake_import(name):
            calls.append(("import", name))
            return SimpleNamespace(
                _kepler=lambda arg, **kwargs: calls.append(("kepler", arg, kwargs))
                or sentinel
            )

        monkeypatch.setattr(importlib, "import_module", fake_import)

        assert pvt.kepler(mean_time_element=False) is sentinel
        assert calls == [
            ("import", "tellurion.ork.element"),
            ("kepler", pvt, {"mean_time_element": False}),
        ]

    def test_equinoctial_method_delegates_to_ork_converter(self, monkeypatch):
        """Test equinoctial() stays as the core wrapper and delegates lazily."""
        pvt = posvel.pvtcart(SATELLITE_STATES[0], None)
        sentinel = object()
        calls = []

        def fake_import(name):
            calls.append(("import", name))
            return SimpleNamespace(
                _equinoctial=lambda arg, **kwargs: calls.append(
                    ("equinoctial", arg, kwargs)
                )
                or sentinel
            )

        monkeypatch.setattr(importlib, "import_module", fake_import)

        assert pvt.equinoctial(mean_time_element=False) is sentinel
        assert calls == [
            ("import", "tellurion.ork.element"),
            ("equinoctial", pvt, {"mean_time_element": False}),
        ]

    def test_circular_method_delegates_to_ork_converter(self, monkeypatch):
        """Test circular() stays as the core wrapper and delegates lazily."""
        pvt = posvel.pvtcart(SATELLITE_STATES[0], None)
        sentinel = object()
        calls = []

        def fake_import(name):
            calls.append(("import", name))
            return SimpleNamespace(
                _circular=lambda arg, **kwargs: calls.append(("circular", arg, kwargs))
                or sentinel
            )

        monkeypatch.setattr(importlib, "import_module", fake_import)

        assert pvt.circular(mean_time_element=False) is sentinel
        assert calls == [
            ("import", "tellurion.ork.element"),
            ("circular", pvt, {"mean_time_element": False}),
        ]


class TestPvtcartFunction:
    """Tests for the pvtcart() helper function."""

    def test_single_state_with_time(self):
        """Test creating from single state with explicit time."""
        state = SATELLITE_STATES[0, :6]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt = posvel.pvtcart(state, time)

        assert isinstance(pvt, posvel.PositionVelocityT)
        assert pvt.isscalar is True

    def test_single_state_time_in_array(self):
        """Test creating from single state with time in array."""
        state = SATELLITE_STATES[0]  # includes time as last element

        pvt = posvel.pvtcart(state, None)

        assert isinstance(pvt, posvel.PositionVelocityT)
        assert_allclose(pvt.time.mjd, SATELLITE_STATES[0, 6], rtol=1e-10)

    def test_multiple_states(self):
        """Test creating from multiple states."""
        pvt = posvel.pvtcart(SATELLITE_STATES, None)

        assert isinstance(pvt, posvel.PositionVelocityT)
        assert len(pvt) == len(SATELLITE_STATES)

    def test_position_only_3_elements(self):
        """Test creating PositionT from 3-element array."""
        pos = SATELLITE_STATES[0, :3]
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pt = posvel.pvtcart(pos, time)

        assert isinstance(pt, posvel.PositionT)
        assert pt.has_velocity is False

    def test_position_only_with_time_column(self):
        """Test creating PositionT from 4-column array."""
        pos_time = np.column_stack([SATELLITE_STATES[:3, :3],
                                     SATELLITE_STATES[:3, 6]])

        pt = posvel.pvtcart(pos_time, None)

        assert isinstance(pt, posvel.PositionT)
        assert len(pt) == 3

    def test_list_input(self):
        """Test creating from list of arrays."""
        states_list = [SATELLITE_STATES[i, :6] for i in range(3)]
        times = astropy.time.Time(SATELLITE_STATES[:3, 6], format='mjd')

        pvt = posvel.pvtcart(states_list, times)

        assert len(pvt) == 3


class TestCoordinateConversion:
    """Tests for coordinate system conversion."""

    def test_cartesian_to_spherical_position_only(self):
        """Test conversion from Cartesian to spherical coordinates for PositionT."""
        pos = SATELLITE_STATES[0, :3] * u.km
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pt = posvel.PositionT(time=time, cartesian=pos)
        sph = pt.spherical

        # Check structure
        assert 'rtasc' in sph.dtype.names
        assert 'decl' in sph.dtype.names
        assert 'distance' in sph.dtype.names

    def test_cartesian_to_spherical_with_velocity(self):
        """Test conversion from Cartesian to spherical coordinates for PositionVelocityT."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)

        sph = pvt.spherical

        # Check structure - should have position and velocity components
        assert 'rtasc' in sph.dtype.names
        assert 'decl' in sph.dtype.names
        assert 'distance' in sph.dtype.names
        assert 'rtasc_r' in sph.dtype.names
        assert 'decl_r' in sph.dtype.names
        assert 'distance_r' in sph.dtype.names


class TestAuxAttributes:
    """Tests for auxiliary attributes."""

    def test_create_with_aux(self):
        """Test creating with auxiliary attributes."""
        state = SATELLITE_STATES[0]
        aux = {'label': 'test_satellite', 'flag': 'active'}
        time = astropy.time.Time(SATELLITE_STATES[0, 6], format='mjd')

        pvt = posvel.pvtcart(state[:6], time)
        pvt.aux = aux

        assert pvt.aux['label'] == 'test_satellite'
        assert pvt.aux['flag'] == 'active'

    def test_aux_in_ephemeris(self):
        """Test that aux attributes appear in ephemeris."""
        states = SATELLITE_STATES[:3]
        pvt = posvel.pvtcart(states, None)
        pvt.aux = {'label': 'sat1 sat2 sat3'}

        ts = pvt.ephemeris()

        assert 'label' in ts.colnames

    def test_aux_concatenation(self):
        """Test auxiliary attribute handling during concatenation."""
        state1 = SATELLITE_STATES[0]
        state2 = SATELLITE_STATES[1]

        pvt1 = posvel.pvtcart(state1, None)
        pvt1.aux = {'label': 'sat1'}

        pvt2 = posvel.pvtcart(state2, None)
        pvt2.aux = {'label': 'sat2'}

        pvt_combined = pvt1.concatenate(pvt2)

        assert 'label' in pvt_combined.aux
        assert 'sat1' in pvt_combined.aux['label']
        assert 'sat2' in pvt_combined.aux['label']


class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_scalar_no_len(self):
        """Test that scalar objects raise TypeError for len()."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)

        with pytest.raises(TypeError, match="has no len"):
            len(pvt)

    def test_scalar_not_subscriptable(self):
        """Test that scalar objects are not subscriptable."""
        state = SATELLITE_STATES[0]
        pvt = posvel.pvtcart(state, None)

        with pytest.raises(TypeError, match="not subscriptable"):
            _ = pvt[0]

    def test_cannot_specify_both_coordinates(self):
        """Test that specifying both cartesian and spherical raises error."""
        pos = np.array([1.0, 2.0, 3.0]) * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        # Create a dummy spherical coordinate
        sph = pos  # This would be real spherical coords in practice

        with pytest.raises(ValueError, match="Cannot specify both"):
            posvel.PositionT(time=time, cartesian=pos, spherical=sph)

    def test_must_specify_coordinates(self):
        """Test that omitting both coordinates raises error."""
        time = astropy.time.Time(60676.0, format='mjd')

        with pytest.raises(ValueError, match="Must specify either"):
            posvel.PositionT(time=time)


class TestLegacyCompatibility:
    """Tests for legacy compatibility features."""

    def test_pvt_alias(self):
        """Test that PVT is an alias for PositionVelocityT."""
        assert posvel.PVT is posvel.PositionVelocityT

    def test_pv_parameter_legacy(self):
        """Test legacy pv parameter in initialization."""
        pos = np.array([5740132.6835, 3314067.15, 0.0]) * u.km
        time = astropy.time.Time(60676.0, format='mjd')

        # Using legacy pv parameter
        pt = posvel.PositionT(time=time, pv=pos)

        assert_allclose(pt.cartesian.to(u.km).value, pos.value, rtol=1e-10)


class TestIntegration:
    """Integration tests with complete workflows."""

    def test_full_satellite_ephemeris_workflow(self):
        """Test complete workflow from data to ephemeris."""
        # Create from test data
        pvt = posvel.pvtcart(SATELLITE_STATES, None)

        # Verify creation
        assert len(pvt) == len(SATELLITE_STATES)

        # Extract subset
        pvt_subset = pvt[2:8]
        assert len(pvt_subset) == 6

        # Create ephemeris
        eph = pvt_subset.ephemeris()
        assert isinstance(eph, TimeSeries)
        assert len(eph) == 6

        # Convert back to array
        arr = pvt_subset.to_array()
        assert arr.shape == (6, 7)

    def test_merge_different_orbits(self):
        """Test merging data from different time periods."""
        # Split data into two parts
        pvt1 = posvel.pvtcart(SATELLITE_STATES[:5], None)
        pvt2 = posvel.pvtcart(SATELLITE_STATES[5:], None)

        # Merge
        pvt_merged = pvt1.merge(pvt2)

        # Should have all data in time order
        assert len(pvt_merged) == len(SATELLITE_STATES)

        # Verify time ordering
        times = pvt_merged.time.mjd
        assert np.all(times[:-1] <= times[1:])

    def test_roundtrip_to_array_and_back(self):
        """Test converting to array and back."""
        pvt1 = posvel.pvtcart(SATELLITE_STATES[:3], None)

        # Convert to array (SI units)
        arr = pvt1.to_array()

        # Create new from array using SI units
        from tellurion.astro import units as tell_units
        pvt2 = posvel.pvtcart(arr, None, specunits=tell_units.siunits)

        # Should match
        assert len(pvt2) == len(pvt1)
        assert_allclose(pvt2.position_vector.to(u.m).value,
                        pvt1.position_vector.to(u.m).value,
                        rtol=1e-10)

class TestMonkeyPatchedMethods:
    """Tests for monkey-patched methods on AstroPy classes."""

    def test_timeseries_pvt_method(self):
        """Test .pvt() method added to TimeSeries."""
        pvt = posvel.pvtcart(SATELLITE_STATES[:3], None)
        ts = pvt.ephemeris()

        # Convert back using monkey-patched method
        pvt_back = ts.pvt()

        assert isinstance(pvt_back, (posvel.PositionT, posvel.PositionVelocityT))
        assert len(pvt_back) == 3

    def test_table_row_pvt_method(self):
        """Test .pvt() method added to Table Row."""
        pvt = posvel.pvtcart(SATELLITE_STATES[:3], None)
        ts = pvt.ephemeris()

        # Get single row and convert
        row = ts[0]
        pvt_single = row.pvt()

        assert isinstance(pvt_single, (posvel.PositionT, posvel.PositionVelocityT))
        assert pvt_single.isscalar is True


class TestQuantityExtension:
    """Tests for Quantity.tovector() extension."""

    def test_scalar_to_vector(self):
        """Test tovector() on scalar Quantity."""
        scalar_q = 5.0 * u.km

        vector_q = scalar_q.tovector()

        assert vector_q.shape == (1,)
        assert vector_q[0] == scalar_q

    def test_vector_unchanged(self):
        """Test tovector() on already-vector Quantity."""
        vector_q = np.array([1.0, 2.0, 3.0]) * u.km

        result = vector_q.tovector()

        assert np.all(result == vector_q)
        assert result.shape == vector_q.shape


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
