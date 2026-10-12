import astropy.time
import astropy.units as u
import numpy as np
import pytest

from tellurion.core import posvel


@pytest.fixture
def times():
    return astropy.time.Time([60676.0, 60676.01], format="mjd")


@pytest.fixture
def position_velocity(times):
    states = np.array([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]])
    return posvel.pvtcart(states, times)


@pytest.fixture
def position_only(times):
    positions = np.array([[1, 2, 3], [4, 5, 6]]) * u.km
    return posvel.PositionT(time=times, cartesian=positions)


def test_cartesian_velocity_uses_custom_columnnames(position_velocity):
    ts = position_velocity.ephemeris(
        elapsed=False, reftime=None, columnnames=["place", "motion"]
    )

    assert ts.colnames == ["time", "place", "motion"]
    assert ts["place"].shape == (2, 3)
    assert ts["motion"].shape == (2, 3)


def test_position_only_custom_names_and_split_components(position_only):
    compact = position_only.ephemeris(
        elapsed=False, reftime=None, columnnames=["coords"]
    )
    split = position_only.ephemeris(
        elapsed=False,
        reftime=None,
        columnnames=["coords"],
        components=True,
    )

    assert compact.colnames == ["time", "coords"]
    assert compact["coords"].shape == (2, 3)
    assert split.colnames == ["time", "coords_x", "coords_y", "coords_z"]


def test_spherical_velocity_uses_custom_columnnames(position_velocity):
    ts = position_velocity.ephemeris(
        elapsed=False,
        reftime=None,
        coordinate_type="spherical",
        columnnames=["sky", "sky_rate"],
    )

    assert ts.colnames == [
        "time",
        "sky_rtasc",
        "sky_decl",
        "sky_distance",
        "sky_rate_rtasc_r",
        "sky_rate_decl_r",
        "sky_rate_distance_r",
    ]
    assert ts["sky_rtasc"].info.format == posvel.pvhelper._pos_format
    assert ts["sky_rate_rtasc_r"].info.format == posvel.pvhelper._vel_format


def test_spherical_position_only_uses_custom_name(position_only):
    ts = position_only.ephemeris(
        elapsed=False,
        reftime=None,
        coordinate_type="spherical",
        columnnames=["sky"],
    )

    assert ts.colnames == ["time", "sky_rtasc", "sky_decl", "sky_distance"]
    assert ts["sky_distance"].info.format == posvel.pvhelper._pos_format


def test_cartesian_velocity_can_split_components(position_velocity):
    ts = position_velocity.ephemeris(
        elapsed=False,
        reftime=None,
        columnnames=["place", "motion"],
        components=True,
    )

    assert ts.colnames == [
        "time",
        "place_x",
        "place_y",
        "place_z",
        "motion_x",
        "motion_y",
        "motion_z",
    ]
    assert ts["place_x"].info.format == posvel.pvhelper._pos_format
    assert ts["motion_x"].info.format == posvel.pvhelper._vel_format


@pytest.mark.parametrize(
    "obj,names,expected",
    [
        ("velocity", ["place"], "2 data-column names"),
        ("velocity", ["time", "place", "motion"], "2 data-column names"),
        ("position", ["place", "motion"], "1 data-column name"),
        ("position", ["time", "place"], "1 data-column name"),
    ],
)
def test_invalid_columnnames_length_raises_value_error(
    obj, names, expected, position_velocity, position_only
):
    instance = position_velocity if obj == "velocity" else position_only

    with pytest.raises(ValueError, match=expected):
        instance.ephemeris(columnnames=names)


def test_dynamic_name_prefixes_data_but_not_time(position_velocity):
    position_velocity.name = "Probe"
    ts = position_velocity.ephemeris(
        elapsed=False, reftime=None, columnnames=["place", "motion"]
    )

    assert ts.colnames == ["time", "Probe place", "Probe motion"]


def test_elapsed_and_reference_time_with_custom_data_names(position_velocity):
    ts = position_velocity.ephemeris(
        columnnames=["place", "motion"], elapsed=True, reftime="epoch"
    )

    assert ts.time[0] == position_velocity.time[0]
    assert "elapsed" in ts.colnames
    assert "from epoch" in ts.colnames


def test_default_column_names_and_elapsed_reference_time(position_velocity):
    ts = position_velocity.ephemeris(elapsed=True, reftime="epoch")

    assert "time" in ts.colnames
    assert "position" in ts.colnames
    assert "velocity" in ts.colnames
    assert "elapsed" in ts.colnames
    assert "from epoch" in ts.colnames
    assert ts["position"].shape == (2, 3)
    assert ts["position"].info.format == posvel.pvhelper._pos_format
    assert ts["velocity"].info.format == posvel.pvhelper._vel_format


def test_dynamic_name_prefixes_default_names_only(position_velocity):
    position_velocity.name = "Satellite"
    ts = position_velocity.ephemeris(elapsed=False, reftime=None)

    assert ts.colnames == ["time", "Satellite position", "Satellite velocity"]


def test_default_position_only_and_spherical_names_remain_unchanged(
    position_only, position_velocity
):
    position_ts = position_only.ephemeris(elapsed=False, reftime=None)
    spherical_ts = position_velocity.ephemeris(
        elapsed=False, reftime=None, coordinate_type="spherical"
    )

    assert position_ts.colnames == ["time", "position_x", "position_y", "position_z"]
    assert spherical_ts.colnames == [
        "time",
        "rtasc",
        "decl",
        "distance",
        "rtasc_r",
        "decl_r",
        "distance_r",
    ]
