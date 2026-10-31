import numpy as np
import pytest
import astropy.units as u

import tellurion as tell


pytestmark = pytest.mark.requires_orekit_data


EARTH_EME2000 = tell.FrameDescriptor("Earth", "EME2000")
SSB_ICRF = tell.FrameDescriptor("Solar System Barycenter", "ICRF")


def test_resolve_supported_inertial_frames():
    assert tell.resolve_orekit_frame(EARTH_EME2000).getName() == "EME2000"
    assert tell.resolve_orekit_frame(SSB_ICRF).getName() == "ICRF"
    assert tell.resolve_orekit_frame(
        tell.FrameDescriptor("earth", "J2000")
    ).getName() == "EME2000"
    assert tell.resolve_orekit_frame(
        tell.FrameDescriptor("SSB", "ICRF")
    ).getName() == "ICRF"


@pytest.mark.parametrize(
    "descriptor",
    [
        tell.FrameDescriptor("Earth", "ICRF"),
        tell.FrameDescriptor("Moon", "EME2000"),
        tell.FrameDescriptor(
            "Earth", "ITRF", tell.FrameKind.BODY_FIXED
        ),
    ],
)
def test_reject_unsupported_frame_descriptors(descriptor):
    with pytest.raises(ValueError, match="Unsupported|not supported"):
        tell.resolve_orekit_frame(descriptor)


def test_transform_pvt_round_trip_scalar():
    state = tell.pvtcart(
        [7000.0, -1200.0, 500.0, 1.0, 7.0, -0.5],
        tell.abstime("2025-01-01T00:00:00"),
        frame=EARTH_EME2000,
    )

    barycentric = tell.transform_pvt(state, SSB_ICRF)
    restored = tell.transform_pvt(barycentric, EARTH_EME2000)

    assert barycentric.frame == SSB_ICRF
    assert np.linalg.norm(barycentric.position_vector.to_value(u.m)) > 1e10
    np.testing.assert_allclose(
        restored.position_vector.to_value(u.m),
        state.position_vector.to_value(u.m),
        rtol=0.0,
        atol=1e-2,
    )
    np.testing.assert_allclose(
        restored.velocity_vector.to_value(u.m / u.s),
        state.velocity_vector.to_value(u.m / u.s),
        rtol=0.0,
        atol=1e-5,
    )


def test_transform_pvt_round_trip_series():
    times = tell.abstime(
        ["2025-01-01T00:00:00", "2025-01-01T00:10:00"]
    )
    states = np.array(
        [
            [7000.0, -1200.0, 500.0, 1.0, 7.0, -0.5],
            [7010.0, -1150.0, 470.0, 0.9, 7.1, -0.4],
        ]
    )
    pvt = tell.pvtcart(states, times, frame=EARTH_EME2000)

    restored = tell.transform_pvt(
        tell.transform_pvt(pvt, SSB_ICRF), EARTH_EME2000
    )

    assert len(restored) == len(pvt)
    assert np.all(restored.time == pvt.time)
    np.testing.assert_allclose(
        restored.position_vector.to_value(u.m),
        pvt.position_vector.to_value(u.m),
        rtol=0.0,
        atol=1e-2,
    )
    np.testing.assert_allclose(
        restored.velocity_vector.to_value(u.m / u.s),
        pvt.velocity_vector.to_value(u.m / u.s),
        rtol=0.0,
        atol=1e-5,
    )


def test_transform_requires_frame_and_epoch():
    pvt = tell.pvtcart(
        [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0],
        tell.abstime("2025-01-01T00:00:00"),
    )
    with pytest.raises(ValueError, match="frame descriptor"):
        tell.transform_pvt(pvt, SSB_ICRF)

    no_epoch = tell.pvtcart(
        [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0],
        None,
        frame=EARTH_EME2000,
    )
    with pytest.raises(ValueError, match="epoch"):
        tell.transform_pvt(no_epoch, SSB_ICRF)
