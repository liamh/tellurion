import numpy as np
import pytest
import astropy.units as u

import tellurion as tell


pytestmark = pytest.mark.requires_orekit_data


EARTH_EME2000 = tell.FrameDescriptor("Earth", "EME2000")
SSB_ICRF = tell.FrameDescriptor("Solar System Barycenter", "ICRF")
EPOCH = tell.abstime("2025-01-01T00:00:00")
LEO = [5740.13268349, 3314.06715, 0.0, -2.75082684, 4.76457184, 5.50165367]


@pytest.fixture
def leo_eci():
    return tell.pvtcart(LEO, EPOCH, frame=EARTH_EME2000)


def test_generator_contents(leo_eci):
    gen = tell.mbprepare(leo_eci, 1 * u.hour, ["Sun", "earth", "Moon"])
    assert gen["frame"] == SSB_ICRF
    assert len(gen["bodies"]) == 3


def test_output_frame_tags(leo_eci):
    gen = tell.mbprepare(leo_eci, 1 * u.hour)
    steps = np.linspace(10 * u.minute, 60 * u.minute, 6)
    bary = tell.mbpropagate(gen, steps)
    eci = tell.mbpropagate(gen, steps, frame=EARTH_EME2000)
    assert bary.frame == SSB_ICRF
    assert eci.frame == EARTH_EME2000
    assert len(bary) == len(eci) == 6


def test_include_init_round_trips_to_input(leo_eci):
    gen = tell.mbprepare(leo_eci, 1 * u.hour)
    out = tell.mbpropagate(
        gen, [30 * u.minute], include_init=True, frame=EARTH_EME2000
    )
    np.testing.assert_allclose(
        out[0].position_vector.to_value(u.km), LEO[:3], atol=1e-6
    )
    np.testing.assert_allclose(
        out[0].velocity_vector.to_value(u.km / u.s), LEO[3:], atol=1e-9
    )


def test_cislunar_matches_two_body_over_short_arc(leo_eci):
    gen = tell.mbprepare(leo_eci, 1 * u.hour)
    mb = tell.mbpropagate(gen, [1 * u.hour], frame=EARTH_EME2000)
    kgen = tell.prepare(leo_eci, 1 * u.hour, propagator="keplerian")
    kep = tell.propagate(kgen, [1 * u.hour], output="pvt")
    final = lambda x: np.ravel(x.position_vector.to_value(u.km))[-3:]
    diff = np.linalg.norm(final(mb) - final(kep))
    assert diff < 1.0


def test_backward_propagation_returns_to_start(leo_eci):
    fwd = tell.mbprepare(leo_eci, 2 * u.hour)
    end = tell.mbpropagate(fwd, [2 * u.hour], frame=EARTH_EME2000)
    back = tell.mbprepare(end, -2 * u.hour)
    start = tell.mbpropagate(back, [-2 * u.hour], frame=EARTH_EME2000)
    np.testing.assert_allclose(
        np.ravel(start.position_vector.to_value(u.km)), LEO[:3], atol=1e-3
    )


def test_interplanetary_heliocentric_arc():
    au = 149597870.7
    vcirc = 29.78
    pvt = tell.pvtcart(
        [au, 0.0, 0.0, 0.0, vcirc, 0.0], EPOCH, frame=SSB_ICRF
    )
    gen = tell.mbprepare(pvt, 30 * u.day, ["Sun", "Jupiter"])
    out = tell.mbpropagate(gen, [30 * u.day])
    r = np.linalg.norm(out.position_vector.to_value(u.km))
    assert abs(r - au) / au < 0.05
    assert out.frame == SSB_ICRF


def test_errors(leo_eci):
    with pytest.raises(ValueError, match="Unknown celestial body"):
        tell.mbprepare(leo_eci, 1 * u.hour, ["Vulcan"])
    with pytest.raises(ValueError, match="At least one"):
        tell.mbprepare(leo_eci, 1 * u.hour, [])
    noframe = tell.pvtcart(LEO, EPOCH)
    with pytest.raises(ValueError, match="frame descriptor"):
        tell.mbprepare(noframe, 1 * u.hour)
