# test_prop.py
"""Propagation tests using shared orbit fixtures."""

import warnings

import astropy.units as u
import numpy as np
import pytest

import tellurion as tell

# test_prop.py
"""Propagation tests using shared orbit fixtures."""


@pytest.mark.parametrize("propagator,state_type", [
    ("sgp4", "mean"),
    ("numerical", "mean"),
    ("keplerian", "mean"),
    ("dsst", "invalid"),
])
def test_unsupported_state_type(propagator, state_type):
    with pytest.raises(ValueError, match="state"):
        tell.prepare(None, 1 * u.hour, propagator=propagator, state_type=state_type)


def test_sgp4_explicit_osculating():
    tle = (
        "1 25544U 98067A   24306.50000000  .00001234  00000-0  23456-4 0  9991",
        "2 25544  51.6436 150.1234 0001234  90.1234 270.0123 15.49567890123456",
    )
    initial = tell.MeanElementSetT({}, tell.abstime("2024-11-01T12:00:00"),
                                  tle, "SGP4", {})
    default = tell.prepare(initial, 10 * u.min)
    explicit = tell.prepare(initial, 10 * u.min, state_type="osculating")
    assert default["state_type"] == explicit["state_type"] == "osculating"
    assert np.linalg.norm((default["final"]["pvt"].position_vector -
                           explicit["final"]["pvt"].position_vector).to_value(u.m)) < 1e-6


class TestKeplerianPropagator:
    """Test Keplerian (two-body) propagator on various orbits."""

    def test_keplerian_leo1_pvt(self, leo1):
        """Propagate LEO circular orbit (PVT)."""
        initstate = leo1["pvt"]
        proptime = 1.0 * u.day

        gen = tell.prepare(initstate, proptime, propagator="keplerian")
        result = tell.propagate(gen, [0.5 * u.day])
        assert result is not None

    def test_keplerian_leo1_kep(self, leo1):
        """Propagate LEO circular orbit (Keplerian elements)."""
        initstate = leo1["kep"]
        proptime = 1.0 * u.day

        gen = tell.prepare(initstate, proptime, propagator="keplerian")
        result = tell.propagate(gen, [0.5 * u.day])
        assert result is not None

    def test_keplerian_all_orbits_pvt(self, all_orbits):
        """Keplerian propagator should work on all orbit types (PVT)."""
        proptime = 2.0 * u.day

        for orbit in all_orbits:
            initstate = orbit["pvt"]
            gen = tell.prepare(initstate, proptime, propagator="keplerian")
            result = tell.propagate(gen, [1.0 * u.day])
            assert result is not None

    def test_keplerian_all_orbits_kep(self, all_orbits):
        """Keplerian propagator should work on all orbit types (Keplerian elements)."""
        proptime = 2.0 * u.day

        for orbit in all_orbits:
            initstate = orbit["kep"]
            gen = tell.prepare(initstate, proptime, propagator="keplerian")
            result = tell.propagate(gen, [1.0 * u.day])
            assert result is not None


class TestBrouwerLyddanePropagator:
    """Test Brouwer-Lyddane (J2-J5 perturbations) propagator."""

    def test_state_type_and_initial_state(self, ell2):
        from tellurion.ork import convert
        from org.orekit.propagation.analytical import BrouwerLyddanePropagator

        forceenv = tell.setgravity(5, 0)
        gens = {
            kind: tell.prepare(ell2["pvt"], 30 * u.min, forceenv=forceenv,
                               propagator="brouwer-lyddane", state_type=kind)
            for kind in ("osculating", "mean")
        }
        initial = tell.propagate(gens["osculating"], 0 * u.s,
                                 include_init=False, output="pvt")
        assert np.linalg.norm((initial.position_vector - ell2["pvt"].position_vector)
                              .to_value(u.m)) < 10

        date = gens["mean"]["epoch"].shiftedBy(600.0)
        physical = gens["mean"]["osculating propfn"](date)
        provider = forceenv["gravity-unnorm"]
        expected = BrouwerLyddanePropagator.computeMeanOrbit(
            physical.getOrbit(), provider, provider.onDate(date), 0.0, 1e-5, 500
        )
        actual = gens["mean"]["propfn"](date)
        assert actual.getPVCoordinates().getPosition().distance(
            expected.getPVCoordinates().getPosition()
        ) < 1e-6
        assert actual.getPVCoordinates().getPosition().distance(
            physical.getPVCoordinates().getPosition()
        ) > 1.0
        assert np.linalg.norm((gens["mean"]["final"]["pvt"].position_vector -
                               convert._pvt(gens["mean"]["final"]["state"]).position_vector)
                              .to_value(u.m)) < 1e-6

    def test_mean_output_uses_osculating_event_states(self, ell2, monkeypatch):
        from tellurion.ork.event import prop as event

        gen = tell.prepare(ell2["pvt"], 30 * u.min, forceenv=tell.setgravity(5, 0),
                           propagator="brouwer-lyddane", state_type="mean")
        captured = []
        original = event._evstates

        def capture(generator, states):
            captured.extend(states if isinstance(states, list) else [states])
            return original(generator, states)

        monkeypatch.setattr(event, "_evstates", capture)
        tell.propagate(gen, 10 * u.min, include_init=False, output="pvt")
        date = gen["epoch"].shiftedBy(600)
        assert captured[0].getPosition().distance(
            gen["osculating propfn"](date).getPosition()
        ) < 1e-6
        assert captured[0].getPosition().distance(gen["propfn"](date).getPosition()) > 1

    def test_bl_ell2_from_kep(self, ell2):
        """Propagate moderately elliptical orbit (Keplerian elements)."""
        initstate = ell2["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_bl_ell2_from_pvt(self, ell2):
        """Propagate moderately elliptical orbit (PVT)."""
        initstate = ell2["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_bl_leo1_from_kep(self, leo1):
        """Propagate circular orbit (Keplerian elements)."""
        initstate = leo1["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_bl_leo1_from_pvt(self, leo1):
        """Propagate circular orbit (PVT)."""
        initstate = leo1["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

# See post https://forum.orekit.org/t/unable-to-compute-brouwer-lyddane-mean-parameters/5375/5?u=liamh
    @pytest.mark.xfail(reason="Orekit BrouwerLyddanePropagator convergence"
                       "failure on leo2 for unknown reasons")
    def test_bl_leo2_from_kep(self, leo2):
        """Propagate moderately elliptical orbit (Keplerian elements)."""
        initstate = leo2["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

# See post https://forum.orekit.org/t/unable-to-compute-brouwer-lyddane-mean-parameters/5375/5?u=liamh
    @pytest.mark.xfail(reason="Orekit BrouwerLyddanePropagator convergence"
                       "failure on leo2 for unknown reasons")
    def test_bl_leo2_from_pvt(self, leo2):
        """Propagate moderately elliptical orbit (PVT)."""
        initstate = leo2["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_bl_ell1_from_kep(self, ell1):
        """Propagate highly elliptical orbit (Keplerian elements)."""
        initstate = ell1["kep"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_bl_ell1_from_pvt(self, ell1):
        """Propagate highly elliptical orbit (PVT)."""
        initstate = ell1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_bl_vang1_from_kep(self, vang1):
        """Propagate Vanguard orbit (Keplerian elements)."""
        initstate = vang1["kep"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_bl_vang1_from_pvt(self, vang1):
        """Propagate Vanguard orbit (PVT)."""
        initstate = vang1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_bl_vs_keplerian_ell2_kep(self, ell2):
        """Brouwer-Lyddane differs from Keplerian (Keplerian elements)."""
        initstate = ell2["kep"]
        proptime = 10.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen_kep = tell.prepare(initstate, proptime, propagator="keplerian")
        gen_bl = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )

        times = [5.0 * u.day]
        result_kep = tell.propagate(gen_kep, times, output="pvt")
        result_bl = tell.propagate(gen_bl, times, output="pvt")

        # Results should differ due to J2-J5 perturbations
        assert not np.allclose(
            result_kep.cartesian["position"].value,
            result_bl.cartesian["position"].value,
        )

    def test_bl_vs_keplerian_ell2_pvt(self, ell2):
        """Brouwer-Lyddane differs from Keplerian (PVT)."""
        initstate = ell2["pvt"]
        proptime = 10.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen_kep = tell.prepare(initstate, proptime, propagator="keplerian")
        gen_bl = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="brouwer-lyddane"
        )

        times = [5.0 * u.day]
        result_kep = tell.propagate(gen_kep, times, output="pvt")
        result_bl = tell.propagate(gen_bl, times, output="pvt")

        # Results should differ due to J2-J5 perturbations
        assert not np.allclose(
            result_kep.cartesian["position"].value,
            result_bl.cartesian["position"].value,
        )


class TestNumericalPropagator:
    """Test numerical integration propagator."""

    def test_numerical_leo1_pvt(self, leo1):
        """Propagate LEO orbit with full gravity model (PVT)."""
        initstate = leo1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(
        reason="Circular orbit (e=0) has singular Jacobian in Keplerian coordinates"
    )
    def test_numerical_leo1_kep(self, leo1):
        """Propagate LEO orbit with full gravity model (Keplerian elements)."""
        initstate = leo1["kep"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_numerical_leo1_circular(self, leo1):
        """Propagate LEO orbit with circular elements (non-singular)."""
        initstate = leo1["pvt"].circular()
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_numerical_leo2_kep(self, leo2):
        """Propagate moderately elliptical orbit (Keplerian elements)."""
        initstate = leo2["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_leo2_pvt(self, leo2):
        """Propagate moderately elliptical orbit (PVT)."""
        initstate = leo2["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_ell2_pvt(self, ell2):
        """Propagate elliptical orbit with full gravity model (PVT)."""
        initstate = ell2["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_ell2_kep(self, ell2):
        """Propagate elliptical orbit with full gravity model (Keplerian elements)."""
        initstate = ell2["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_geo1_pvt(self, geo1):
        """Propagate GEO orbit with full gravity model (PVT)."""
        initstate = geo1["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(
        reason="Circular orbit (e=0) has singular Jacobian in Keplerian coordinates"
    )
    def test_numerical_geo1_kep(self, geo1):
        """Propagate GEO orbit with full gravity model (Keplerian elements)."""
        initstate = geo1["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(
        reason="Circular elements cannot represent equatorial orbit (i=0)"
    )
    def test_numerical_geo1_circular(self, geo1):
        """Propagate GEO orbit with circular elements (fails for equatorial)."""
        initstate = geo1["pvt"].circular()
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_geo1_equinoctial(self, geo1):
        """Propagate GEO orbit with equinoctial elements (handles equatorial)."""
        initstate = geo1["pvt"].equinoctial()
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_numerical_vs_keplerian(self, leo1):
        """Numerical should differ from Keplerian due to higher-order gravity."""
        initstate = leo1["pvt"]
        proptime = 5.0 * u.day

        gen_kep = tell.prepare(initstate, proptime, propagator="keplerian")
        gen_num = tell.prepare(
            initstate,
            proptime,
            forceenv=tell.setgravity(20, 20),
            propagator="numerical",
        )

        times = [2.0 * u.day]
        result_kep = tell.propagate(gen_kep, times, output="pvt")
        result_num = tell.propagate(gen_num, times, output="pvt")

        # Results should differ due to higher-order gravity terms
        assert not np.allclose(
            result_kep.cartesian["position"].value,
            result_num.cartesian["position"].value,
        )


class TestEclipseDetection:
    """Test eclipse detection with various propagators."""

    def test_keplerian_leo1_eclipse(self, leo1):
        """Keplerian propagator with eclipse detection (LEO)."""
        initstate = leo1["pvt"]
        proptime = 8.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}

        gen = tell.prepare(initstate, proptime, events=events, propagator="keplerian")

        # Check that sun transitions were recorded
        assert "sun transition" in gen
        sun_trans = gen["sun transition"]
        assert len(sun_trans) > 0, "Should detect eclipse transitions in LEO"

        # Propagate and verify ephemeris has sunlight status
        result = tell.propagate(
            gen, np.linspace(1.0 * u.hour, 6.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_keplerian_vang1_eclipse(self, vang1):
        """Keplerian propagator with eclipse detection (Vanguard)."""
        initstate = vang1["pvt"]
        proptime = 24.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}

        gen = tell.prepare(initstate, proptime, events=events, propagator="keplerian")

        # Vanguard has high eccentricity, may have fewer eclipse events
        assert "sun transition" in gen
        result = tell.propagate(
            gen, np.linspace(4.0 * u.hour, 20.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_bl_leo1_eclipse(self, leo1):
        """Brouwer-Lyddane propagator with eclipse detection (LEO)."""
        initstate = leo1["pvt"]
        proptime = 8.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="brouwer-lyddane",
        )

        assert "sun transition" in gen
        sun_trans = gen["sun transition"]
        assert len(sun_trans) > 0, "Should detect eclipse transitions in LEO"

        result = tell.propagate(
            gen, np.linspace(1.0 * u.hour, 6.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_bl_vang1_eclipse(self, vang1):
        """Brouwer-Lyddane propagator with eclipse detection (Vanguard)."""
        initstate = vang1["pvt"]
        proptime = 24.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(5, 0)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="brouwer-lyddane",
        )

        assert "sun transition" in gen
        result = tell.propagate(
            gen, np.linspace(4.0 * u.hour, 20.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_numerical_leo1_eclipse(self, leo1):
        """Numerical propagator with eclipse detection (LEO)."""
        initstate = leo1["pvt"]
        proptime = 8.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="numerical",
        )

        assert "sun transition" in gen
        sun_trans = gen["sun transition"]
        assert len(sun_trans) > 0, "Should detect eclipse transitions in LEO"

        result = tell.propagate(
            gen, np.linspace(1.0 * u.hour, 6.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_numerical_ell2_eclipse(self, ell2):
        """Numerical propagator with eclipse detection (elliptical)."""
        initstate = ell2["pvt"]
        proptime = 24.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="numerical",
        )

        assert "sun transition" in gen
        result = tell.propagate(
            gen, np.linspace(4.0 * u.hour, 20.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_eclipse_sunlight_values(self, leo1):
        """Verify eclipse events produce expected sunlight status values."""
        initstate = leo1["pvt"]
        proptime = 8.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}

        gen = tell.prepare(initstate, proptime, events=events, propagator="keplerian")
        result = tell.propagate(
            gen, np.linspace(1.0 * u.hour, 6.0 * u.hour, 16), output="pvt"
        )

        # Sunlight status should be one of: 'u' (umbra), 'p' (penumbra), 's' (full sun)
        assert "eclipse" in result.aux and isinstance(result.aux["eclipse"], str)

        # Extract distinct non-space characters
        sunlight_statuses = set(result.aux["eclipse"].replace(" ", ""))
        valid_statuses = {"u", "p", "s"}

        assert sunlight_statuses.issubset(valid_statuses), (
            f"Unexpected sunlight statuses: {sunlight_statuses - valid_statuses}"
        )


class TestDSSSTPropagator:
    """Test DSST (Draper Semi-analytical Satellite Theory) propagator."""

    def test_mean_early_stop_matches_osculating(self):
        from tellurion.ork import convert

        descending = tell.kepler(
                tell.allplane(
                    {
                        "altper": 160 * u.km,
                        "altapo": 20250 * u.km,
                        "inc": 28.5 * u.deg,
                        "argper": 0 * u.deg,
                        "raan": 0 * u.deg,
                        "ma": 180 * u.deg,
                    }
                ),
                tell.abstime("2026-01-01 08:05:00"),
        ).pvt()
        proptime = 4 * u.hour
        events = {"altitude": 5000 * u.km, "eclipse": [], "visibility": []}
        gens = {
            kind: tell.prepare(
                descending,
                proptime,
                events=events,
                forceenv=tell.setgravity(2, 0),
                propagator="dsst",
                state_type=kind,
            )
            for kind in ("mean", "osculating")
        }

        mean_final = gens["mean"]["final"]["state"]
        osculating_final = gens["osculating"]["final"]["state"]
        assert mean_final.getDate().equals(osculating_final.getDate())
        assert gens["mean"]["mindate"].equals(gens["mean"]["epoch"])
        assert gens["mean"]["maxdate"].equals(mean_final.getDate())
        assert gens["osculating"]["maxdate"].equals(osculating_final.getDate())

        mean_pvt = gens["mean"]["final"]["pvt"]
        assert np.linalg.norm(
            (
                mean_pvt.position_vector
                - convert._pvt(mean_final).position_vector
            ).to_value(u.m)
        ) < 1e-6
        assert mean_final.getDate().durationFrom(gens["mean"]["epoch"]) < (
            proptime.to_value(u.s)
        )

    def test_state_types_and_physical_events(self, leo1, monkeypatch):
        from tellurion.ork.event import prop as event

        forceenv = tell.setgravity(2, 0)
        events = {"altitude": 125 * u.km, "eclipse": [True, True], "visibility": []}
        mean = tell.prepare(leo1["pvt"], 100 * u.min, events=events,
                            forceenv=forceenv, propagator="dsst", state_type="mean")
        osc = tell.prepare(leo1["pvt"], 100 * u.min, events=events,
                           forceenv=forceenv, propagator="dsst", state_type="osculating")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", FutureWarning)
            default = tell.prepare(leo1["pvt"], 100 * u.min, events=events,
                                   forceenv=forceenv, propagator="dsst")
        assert not any(isinstance(w.message, FutureWarning) and
                       "state_type" in str(w.message) for w in caught)
        assert mean["state_type"] == "mean"
        assert osc["state_type"] == "osculating"
        assert default["state_type"] == "mean"
        date = mean["epoch"].shiftedBy(600)
        mean_state = mean["propfn"](date)
        physical = mean["osculating propfn"](date)
        assert mean_state.getPosition().distance(physical.getPosition()) > 1
        assert mean_state.getPosition().distance(default["propfn"](date).getPosition()) < 1e-4
        assert physical.getPosition().distance(osc["propfn"](date).getPosition()) < 1e-4
        initial = tell.propagate(osc, 0 * u.s, include_init=False, output="pvt")
        assert np.linalg.norm((initial.position_vector - leo1["pvt"].position_vector)
                              .to_value(u.m)) < 0.1

        captured = []
        original = event._evstates

        def capture(generator, states):
            captured.extend(states if isinstance(states, list) else [states])
            return original(generator, states)

        monkeypatch.setattr(event, "_evstates", capture)
        mean_pvt = tell.propagate(mean, [10 * u.min, 40 * u.min],
                                  include_init=False, output="pvt")
        osc_pvt = tell.propagate(osc, [10 * u.min, 40 * u.min],
                                 include_init=False, output="pvt")
        assert captured[0].getPosition().distance(physical.getPosition()) < 1e-6
        assert mean_pvt.aux == osc_pvt.aux
        assert mean["sun transition"].time.size == osc["sun transition"].time.size

    @pytest.mark.parametrize("state_type", ["mean", "osculating"])
    def test_dsst_zero_degree(self, leo1, state_type):
        gen = tell.prepare(leo1["pvt"], 20 * u.min, forceenv=tell.setgravity(0, 0),
                           propagator="dsst", state_type=state_type)
        assert tell.propagate(gen, 10 * u.min, include_init=False, output="pvt") is not None

    def test_dsst_osculating_against_numerical(self, leo1):
        forceenv = tell.setgravity(2, 0)
        dsst = tell.prepare(leo1["pvt"], 1 * u.hour, forceenv=forceenv,
                            propagator="dsst", state_type="osculating")
        numerical = tell.prepare(leo1["pvt"], 1 * u.hour, forceenv=forceenv,
                                 propagator="numerical")
        actual = tell.propagate(dsst, 1 * u.hour, include_init=False, output="pvt")
        reference = tell.propagate(numerical, 1 * u.hour, include_init=False, output="pvt")
        assert np.linalg.norm((actual.position_vector - reference.position_vector)
                              .to_value(u.m)) < 1000

    def test_dsst_drag_force(self, leo1):
        forceenv = tell.dragforce(tell.setgravity(2, 0))
        gen = tell.prepare(leo1["pvt"], 10 * u.min, forceenv=forceenv,
                           propagator="dsst", state_type="osculating")
        assert tell.propagate(gen, 5 * u.min, include_init=False, output="pvt") is not None

    def test_dsst_leo1_pvt(self, leo1):
        """Propagate LEO orbit with DSST (PVT)."""
        initstate = leo1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(reason="Circular orbit (e=0) has singular"
                       "Jacobian in Keplerian coordinates")
    def test_dsst_leo1_kep(self, leo1):
        """Propagate LEO orbit with DSST (Keplerian elements)."""
        initstate = leo1["kep"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_dsst_leo1_circular(self, leo1):
        """Propagate LEO orbit with DSST and circular elements (non-singular)."""
        initstate = leo1["pvt"].circular()
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [2.0 * u.day])
        assert result is not None

    def test_dsst_ell2_pvt(self, ell2):
        """Propagate elliptical orbit with DSST (PVT)."""
        initstate = ell2["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_dsst_ell2_kep(self, ell2):
        """Propagate elliptical orbit with DSST (Keplerian elements)."""
        initstate = ell2["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None


    def test_dsst_geo1_pvt(self, geo1):
        """Propagate GEO orbit with DSST (PVT)."""
        initstate = geo1["pvt"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(reason="Circular orbit (e=0) has singular"
                       "Jacobian in Keplerian coordinates")
    def test_dsst_geo1_kep(self, geo1):
        """Propagate GEO orbit with DSST (Keplerian elements)."""
        initstate = geo1["kep"]
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    @pytest.mark.xfail(
        reason="Circular elements cannot represent equatorial orbit (i=0)"
    )
    def test_dsst_geo1_circular(self, geo1):
        """Propagate GEO orbit with DSST and circular elements
        (fails for equatorial)."""
        initstate = geo1["pvt"].circular()
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_dsst_geo1_equinoctial(self, geo1):
        """Propagate GEO orbit with DSST and equinoctial elements
        (handles equatorial)."""
        initstate = geo1["pvt"].equinoctial()
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [3.0 * u.day])
        assert result is not None

    def test_dsst_vang1_pvt(self, vang1):
        """Propagate Vanguard orbit with DSST (PVT)."""
        initstate = vang1["pvt"]
        proptime = 10.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [5.0 * u.day])
        assert result is not None

    def test_dsst_vang1_kep(self, vang1):
        """Propagate Vanguard orbit with DSST (Keplerian elements)."""
        initstate = vang1["kep"]
        proptime = 10.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [5.0 * u.day])
        assert result is not None

    def test_dsst_vs_numerical(self, leo1):
        """DSST and numerical propagators should give similar results (LEO)."""
        initstate = leo1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen_num = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        gen_dsst = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="dsst"
        )

        times = [2.0 * u.day]
        result_num = tell.propagate(gen_num, times, output="pvt")
        result_dsst = tell.propagate(gen_dsst, times, output="pvt")

        # Results should be reasonably close (within a few km for LEO)
        # DSST is semi-analytical so not exact match, but much closer than Keplerian
        position_diff = np.linalg.norm(
            result_num.cartesian["position"].value
            - result_dsst.cartesian["position"].value
        )
        assert position_diff < 1000.0, (
            f"DSST vs numerical position diff: {position_diff} m"
        )

    def test_dsst_vs_numerical_ell2(self, ell2):
        """DSST and numerical propagators should give similar results (elliptical)."""
        initstate = ell2["pvt"]
        proptime = 10.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen_num = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="numerical"
        )
        gen_dsst = tell.prepare(
            initstate, proptime, forceenv=forceenv, propagator="dsst"
        )

        times = [5.0 * u.day]
        result_num = tell.propagate(gen_num, times, output="pvt")
        result_dsst = tell.propagate(gen_dsst, times, output="pvt")

        position_diff = np.linalg.norm(
            result_num.cartesian["position"].value
            - result_dsst.cartesian["position"].value
        )
        assert position_diff < 5000.0, (
            f"DSST vs numerical position diff: {position_diff} m"
        )

    def test_dsst_leo1_eclipse(self, leo1):
        """DSST propagator with eclipse detection (LEO)."""
        initstate = leo1["pvt"]
        proptime = 8.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="dsst",
        )

        assert "sun transition" in gen
        sun_trans = gen["sun transition"]
        assert len(sun_trans) > 0, "Should detect eclipse transitions in LEO"

        result = tell.propagate(
            gen, np.linspace(1.0 * u.hour, 6.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_dsst_vang1_eclipse(self, vang1):
        """DSST propagator with eclipse detection (Vanguard)."""
        initstate = vang1["pvt"]
        proptime = 24.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="dsst",
        )

        assert "sun transition" in gen
        result = tell.propagate(
            gen, np.linspace(4.0 * u.hour, 20.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_dsst_ell2_eclipse(self, ell2):
        """DSST propagator with eclipse detection (elliptical)."""
        initstate = ell2["pvt"]
        proptime = 24.0 * u.hour
        events = {"altitude": 125.0 * u.km, "eclipse": [True, True], "visibility": []}
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(
            initstate,
            proptime,
            events=events,
            forceenv=forceenv,
            propagator="dsst",
        )

        assert "sun transition" in gen
        result = tell.propagate(
            gen, np.linspace(4.0 * u.hour, 20.0 * u.hour, 16), output="pvt"
        )
        assert result is not None
        assert (
            "eclipse" in result.aux
            and isinstance(result.aux["eclipse"], str)
            and len(result.aux["eclipse"]) == 33
        )

    def test_dsst_no_eclipse_events(self, leo1):
        """DSST propagator without eclipse detection."""
        initstate = leo1["pvt"]
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime, forceenv=forceenv, propagator="dsst")
        result = tell.propagate(gen, [2.0 * u.day], output="pvt")

        assert result is not None
        # Without eclipse events, aux should not have eclipse key
        assert "eclipse" not in result.aux


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
