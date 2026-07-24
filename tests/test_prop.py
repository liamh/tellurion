# test_prop.py
"""Propagation tests using shared orbit fixtures."""

import astropy.units as u
import numpy as np
import pytest
import tellurion as tell


# test_prop.py
"""Propagation tests using shared orbit fixtures."""

import astropy.units as u
import numpy as np
import pytest
import tellurion as tell


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
