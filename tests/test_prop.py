# test_prop.py
"""Propagation tests using shared orbit fixtures."""

import numpy as np
import astropy.units as u
import pytest
import tellurion as tell


class TestKeplerianPropagator:
    """Test Keplerian (two-body) propagator on various orbits."""

    def test_keplerian_leo_circular(self, leo1):
        """Propagate ISS-like orbit."""
        initstate = leo1['pvt']
        proptime = 1.0 * u.day

        gen = tell.prepare(initstate, proptime, propagator='keplerian')
        result = tell.propagate(gen, [0.5*u.day])

        assert result is not None
        # Add specific assertions

    def test_keplerian_all_orbits(self, all_orbits):
        """Keplerian propagator should work on all orbit types."""
        proptime = 2.0 * u.day

        for orbit in all_orbits:
            initstate = orbit['pvt']
            gen = tell.prepare(initstate, proptime, propagator='keplerian')
            result = tell.propagate(gen, [1.0*u.day])
            assert result is not None


class TestBrouwerLyddanePropagator:
    """Test Brouwer-Lyddane (J2 perturbations) propagator."""

    def test_bl_leo_circular(self, leo1):
        """Propagate ISS-like orbit with J2 effects."""
        initstate = leo1['pvt']
        proptime = 7.0 * u.day
        forceenv = tell.setgravity(2,0)

        gen = tell.prepare(initstate, proptime, 
                          forceenv=forceenv,
                          propagator='brouwer-lyddane')
        result = tell.propagate(gen, [3.0*u.day])

        assert result is not None

    def test_bl_vs_keplerian(self, leo1):
        """Brouwer-Lyddane should differ from Keplerian over long periods."""
        initstate = leo1['pvt']
        proptime = 30.0 * u.day
        forceenv = tell.setgravity(2,0)

        gen_kep = tell.prepare(initstate, proptime, propagator='keplerian')
        gen_bl = tell.prepare(initstate, proptime, 
                             forceenv=forceenv,
                             propagator='brouwer-lyddane')

        times = [10.0*u.day]
        result_kep = tell.propagate(gen_kep, times)
        result_bl = tell.propagate(gen_bl, times)

        # Results should differ due to J2
        assert not np.allclose(result_kep.cartesian['position'].value,
                              result_bl.cartesian['position'].value)


class TestNumericalPropagator:
    """Test numerical integration propagator."""

    def test_numerical_leo(self, leo1):
        """Propagate with full gravity model."""
        initstate = leo1['pvt']
        proptime = 5.0 * u.day
        forceenv = tell.setgravity(20, 20)

        gen = tell.prepare(initstate, proptime,
                          forceenv=forceenv,
                          propagator='numerical')
        result = tell.propagate(gen, [2.0*u.day])

        assert result is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
