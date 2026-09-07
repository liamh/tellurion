# test_iod.py
"""Tests for Lambert's problem IOD solver using orbit fixtures."""

import astropy.units as u
import numpy as np
import pytest

import tellurion as tell


class TestLambertSolver:
    """Test IodLambert solver against Keplerian propagated orbits."""

    @pytest.mark.parametrize("tof", [15 * u.min, 30 * u.min, 45 * u.min])
    def test_lambert_leo1(self, leo1, tof):
        """Propagate LEO circular orbit, solve Lambert, and verify initial velocity."""
        pvt_init = leo1["pvt"]

        gen = tell.prepare(pvt_init, tof, propagator="keplerian")
        pvt_final = tell.propagate(gen, [tof], output="pvt")[-1]

        p1 = pvt_init.position
        p2 = pvt_final.position

        pvt1_sol, pvt2_sol = tell.lambert(p1, p2)

        v_expected = pvt_init.velocity_vector.si.value
        v_solved = pvt1_sol.velocity_vector.si.value

        np.testing.assert_allclose(v_solved, v_expected, rtol=1e-5, atol=1e-3)

    @pytest.mark.parametrize("tof", [2 * u.hour, 6 * u.hour, 10 * u.hour])
    def test_lambert_geo1(self, geo1, tof):
        """Propagate GEO orbit, solve Lambert, and verify initial velocity."""
        pvt_init = geo1["pvt"]

        gen = tell.prepare(pvt_init, tof, propagator="keplerian")
        pvt_final = tell.propagate(gen, [tof], output="pvt")[-1]

        p1 = pvt_init.position
        p2 = pvt_final.position

        pvt1_sol, pvt2_sol = tell.lambert(p1, p2)

        v_expected = pvt_init.velocity_vector.si.value
        v_solved = pvt1_sol.velocity_vector.si.value

        np.testing.assert_allclose(v_solved, v_expected, rtol=1e-5, atol=1e-3)

    @pytest.mark.parametrize("tof", [1 * u.hour, 3 * u.hour, 5 * u.hour])
    def test_lambert_ell1(self, ell1, tof):
        """Propagate GTO elliptical orbit, solve Lambert,
        and verify initial velocity."""
        pvt_init = ell1["pvt"]

        gen = tell.prepare(pvt_init, tof, propagator="keplerian")
        pvt_final = tell.propagate(gen, [tof], output="pvt")[-1]

        p1 = pvt_init.position
        p2 = pvt_final.position

        pvt1_sol, pvt2_sol = tell.lambert(p1, p2)

        v_expected = pvt_init.velocity_vector.si.value
        v_solved = pvt1_sol.velocity_vector.si.value

        np.testing.assert_allclose(v_solved, v_expected, rtol=1e-5, atol=1e-3)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
