# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Lambert transfer tutorial
#
# This notebook demonstrates practical use of `tell.lambert` for a single-revolution transfer,
# validation against a propagated endpoint, and short-way versus long-way solutions.
#

# %% [markdown]
# ## Setup
#

# %%
import numpy as np
import astropy.units as u

import tellurion as tell


# %% [markdown]
# ## Build an initial state and a target endpoint
#
# We build a reference endpoint by propagating a known initial state for a fixed time of flight.
#

# %%
pvt0 = tell.pvtcart(
    [5740.1326835, 3314.06715, 0.0, -2.7508268, 4.7645718, 5.5016537],
    tell.abstime("2025-01-01T00:00:00"),
)

tof = 45.0 * u.minute
truth_gen = tell.prepare(pvt0, tof, propagator="keplerian")
pvtf_truth = tell.propagate(truth_gen, [tof], output="pvt")[-1]


# %% [markdown]
# ## Solve Lambert's problem (single revolution)
#

# %%
pvt1_short, pvt2_short = tell.lambert(
    pvt0.position,
    pvtf_truth.position,
    shortway=True,
    n_rev=0,
)

print("Solved initial velocity [m/s]:", pvt1_short.velocity_vector.si.value)
print("Solved final velocity [m/s]:", pvt2_short.velocity_vector.si.value)


# %% [markdown]
# ## Validate against propagated endpoint
#
# Propagate the solved initial Lambert state and compare final position with the reference endpoint.
#

# %%
val_gen = tell.prepare(pvt1_short, tof, propagator="keplerian")
pvtf_check = tell.propagate(val_gen, [tof], output="pvt")[-1]

pos_err_m = np.linalg.norm(
    pvtf_check.position_vector.si.value - pvtf_truth.position_vector.si.value
)
print(f"Endpoint position error: {pos_err_m:.6f} m")


# %% [markdown]
# ## Short-way vs long-way comparison
#
# Different transfer branches usually imply different velocity solutions.
#

# %%
pvt1_long, pvt2_long = tell.lambert(
    pvt0.position,
    pvtf_truth.position,
    shortway=False,
    n_rev=0,
)

dv_init = np.linalg.norm(
    pvt1_long.velocity_vector.si.value - pvt1_short.velocity_vector.si.value
)
print(f"Initial-velocity branch difference: {dv_init:.6f} m/s")


# %% [markdown]
# ## Edge cases and limitations
#
# - `tell.lambert` requires two `PositionT` inputs with valid epochs.
# - For `n_rev > 0`, the current Orekit branch available through Tellurion does not return
#   both short-period and long-period multi-revolution solutions.
# - Near-degenerate geometries (very short TOF, nearly collinear points, inconsistent geometry)
#   can fail to converge or produce physically unsuitable transfers.
#
