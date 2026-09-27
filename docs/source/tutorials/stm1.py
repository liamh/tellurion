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
# # State transition matrix (STM) tutorial
#
# This notebook demonstrates how to request and inspect STM output from propagation,
# then perform a basic finite-difference sensitivity check.
#

# %% [markdown]
# ## Setup
#

# %%
import numpy as np
import astropy.units as u

import tellurion as tell


# %% [markdown]
# ## Prepare propagation with STM enabled
#
# Set `events['stm'] = True` and use a numerical propagator with a force model.
#

# %%
initstate = tell.pvtcart(
    [5740.1326835, 3314.06715, 0.0, -2.7508268, 4.7645718, 5.5016537],
    tell.abstime("2025-01-01T00:00:00"),
)

proptime = 2.0 * u.hour
forceenv = tell.setgravity(20, 20)
events = {
    "altitude": 125.0 * u.km,
    "eclipse": [],
    "visibility": [],
    "stm": True,
}

gen = tell.prepare(
    initstate,
    proptime,
    events=events,
    forceenv=forceenv,
    propagator="numerical",
)

_ = tell.propagate(gen, [proptime], output="pvt")


# %% [markdown]
# ## Inspect STM output
#
# The final STM is stored in `gen['final']['stm']`.
#

# %%
phi = np.array(gen["final"]["stm"])
print("STM shape:", phi.shape)
print("Top-left 3x3 block:\n", phi[:3, :3])


# %% [markdown]
# ## Sensitivity check with finite differences
#
# We perturb the initial x-position by 1 meter, propagate both cases, and compare:
#
# - Finite-difference change in final Cartesian state
# - Linearized prediction from `phi @ delta_x0`
#

# %%
delta = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

cart0 = initstate.cartesian
cart_pert = np.concatenate(
    [
        cart0["position"].si.value + np.array([1.0, 0.0, 0.0]),
        cart0["velocity"].si.value,
    ]
)
init_pert = tell.pvtcart(cart_pert, initstate.time)

gen_nom = tell.prepare(initstate, proptime, events=events, forceenv=forceenv, propagator="numerical")
gen_pert = tell.prepare(init_pert, proptime, events=events, forceenv=forceenv, propagator="numerical")

final_nom = tell.propagate(gen_nom, [proptime], output="pvt")[-1]
final_pert = tell.propagate(gen_pert, [proptime], output="pvt")[-1]

dx_fd = np.concatenate(
    [
        final_pert.position_vector.si.value - final_nom.position_vector.si.value,
        final_pert.velocity_vector.si.value - final_nom.velocity_vector.si.value,
    ]
)
dx_lin = phi @ delta

print("Finite-difference state delta:")
print(dx_fd)
print("\nSTM linearized prediction:")
print(dx_lin)
print("\nDifference norm:", np.linalg.norm(dx_fd - dx_lin))


# %% [markdown]
# ## Practical guidance
#
# - STM is most useful for local sensitivity analysis, targeting, covariance transport,
#   and debugging trajectory response to initial-condition perturbations.
# - Enabling STM adds computational overhead; use it when sensitivity information is needed,
#   and disable it for pure ephemeris production runs.
#
