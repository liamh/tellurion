# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Run this tutorial in Colab
#
# [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](
# https://colab.research.google.com/github/liamh/tellurion/blob/colab-notebooks/orbital_state_conversion.ipynb
# )
#
# If you want to edit and save, use **File → Save a copy in Drive**.

# %%
import sys
import subprocess

if "google.colab" in sys.modules:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "-q",
            "tellurion[data] @ git+https://github.com/liamh/tellurion.git@develop",
        ],
        check=True,
    )

# %% [markdown]
# # Representation and conversion of orbital state
#
# This tutorial shows how to create a {py:class}`~tellurion.PositionVelocityT`
# in Cartesian or spherical form and convert it to {py:class}`~tellurion.ElementSetT`,
# and vice versa.
#
# ## Setup
#
# We import Tellurion plus a few helper packages. {py:mod}`munch` is used as a
# convenient container for grouping related values in each example section.

# %%
import numpy as np
import astropy.units as u
from munch import Munch
import tellurion as tell

# %% [markdown]
# ## Position, velocity, and time (PVT)
#
# We define a low-Earth-orbit satellite from Cartesian position/velocity and an epoch.
# The state vector is a length-6 NumPy array:
# `[x, y, z, vx, vy, vz]`.
#
# For epoch we use {py:func}`~tellurion.abstime` with argument `0`, which means "now".

# %%
leo2 = Munch()
leo2.pv_cartesian = np.array(
    [4542.2829, 5170.0572, 565.0922, -5.2368, 4.1992, 3.5743]
)
leo2.epoch = tell.abstime(0)

# %% [markdown]
# Create a {py:class}`~tellurion.PositionVelocityT` with
# {py:func}`~tellurion.pvtcart`.
#
# The default preferred units are kilometers for position and kilometers/second
# for velocity (see `tell.prefunits` later).

# %%
leo2.pvt = tell.pvtcart(leo2.pv_cartesian, leo2.epoch)
leo2.pvt.cartesian

# %%
type(leo2.pvt.cartesian) is u.Quantity

# %%
leo2.pvt.cartesian.unit

# %%
leo2.pvt.cartesian.dtype

# %% [markdown]
# You can extract:
# - `position` (3-vector),
# - `velocity` (3-vector),
# - or the full numerical array.

# %%
leo2.pvt.cartesian["position"]

# %%
leo2.pvt.cartesian["velocity"]

# %%
leo2.pvt.cartesian.to_array()

# %% [markdown]
# Spherical coordinates are available from
# {py:attr}`tellurion.PositionBase.spherical`.
#
# For readability, use {py:func}`~tellurion.quantity_to_dict` (also available
# as `.to_dict()` on returned structured quantities).

# %%
leo2.pvt.spherical

# %%
leo2.pvt.spherical.to_dict()

# %%
leo2.pvt.spherical["distance"]

# %% [markdown]
# ## Quantities and units
#
# Units are provided by AstroPy (`astropy.units`).

# %%
157.0 * u.meter

# %%
leo2.pvt.cartesian.si

# %% [markdown]
# The package-wide preferred units can be inspected here:

# %%
tell.prefunits

# %% [markdown]
# ## Orbital elements
#
# ### Define an `ElementSetT`
#
# Create an element set directly with constructor functions like
# {py:func}`~tellurion.kepler`.

# %%
kep = tell.kepler(
    {
        "sma": 8000.0 * u.km,
        "ecc": 0.1,
        "inc": 42.0 * u.deg,
        "argper": 66.0 * u.deg,
        "raan": 217.4 * u.deg,
        "ma": 7.25 * u.deg,
    },
    tell.abstime("2023-09-14T08:30:00"),
)
kep.elements.to_dict()

# %%
kep.time

# %% [markdown]
# Convenience helpers can simplify common setups.
#
# Here:
# - {py:func}`~tellurion.allplane` converts perigee/apogee altitudes into
#   semimajor axis and eccentricity,
# - {py:func}`~tellurion.sma` converts a period-like input into orbital size
#   (here used with `altitude=True`).

# %%
ell1 = Munch()
ell1.els = tell.allplane(
    {
        "altper": 350 * u.km,
        "altapo": tell.sma(1.0, altitude=True),
        "inc": 0.0 * u.deg,
        "argper": 120.0 * u.deg,
        "raan": 0.0 * u.deg,
        "ma": 90.0 * u.deg,
    }
)
ell1.est = tell.kepler(ell1.els, tell.abstime("2026-01-01 05:55:00"))
ell1.est.elements.to_dict()

# %% [markdown]
# ### Convert element set → PVT
#
# Convert an {py:class}`~tellurion.ElementSetT` to
# {py:class}`~tellurion.PositionVelocityT` with
# {py:meth}`tellurion.ElementSetT.pvt`.

# %%
ell1.pvt = ell1.est.pvt()
ell1.pvt.cartesian.to_dict()

# %% [markdown]
# ### Convert PVT → element sets
#
# Starting from {py:class}`~tellurion.PositionVelocityT`, compute elements with:
# - {py:meth}`tellurion.PositionVelocityT.kepler`,
# - {py:meth}`tellurion.PositionVelocityT.equinoctial`,
# - {py:meth}`tellurion.PositionVelocityT.circular`.

# %%
leo2.pvt.kepler()

# %%
leo2.pvt.kepler().elements.to_dict()

# %%
leo2.pvt.kepler().elements["sma"]

# %%
leo2.pvt.equinoctial().elements.to_dict()

# %% [markdown]
# ## Individual orbital elements
#
# You can access specific elements directly (or from converted element sets)
# without always materializing every representation.
#
# Common Keplerian element abbreviations:
#
# | Abbreviation | Name | Physical dimension |
# | --- | --- | --- |
# | sma | semimajor axis | length |
# | ecc | eccentricity | dimensionless |
# | inc | inclination | angle |
# | argper | argument of perigee | angle |
# | radper | radius of perigee | length |
# | radapo | radius of apogee | length |
# | altper | altitude of perigee | length |
# | altapo | altitude of apogee | length |
# | raan | right ascension of ascending node | angle |
# | ta | true anomaly | angle |
# | ma | mean anomaly | angle |
# | memo | mean motion | angular speed |
# | period | orbital period | time |

# %% [markdown]
# ## Ephemerides
#
# Multiple states can be stacked and converted to an ephemeris table
# (an AstroPy time-series-like table in this interface).
#
# In practice this usually comes from propagation. Here we use a prepared array.

# %%
eph1 = Munch()
eph1.arr = np.array(
    [
        [5740.1327, 3314.0672, -0.0, -2.7508, 4.7646, 5.5017],
        [4581.8158, 4512.2632, 1616.8266, -4.8914, 3.1416, 5.1664],
        [2865.4996, 5161.0454, 3036.8466, -6.4324, 1.1404, 4.2038],
        [801.2655, 5183.5367, 4088.4417, -7.1880, -0.9901, 2.7365],
        [-1359.9604, 4580.2090, 4646.5576, -7.0741, -2.9891, 0.9484],
        [-3358.3326, 3427.3201, 4647.3126, -6.1152, -4.6175, -0.9413],
        [-4956.7920, 1865.8682, 4094.2852, -4.4363, -5.6868, -2.7068],
        [-5968.5036, 83.2918, 3056.3846, -2.2432, -6.0784, -4.1424],
        [-6277.2076, -1709.2532, 1658.3471, 0.2042, -5.7536, -5.0849],
        [-5848.9770, -3301.2156, 65.5519, 2.6219, -4.7547, -5.4287],
        [-4734.9345, -4506.1030, -1534.9324, 4.7319, -3.1982, -5.1357],
        [-3065.1825, -5182.0380, -2955.1853, 6.2903, -1.2623, -4.2384],
        [-1034.6459, -5247.7284, -4027.3431, 7.1126, 0.8306, -2.8370],
    ]
)

eph1.times = tell.abstime(
    [
        "2025-01-01T00:00:00.000",
        "2025-01-01T00:05:00.000",
        "2025-01-01T00:10:00.000",
        "2025-01-01T00:15:00.000",
        "2025-01-01T00:20:00.000",
        "2025-01-01T00:25:00.000",
        "2025-01-01T00:30:00.000",
        "2025-01-01T00:35:00.000",
        "2025-01-01T00:40:00.000",
        "2025-01-01T00:45:00.000",
        "2025-01-01T00:50:00.000",
        "2025-01-01T00:55:00.000",
        "2025-01-01T01:00:00.000",
    ]
)

eph1.pvts = tell.pvtcart(eph1.arr, eph1.times)
eph1.eph = eph1.pvts.ephemeris()
eph1.eph

# %% [markdown]
# As with a single PVT, component columns can be extracted; for ephemerides,
# these are 2D arrays over time.

# %%
eph1.eph["position"]

# %% [markdown]
# Compute time-series element columns from an ephemeris with
# {py:func}`~tellurion.tselements`.

# %%
eph1.altperapo = tell.tselements(eph1.eph, ["altper", "altapo"])
eph1.altperapo

# %%
eph1.altperapo["altper"]
