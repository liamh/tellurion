# %% [markdown]
# # Representation and conversion of orbital state
# This tutorial will show how to create a `PositionVelocityT` in Cartesian or spherical form and convert it to `ElementSetT`, and vice versa
#
# ## Setup
#
# First, some setup: we load tellurion and other packages.  Munch is a
# convenience for grouping calculations together.

# %%
import numpy as np
import astropy.units as u
from munch import Munch
import tellurion as tell

# %% [markdown]
# ## Position, velocity, time (PVT) in Cartesian and spherical coordinates
# Let's define a low-earth orbiting satellite by its Cartesian position and velocity, and then its epoch time. The position and velocity are given as length-6 {class}`numpy.array`. The epoch time doesn't matter for this tutorial, so we will make it whatever time it is now by giving {func}`~tellurion.abstime` the argument 0, meaning make an absolute time that is the same as the current time, i.e., now.
#
# <!--
# def keppvt(dictels, time):
#     """From the dictionary of element values, compute the Kepler
#     elements (sma, ecc, inc, argper, raan, and ma or ta), the PVT
#     transformation of those elements, and the spherical coordinates."""
#     d = {'els': dictels}
#     d['kep'] = tell.kepler(tell.allplane(d['els']), time)
#     d['pvt'] = tell.pvt(d['kep'])
#     d['pvt'].spherical
#     return d
# leo2 = keppvt({"altper": 525.0*u.km, "altapo": 555.0*u.km,
#                        "inc":28.5*u.deg, "argper": 40.0*u.deg,
#                        "raan": 40.0*u.deg, "ma": -30.0*u.deg},
#                       tell.abstime('2026-01-01 09:30:00'))
# leo2["pvt"].cartesian.to_array()
# -->

# %%
leo2 = Munch()
leo2.pv_cartesian = np.array([4542.2829, 5170.0572,  565.0922,   -5.2368,    4.1992,    3.5743])
leo2.epoch = tell.abstime(0)

# %% [markdown]
# To make a {class}`tell.PositionVelocityT <tellurion.PositionVelocityT>` (`PVT` for short), use the function {func}`tell.pvtcart() <tellurion.pvtcart>`. The [default units](#quantities-and-units) are kilometers for position and kilometers/second for velocity. To check we have the Cartesian orbital state correct, show the `.cartesian` property. This object has AstroPy's [structured units](https://docs.astropy.org/en/stable/units/structured_units.html), that is, the units are mixed.

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
# The individual components can be extracted by using the names `'position'` and `'velocity'`, or the entire vector as an `np.array`.

# %%
leo2.pvt.cartesian['position']

# %%
leo2.pvt.cartesian['velocity']

# %%
leo2.pvt.cartesian.to_array() # Make a numerical array

# %% [markdown]
# The spherical coordinates (right ascension, declination, and geocentric distance, and their rates) are found in the `.spherical` property. To make it easier to see, the function `.to_dict()` will convert the result into a dictionary, or the individual parts can be extracted directly.

# %%
leo2.pvt.spherical

# %%
leo2.pvt.spherical.to_dict()

# %%
leo2.pvt.spherical['distance']


# %% [markdown]
# <a id="quantities-and-units"></a>
# ## Quantities and units
#
# Quantities and units are defined using [AstroPy](https://docs.astropy.org/en/stable/units/index.html); the `import astropy.units as u` at the beginning permits specification of units. Quantities with units may be created or converted.

# %%
157.0*u.meter

# %%
leo2.pvt.cartesian.si

# %% [markdown]
# Default units are used throughout, so for most cases, values can be given as plain numbers

# %%
tell.prefunits

# %% [markdown]
# ## Orbital elements
# The available Keplerian orbital elements are:
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
# | raan | right ascension of the ascending node | angle |
# | ta | true anomaly | angle |
# | ma | mean anomaly | angle |
# | memo | mean motion | angular speed |
# | period | orbital period | time |
#

# %% [markdown]
# ### Defining an ElementSetT
# To define an element set directly from elements, use `tell.kepler()`.
# This is geosynchronous transfer orbit; using the [`tell.allplane()`](#tellurion.allplane) function that
# converts the altitudes of perigee and apogee to semimajor axis and eccentricity, and `tell.sma()` converts multiple quantities into a radial distance, in this case, the number of sidereal days into the radius of a circular orbit with that period.

# %%
ell1 = Munch()
ell1.els = tell.allplane({"altper": 350*u.km, "altapo": tell.sma(1.0, True),
               "inc":0.0*u.deg, "argper": 120.0*u.deg,
               "raan": 0.0*u.deg, "ma": 90.0*u.deg})
ell1.est = tell.kepler(ell1.els, tell.abstime('2026-01-01 05:55:00'))
# Show the elements as a dictionary
ell1.est.elements.to_dict()

# %% [markdown]
# This element set can be converted to a PVT

# %%
ell1.pvt = ell1.est.pvt()
ell1.pvt.cartesian.to_dict()

# %% [markdown]
# ### Convert from PVT
# Orbital elements can be generated from a PVT. The can be displayed more readably as a dictionary, or individual elements extracted

# %%
# This produces an object of class ElementSetT
leo2.pvt.kepler()

# %%
leo2.pvt.kepler().elements.to_dict()

# %%
leo2.pvt.kepler().elements['sma']

# %% [markdown]
# ## Ephemerides
# Multiple states can be stacked and then converted to an ephemeris table, an instance of an [AstroPy time series](https://docs.astropy.org/en/stable/timeseries/index.html). These usually come from orbit propagation; the example here is constructed a data array that originated in a propagation.

# %%
eph1 = Munch()
eph1.arr = np.array([[ 5740.1327,  3314.0672,    -0.    ,    -2.7508,     4.7646,     5.5017],
                [ 4581.8158,  4512.2632,  1616.8266,    -4.8914,     3.1416,     5.1664],
                [ 2865.4996,  5161.0454,  3036.8466,    -6.4324,     1.1404,     4.2038],
                [  801.2655,  5183.5367,  4088.4417,    -7.188 ,    -0.9901,     2.7365],
                [-1359.9604,  4580.209 ,  4646.5576,    -7.0741,    -2.9891,     0.9484],
                [-3358.3326,  3427.3201,  4647.3126,    -6.1152,    -4.6175,    -0.9413],
                [-4956.792 ,  1865.8682,  4094.2852,    -4.4363,    -5.6868,    -2.7068],
                [-5968.5036,    83.2918,  3056.3846,    -2.2432,    -6.0784,    -4.1424],
                [-6277.2076, -1709.2532,  1658.3471,     0.2042,    -5.7536,    -5.0849],
                [-5848.977 , -3301.2156,    65.5519,     2.6219,    -4.7547,    -5.4287],
                [-4734.9345, -4506.103 , -1534.9324,     4.7319,    -3.1982,    -5.1357],
                [-3065.1825, -5182.038 , -2955.1853,     6.2903,    -1.2623,    -4.2384],
                [-1034.6459, -5247.7284, -4027.3431,     7.1126,     0.8306,    -2.837 ]])
eph1.times = tell.abstime(['2025-01-01T00:00:00.000',
                      '2025-01-01T00:05:00.000',
                      '2025-01-01T00:10:00.000',
                      '2025-01-01T00:15:00.000',
                      '2025-01-01T00:20:00.000',
                      '2025-01-01T00:25:00.000',
                      '2025-01-01T00:30:00.000',
                      '2025-01-01T00:35:00.000',
                      '2025-01-01T00:40:00.000',
                      '2025-01-01T00:45:00.000',
                      '2025-01-01T00:50:00.000',
                      '2025-01-01T00:55:00.000',
                      '2025-01-01T01:00:00.000'])
eph1.pvts = tell.pvtcart(eph1.arr, eph1.times)
eph1.eph = eph1.pvts.ephemeris()
eph1.eph

# %% [markdown]
# As with a single PVT, components can be extracted, but this time we get a two-dimensional array.

# %%
eph1.eph['position']

# %% [markdown]
# Element values may be computed from an ephemeris table with
# `tell.tselements()`; for example, the altitudes of perigee and apogee. From this time series, an individual column may be extracted.

# %%
eph1.altperapo = tell.tselements(eph1.eph, ['altper','altapo'])
eph1.altperapo

# %%
eph1.altperapo['altper']
