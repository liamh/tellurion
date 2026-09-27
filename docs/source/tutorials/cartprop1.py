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

# %% [markdown] editable=true slideshow={"slide_type": ""}
# # Demonstration 

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Setup

# %% [markdown] editable=true slideshow={"slide_type": ""}
# First, it is necessary to import packages. A `Munch` is a `Dict` that acts like a structure, so we can use dot notation. We'll use this to demonstrate different functions and keep the results organized.

# %% editable=true slideshow={"slide_type": ""}
import numpy as np
import astropy.units as u
from munch import Munch

import tellurion as tell
# Orekit-backed functions are available on tell.*

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Next we define a specific epoch time `newyear` and a sequence of time steps every 5 minutes for one hour `prop5m1h`

# %% editable=true slideshow={"slide_type": ""}
newyear = tell.abstime('2025-01-01T00:00:00')
prop5m1h = np.linspace(5.0*u.minute, 60.0*u.minute, 12) # Step every 5 minutes for an hour

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## State and elements
#

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Create the tree structure `demoa` to hold values; `demoa.init` will hold initial values.

# %% editable=true slideshow={"slide_type": ""}
demoa = Munch()
demoa.init = Munch()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Create a position, velocity and time, representing the __[orbital state vector](https://en.wikipedia.org/wiki/Orbital_state_vectors)__

# %% editable=true slideshow={"slide_type": ""}
demoa.init.pv = [5740.13268349, 3314.06715, 0., -2.75082684, 4.76457184, 5.50165367]
demoa.init.pvt = tell.pvtcart(demoa.init.pv, newyear)
demoa.init.pvt

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Transform this into a Kepler element set and display it as a dictionary of dimensioned quantities

# %% editable=true slideshow={"slide_type": ""}
demoa.init.kep = demoa.init.pvt.kepler()  # Convert PVT to Kepler elements
demoa.init.kep.elements.to_dict()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# It is often more useful to know the altitudes of perigee and apogee than just the semimajor axis and inclination; this is accomplished with the function `allplane()`, which includes other useful information such as the mean motion

# %% editable=true slideshow={"slide_type": ""}
tell.allplane(demoa.init.kep.elements.to_dict())

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## Propagation

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### Analytical Kepler propagator

# %% [markdown] editable=true slideshow={"slide_type": ""}
# In this section, we propagate the previous state `demoa.init.pvt` using an analytical two-body propagator and store the computation in `demoa.propa`. The ephemeris table includes the initial state (this is the default but is shown explicitly here) and shows both the time from epoch and the elapsed time from the previous time step.

# %% editable=true slideshow={"slide_type": ""}
demoa.propa = Munch()
demoa.propa.gen = tell.prepare(demoa.init.pvt, 1*u.day, forceenv=tell.kepleranalytic())
demoa.propa.ephem = tell.propagate(demoa.propa.gen, prop5m1h, include_init=True)
demoa.propa.ephem

# %% [markdown] editable=true slideshow={"slide_type": ""}
# To get the full vectors, get the `'position'` or `'velocity'` columns explicitly (the suppression of the middle component of the vectors is a limitation of AstroPy that will [fixed soon](https://github.com/astropy/astropy/pull/19123)).

# %% editable=true slideshow={"slide_type": ""}
demoa.propa.ephem[3]['position']

# %% [markdown]
# Find the altitudes of perigee and apogee; because a two-body propagator is used, they only change due to the slow rotation of the coordinate system.

# %%
demoa.propa.altperapo = tell.tselements(demoa.propa.ephem, ["altper","altapo"])
demoa.propa.altperapo

# %% [markdown]
# ### Discrete states and events

# %% [markdown]
# *Discrete states* have a discrete number of values (usually two or three), usually non-numeric, that are functions of the state. For example, whether a spacecraft is eclipsed or sunlit. *Events* are transitions between states. In the example of eclipsing, there are three states, designated `u` for umbra (completely eclipsed), `p` for penumbra (partially eclipsed), and `s` for sunlit. Transitions are designated with two letters, for example, `pu` means going from penumbra into umbra. 

# %%
demoa.propa.eclipse = Munch()
demoa.propa.eclipse.events = {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []}
demoa.propa.eclipse.genevpvt = tell.prepare(demoa.init.pvt, 8*u.hour, demoa.propa.eclipse.events, output='pvt')
demoa.propa.eclipse.suntrans = demoa.propa.eclipse.genevpvt['sun transition'].ephemeris()


# %%
demoa.propa.eclipse = Munch()
demoa.propa.eclipse.events = {'altitude': 125.0*u.km, 'eclipse': [True, True], 'visibility': []}
demoa.propa.eclipse.genevpvt = tell.prepare(demoa.init.pvt, 8*u.hour, demoa.propa.eclipse.events, forceenv=tell.kepleranalytic(), output='pvt')
demoa.propa.eclipse.suntrans = demoa.propa.eclipse.genevpvt['sun transition'].ephemeris()
demoa.propa.eclipse.suntrans

# %% [markdown]
# The `[True, True]` argument in the `events` dictionary for `'eclipse'` indicates we would like to see both umbra and penumbra events. Notice the time in umbra is around 8 or 9 seconds, as this is a low-earth orbit. 

# %% [markdown]
# ### Numerical propagator

# %% [markdown] editable=true slideshow={"slide_type": ""}
# The previous examples, in `demoa.propa`, used an analytic two-body (Kepler) propagator. To do propagation with perturbations, we usually need a numerical propagator. First, we start with with just the two-body force; propagation results will be saved in `demoa.propn`.

# %% editable=true slideshow={"slide_type": ""}
demoa.propn = Munch()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# Make the propagation generator, which is necessary to create an ephemeris. A maximum time must be specified; we will make it 1 day.

# %%
demoa.propn.gen = tell.prepare(demoa.init.pvt, 1*u.day)

# %% [markdown]
# Propagate for 1 hour in 5 minute steps, and include the initial state in the ephemeris table; note the middle (`Y`) value for the position and velocity vectors has been elided, but it is present 

# %%
demoa.propn.ephem = tell.propagate(demoa.propn.gen, prop5m1h, True)
demoa.propn.ephem

# %% [markdown]
# Let's look at the difference between the and numerical results using the function `magdiff()` that computes the magnitude of the difference between two vectors.

# %%
demoa.propa.andiff = tell.magdiff(demoa.propa.ephem['position'], demoa.propn.ephem['position'])
demoa.propa.andiff


# %% [markdown]
# At each of the 13 time steps, the position difference is a maximum of about 5 microns, which is very good agreement.

# %% [markdown]
# ## Data at single points

# %% [markdown]
# This section shows how to get single-point propagation data. In the first computation, propagate to 12 hours with the previously-defined generator, which is within the maximum allowable time of 1 day for this generator.

# %% [markdown]
# Find the state at step 5. The function `.pvt()` converts the ephemeris table into a `PositionVelocityT` object, and `.to_dict()` in the last step is to make the elements result easier to read.

# %%
pvt5 = demoa.propa.ephem.pvt()[5] # demoa.propa.ephem[5].pvt() would work equally well
pvt5.cartesian.to_dict()

# %% [markdown]
# Find the state at 45 minutes using the `loc` attribute of tables

# %%
demoa.propa.ephem.loc[newyear+45*u.min]

# %%
demoa.propa.at45min = demoa.propa.ephem.loc[newyear+45*u.min]
(demoa.propa.at45min['time'], demoa.propa.at45min['position'], demoa.propa.at45min['velocity'])


# %% [markdown]
# Or, you could make a `PositionVelocityT` at that time

# %%
demoa.propa.at45min.pvt().cartesian

# %% [markdown]
# ### Time series selection at a single time point

# %% [markdown]
# This section shows different ways of getting single points from time series

# %%
demoa.tssel = Munch()
demoa.tssel.pvt15m = demoa.propn.ephem[3].pvt() # PVT for 15min by index
demoa.tssel.pvt35m = demoa.propn.ephem.loc['2025-01-01 00:35:00']
demoa.tssel.pvt45m = demoa.propn.ephem.loc[tell.abstime(45*u.min, newyear)].pvt() # PVT for 45min by relative time
demoa.tssel.pvt1h = demoa.propn.ephem[-1].pvt() # PVT at the end of the ephemeris
demoa.tssel.kep1h = demoa.tssel.pvt1h.kepler()  # Convert PVT to Kepler elements


# %%
demoa.tssel.pvt15m.cartesian

# %%
demoa.tssel.pvt35m

# %%
demoa.tssel.pvt45m.cartesian

# %%
demoa.tssel.pvt1h.time

# %%
demoa.tssel.kep1h.elements.to_dict()
