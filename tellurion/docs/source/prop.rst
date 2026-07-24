.. _propagation:

*************************************************************
Orbital Propagation (`tellurion.ork.prop`)
*************************************************************

Introduction
============

The :mod:`tellurion.ork.prop` module provides satellite orbital propagation using
multiple analytic and numerical methods. It supports:

* **Analytic propagators**: Keplerian (two-body), Brouwer-Lyddane (J2-J5 perturbations)
* **Numerical integration**: Full gravity model with optional atmospheric drag
* **Event detection**: Eclipse, altitude termination, ground station visibility
* **Flexible input**: Keplerian, circular, equinoctial, or Cartesian coordinates
* **Two-step workflow**: Prepare propagation, then propagate to arbitrary times

Getting Started
===============

Basic Usage
-----------

Propagate a satellite orbit using the two-step process::

    import astropy.units as u
    import tellurion as tell

    # Initial state as position-velocity-time
    initstate = tell.pvtcart([5740.1326835, 3314.06715, 0.,
                              -2.7508268, 4.7645718, 5.5016537],
                             tell.abstime('2025-01-01T00:00:00'))

    # Step 1: Prepare propagation
    proptime = 1.0 * u.day
    gen = tell.prepare(initstate, proptime, propagator='keplerian')

    # Step 2: Propagate to specific times
    times = [0.5*u.day, 0.75*u.day, 1.0*u.day]
    ephemeris = tell.propagate(gen, times, output='et')

    print(ephemeris)

Choosing a Propagator
---------------------

The ``prepare()`` function automatically selects a propagator based on the
initial state type and force environment::

    # Auto-selection (recommended)
    gen = tell.prepare(initstate, proptime)

    # Explicit selection
    gen = tell.prepare(initstate, proptime, propagator='keplerian')
    gen = tell.prepare(initstate, proptime, propagator='brouwer-lyddane',
                      forceenv=tell.setgravity(5, 0))
    gen = tell.prepare(initstate, proptime, propagator='numerical',
                      forceenv=tell.setgravity(20, 20))

Propagators
===========

Keplerian (Two-Body)
--------------------

The simplest analytic propagator using two-body dynamics (no perturbations).

**When to use:**

* Quick estimates
* No perturbations needed
* All orbit types (LEO, GEO, highly elliptical)
* Circular and equatorial orbits

**Characteristics:**

* Fastest computation
* No force environment required
* Osculating elements match the orbit at initial epoch only

**Example:**

.. code-block:: python

    gen = tell.prepare(initstate, proptime, propagator='keplerian')
    result = tell.propagate(gen, [0.5*u.day])

Brouwer-Lyddane (J2-J5 Perturbations)
-------------------------------------

Analytic propagator with mean element theory accounting for second through
fifth-order zonal harmonics (J2-J5). Specifically designed to handle circular
and equatorial orbits.

**When to use:**

* Medium-accuracy propagation needed
* Perturbation effects important (J2 precession, etc.)
* Fast computation essential
* Circular or equatorial orbits (key advantage over Keplerian)

**Characteristics:**

* Faster than numerical integration
* Accounts for J2-J5 perturbations
* Uses osculating-to-mean element conversion
* Much more accurate than Keplerian for perturbed orbits
* Robust handling of singular coordinate cases

**Advantages over Keplerian:**

* Includes J2-J5 perturbations (precession, etc.)
* Handles circular orbits (e=0) naturally
* Handles equatorial orbits (i=0) naturally

**Example with Elliptical Orbit:**

.. code-block:: python

    forceenv = tell.setgravity(5, 0)  # J2-J5 gravity model
    gen = tell.prepare(initstate, proptime, 
                      forceenv=forceenv, 
                      propagator='brouwer-lyddane')
    result = tell.propagate(gen, [0.5*u.day])

**Handling Circular Orbits:**

For circular orbits, pass initial state in circular element form::

    # LEO circular orbit - convert to circular elements first
    initstate_circ = leo_pvt.circular()
    gen = tell.prepare(initstate_circ, proptime, 
                      forceenv=forceenv,
                      propagator='brouwer-lyddane')

**Handling Equatorial Orbits:**

For equatorial orbits, pass initial state in equinoctial element form::

    # GEO equatorial orbit - convert to equinoctial elements
    initstate_eq = geo_pvt.equinoctial()
    gen = tell.prepare(initstate_eq, proptime,
                      forceenv=forceenv,
                      propagator='brouwer-lyddane')

Numerical Integration
---------------------

Full numerical integration with configurable gravity model and optional
atmospheric drag.

**When to use:**

* High accuracy required
* Atmospheric drag significant (low altitude orbits)
* Custom force models
* All orbit types including circular and equatorial

**Characteristics:**

* Highest accuracy
* Slowest computation
* Configurable gravity degree/order
* Optional drag force model
* Handles all singular coordinate cases

**Example:**

.. code-block:: python

    forceenv = tell.setgravity(20, 20)  # 20x20 gravity model
    gen = tell.prepare(initstate, proptime,
                      forceenv=forceenv,
                      propagator='numerical')
    result = tell.propagate(gen, [0.5*u.day])

**With Atmospheric Drag:**

.. code-block:: python

    forceenv = tell.setgravity(20, 20)
    forceenv = tell.dragforce(forceenv, cd=2.2, area=1.0, mass=1000.0)
    gen = tell.prepare(initstate, proptime,
                      forceenv=forceenv,
                      propagator='numerical')
    result = tell.propagate(gen, [0.5*u.day])

Input Coordinate Systems
========================

All propagators accept initial states in multiple coordinate representations:

Cartesian (Position-Velocity)
------------------------------

Traditional Cartesian coordinates (x, y, z, vx, vy, vz)::

    pvt = tell.pvtcart([5740.1, 3314.1, 0., -2.75, 4.76, 5.50],
                       tell.abstime('2025-01-01T00:00:00'))
    gen = tell.prepare(pvt, 1.0*u.day)

Keplerian Elements
------------------

Six classical orbital elements (sma, ecc, inc, argper, raan, ma/ta)::

    kep = tell.kepler({'sma': 6720*u.km, 'ecc': 0.001,
                       'inc': 51.6*u.deg, 'argper': 0*u.deg,
                       'raan': 0*u.deg, 'ma': 0*u.deg},
                      tell.abstime('2025-01-01T00:00:00'))
    gen = tell.prepare(kep, 1.0*u.day)

Circular Elements
-----------------

For non-equatorial circular orbits (sma, cex, cey, inc, raan, mla/tla)::

    circ = leo_pvt.circular()
    gen = tell.prepare(circ, 1.0*u.day)

Equinoctial Elements
--------------------

For all orbits including equatorial (sma, ex, ey, hx, hy, ml/tl)::

    eq = geo_pvt.equinoctial()
    gen = tell.prepare(eq, 1.0*u.day)

Event Detection
===============

Detect and record events during propagation using the ``events`` dictionary.

Eclipse Detection
-----------------

Detect when the satellite enters Earth's shadow::

    events = {
        'altitude': 125.0*u.km,
        'eclipse': [True, True],  # [detect umbra, detect penumbra]
        'visibility': []
    }
    gen = tell.prepare(initstate, proptime, events=events)

The ephemeris will include an ``'eclipse'`` column with status at each time:

* ``'s'`` — full sunlight
* ``'p'`` — penumbra (partial eclipse)
* ``'u'`` — umbra (total eclipse)
* Space characters separate time steps

**Example:**

    result = tell.propagate(gen, times, output='pvt')
    print(result.aux['eclipse'])  # 's p p u u p p s s'

The generator's ``'sun transition'`` entry contains the times of eclipse
transitions::

    sun_trans = gen['sun transition']
    print(sun_trans.ephemeris())

Altitude Termination
--------------------

Stop propagation if altitude drops below a threshold::

    events = {'altitude': 125.0*u.km, 'eclipse': [], 'visibility': []}
    gen = tell.prepare(initstate, proptime, events=events)

Propagation will terminate when the satellite reaches the specified altitude.

Ground Station Visibility
--------------------------

Detect visibility from one or more Earth locations::

    from astropy.coordinates import EarthLocation
    import astropy.units as u

    # Define ground stations
    stations = [
        {'name': 'Station 1', 'lon': 0.0*u.deg, 'lat': 45.0*u.deg, 'elevation': 0*u.m},
        {'name': 'Station 2', 'lon': 180.0*u.deg, 'lat': -45.0*u.deg, 'elevation': 0*u.m}
    ]

    events = {
        'altitude': 125.0*u.km,
        'eclipse': [],
        'visibility': stations
    }
    gen = tell.prepare(initstate, proptime, events=events)

Propagation Output
==================

The ``prepare()`` function returns a generator dictionary with propagation
setup. The ``propagate()`` function then computes results.

Output Formats
--------------

Control output format with the ``output`` parameter::

    # Ephemeris table (default)
    ephem = tell.propagate(gen, times, output='et')

    # Position-Velocity-Time object
    pvt = tell.propagate(gen, times, output='pvt')

    # Orekit SpacecraftState objects
    states = tell.propagate(gen, times, output='ss')

Ephemeris Table (ET)
~~~~~~~~~~~~~~~~~~~~

AstroPy TimeSeries with columns for position, velocity, time, and any
requested auxiliary data::

    ephem = tell.propagate(gen, times, output='et')
    print(ephem)
    print(ephem['position'])  # Cartesian position
    print(ephem['velocity'])  # Cartesian velocity
    print(ephem['sunlight']) # Eclipse status (if detected)

Position-Velocity (PVT)
~~~~~~~~~~~~~~~~~~~~~~~

Native ``PositionVelocityT`` object with methods for coordinate conversion::

    pvt = tell.propagate(gen, times, output='pvt')
    print(pvt.cartesian)
    print(pvt.spherical)
    kep = pvt.kepler()

SpacecraftState (SS)
~~~~~~~~~~~~~~~~~~~~

Orekit ``SpacecraftState`` objects (advanced use)::

    states = tell.propagate(gen, times, output='ss')

Propagation Example: Complete Workflow
=======================================

Here's a complete example comparing propagators on a LEO orbit with
eclipse detection::

    import numpy as np
    import astropy.units as u
    import tellurion as tell

    # Define initial state
    initstate = tell.pvtcart([5740.1326835, 3314.06715, 0.,
                              -2.7508268, 4.7645718, 5.5016537],
                             tell.abstime('2025-01-01T00:00:00'))

    proptime = 8.0 * u.hour
    times = np.linspace(1.0*u.hour, 6.0*u.hour, 32)

    # Propagation events
    events = {
        'altitude': 125.0*u.km,
        'eclipse': [True, True],
        'visibility': []
    }

    # Method 1: Keplerian (fast, no perturbations)
    gen_kep = tell.prepare(initstate, proptime, events=events,
                          propagator='keplerian')
    result_kep = tell.propagate(gen_kep, times, output='et')

    # Method 2: Brouwer-Lyddane (medium speed, J2-J5)
    forceenv_bl = tell.setgravity(5, 0)
    gen_bl = tell.prepare(initstate, proptime, events=events,
                         forceenv=forceenv_bl,
                         propagator='brouwer-lyddane')
    result_bl = tell.propagate(gen_bl, times, output='et')

    # Method 3: Numerical (slow, high accuracy)
    forceenv_num = tell.setgravity(20, 20)
    gen_num = tell.prepare(initstate, proptime, events=events,
                          forceenv=forceenv_num,
                          propagator='numerical')
    result_num = tell.propagate(gen_num, times, output='et')

    # Compare results
    print("Keplerian position at t=4 hours:")
    print(result_kep[16]['position'])
    print("\nBrouwer-Lyddane position at t=4 hours:")
    print(result_bl[16]['position'])
    print("\nNumerical position at t=4 hours:")
    print(result_num[16]['position'])

    # Check eclipse status
    print("\nEclipse transitions (Keplerian):")
    print(gen_kep['sun transition'].ephemeris())

API Reference
=============

.. autoapi:: tellurion.core.prop
   :noindex:

.. seealso::

   | :ref:`element` — orbital element representations
   | :ref:`posvel` — position-velocity-time objects
   | :ref:`force` — force environment configuration
