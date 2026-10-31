.. _propagation:

*************************************************************
Propagation
*************************************************************

Introduction
============

The :mod:`tellurion` module provides satellite orbital propagation using
multiple analytic and numerical methods. It supports:

* **Analytic propagators**: Keplerian (two-body), Brouwer-Lyddane (J2-J5 perturbations)
* **Semi-analytical propagator**: DSST (Draper Semi-analytical Satellite Theory) with full gravity model
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

    ephemeris.pprint_all(max_width=-1)

Choosing a Propagator
---------------------

The :func:`tellurion.prepare` function automatically selects a propagator based on the
initial state type and force environment::

    # Auto-selection (recommended)
    gen = tell.prepare(initstate, proptime)

    # Explicit selection
    gen = tell.prepare(initstate, proptime, propagator='keplerian')
    gen = tell.prepare(initstate, proptime, propagator='brouwer-lyddane',
                      forceenv=tell.setgravity(5, 0))
    gen = tell.prepare(initstate, proptime, propagator='dsst',
                      forceenv=tell.setgravity(20, 20))
    gen = tell.prepare(initstate, proptime, propagator='numerical',
                      forceenv=tell.setgravity(20, 20))

.. _propagators:

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
* The default force environment is sufficient, so ``forceenv`` need not be explicitly specified
* Osculating elements match the orbit at initial epoch only

**Example:**

.. code-block:: python

    gen = tell.prepare(initstate, proptime, propagator='keplerian')
    result = tell.propagate(gen, [0.5*u.day])

Brouwer-Lyddane (zonal perturbations)
-------------------------------------

Analytic propagator with mean element theory accounting for second through
fifth-order zonal harmonics (:math:`J_2` through :math:`J_5`). Specifically designed to handle circular
and equatorial orbits.

**When to use:**

* Medium-accuracy propagation needed
* Perturbation effects important (:math:`J_2` precession of the node, etc.)
* Fast computation essential
* Circular or equatorial orbits (key advantage over Keplerian)

**Characteristics:**

* Faster than numerical integration
* Accounts for :math:`J_2` through :math:`J_5` zonal harmonic perturbations
* Uses osculating-to-mean element conversion
* Much more accurate than Keplerian for perturbed orbits
* Robust handling of singular coordinate cases

**Advantages over Keplerian:**

* Includes :math:`J_2` through :math:`J_5` geopotential perturbations
* Handles circular orbits (:math:`e=0`) naturally
* Handles equatorial orbits (:math:`i=0`) naturally

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

DSST (Semi-Analytical)
----------------------

Semi-analytical propagator using Draper Semi-analytical Satellite
Theory.  Combines numerical integration with perturbation averaging
for efficient long-term propagation with full gravity model. Note: the
DSST model will sometimes crash; the location of this error has not
been determined.

**When to use:**

* Long-duration propagation (days to weeks)
* Full gravity model needed (not just :math:`J_2`--:math:`J_5`)
* Better speed/accuracy tradeoff than pure numerical integration
* All orbit types including circular and equatorial

**Characteristics:**

* Faster than numerical integration (2-10x speedup)
* Full spherical harmonic gravity model support
* Accounts for zonal and tesseral harmonics
* Semi-analytical (not pure analytic)
* Handles all singular coordinate cases
* Excellent for long-term propagation

**Advantages over Numerical Integration:**

* Significantly faster (2-10x depending on gravity degree/order)
* Can handle extended propagation periods
* Still uses full gravity model
* Smaller numerical errors due to averaging over orbital period

**Advantages over Brouwer-Lyddane:**

* Supports full tesseral harmonics
* More accurate for long-term propagation
* Better for higher-degree gravity models

**Example:**

.. code-block:: python

    forceenv = tell.setgravity(20, 20)  # Full 20x20 gravity model
    gen = tell.prepare(initstate, proptime,
                      forceenv=forceenv,
                      propagator='dsst')
    result = tell.propagate(gen, [2.0*u.day])

**Long-Term Propagation Example:**

.. code-block:: python

    # Propagate for one month with DSST
    proptime = 30.0 * u.day
    forceenv = tell.setgravity(20, 20)
    
    gen = tell.prepare(initstate, proptime,
                      forceenv=forceenv,
                      propagator='dsst')
    
    # Request positions at multiple times
    times = np.linspace(1.0*u.day, 30.0*u.day, 30)
    result = tell.propagate(gen, times, output='et')

**Comparison: DSST vs Numerical for LEO (50 day propagation):**

.. code-block:: python

    import time
    
    proptime = 50.0 * u.day
    forceenv = tell.setgravity(20, 20)
    times = np.linspace(1.0*u.day, 50.0*u.day, 50)
    
    # Numerical propagation
    t0 = time.time()
    gen_num = tell.prepare(initstate, proptime, forceenv=forceenv,
                          propagator='numerical')
    result_num = tell.propagate(gen_num, times, output='et')
    t_num = time.time() - t0
    
    # DSST propagation
    t0 = time.time()
    gen_dsst = tell.prepare(initstate, proptime, forceenv=forceenv,
                           propagator='dsst')
    result_dsst = tell.propagate(gen_dsst, times, output='et')
    t_dsst = time.time() - t0
    
    print(f"Numerical: {t_num:.2f}s")
    print(f"DSST: {t_dsst:.2f}s")
    print(f"Speedup: {t_num/t_dsst:.1f}x")

**Handling Circular Orbits:**

For circular orbits, convert to circular elements first::

    initstate_circ = leo_pvt.circular()
    gen = tell.prepare(initstate_circ, proptime,
                      forceenv=forceenv,
                      propagator='dsst')

**Handling Equatorial Orbits:**

For equatorial orbits, convert to equinoctial elements::

    initstate_eq = geo_pvt.equinoctial()
    gen = tell.prepare(initstate_eq, proptime,
                      forceenv=forceenv,
                      propagator='dsst')

Numerical Integration
---------------------

Full numerical integration with configurable gravity model and optional
atmospheric drag.

**When to use:**

* High accuracy required (short propagation periods)
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

    forceenv = tell.setgravity(20, 20, mass=1000.0)
    forceenv = tell.dragforce(forceenv, dragcoef=2.2, dragarea=1.0)
    gen = tell.prepare(initstate, proptime,
                      forceenv=forceenv,
                      propagator='numerical')
    result = tell.propagate(gen, [0.5*u.day])

Force Model Recipes
===================

Use these practical force-environment recipes as starting points.

Two-Body Baseline (fastest)
---------------------------

Use analytic Keplerian propagation with no perturbation model::

    gen = tell.prepare(initstate, 1.0*u.day, propagator='keplerian')

J2-J5 Analytic Screening
------------------------

Capture dominant zonal perturbations with Brouwer-Lyddane::

    forceenv = tell.setgravity(5, 0)  # J2-J5 zonal terms
    gen = tell.prepare(initstate, 2.0*u.day,
                      forceenv=forceenv,
                      propagator='brouwer-lyddane')

DSST Long-Arc Tradeoff
----------------------

Use DSST for longer arcs with full gravity and better speed than full
numerical integration::

    forceenv = tell.setgravity(20, 20)
    gen = tell.prepare(initstate, 15.0*u.day,
                      forceenv=forceenv,
                      propagator='dsst')

Numerical + Drag (LEO)
----------------------

Use numerical propagation when atmospheric drag is mission-significant::

    forceenv = tell.setgravity(20, 20, mass=1000.0)
    forceenv = tell.dragforce(forceenv, atmdensname='hp',
                             dragcoef=2.2, dragarea=1.0)
    gen = tell.prepare(initstate, 1.0*u.day,
                      forceenv=forceenv,
                      propagator='numerical')

``atmdensname`` can be ``'hp'`` (Harris-Priester), ``'dtm'`` (DTM2000), or
``'msis'`` (NRLMSISE00).

.. _propcomp:

Propagator Comparison
=====================

Quick reference for choosing a propagator:

.. list-table::
   :header-rows: 1
   :widths: 20, 15, 15, 15, 15, 20

   * - Propagator
     - Speed
     - Accuracy
     - Max Gravity
     - Best For
     - Duration

   * - Keplerian
     - ★★★★★
     - ★★☆☆☆
     - None
     - Quick estimates
     - Any

   * - Brouwer-Lyddane
     - ★★★★☆
     - ★★★☆☆
     - J2-J5 only
     - Medium accuracy, mean elements
     - Days

   * - DSST
     - ★★★☆☆
     - ★★★★☆
     - Full (20x20+)
     - Long-term, full gravity
     - Weeks-months

   * - Numerical
     - ★★☆☆☆
     - ★★★★★
     - Full (20x20+)
     - High accuracy, drag
     - Hours-days

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

.. _event:

Discrete state and event
========================

An *event* is a transition of some *discrete state* from one value to another. For example, a transition from full sunlight to penumbra is an eclipse event. Tellurion can record several discrete states for any ephemeris, and create an ephemeris table of the transition events for that class of state.

To detect and record states and events during propagation, using the ``events`` dictionary when calling
:func:`~tellurion.prepare`.

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
    sun_trans.ephemeris().pprint_all(max_width=-1)

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

State Transition Matrix (STM)
-----------------------------

Request the final state transition matrix by setting ``'stm': True`` in the
``events`` dictionary::

    events = {
        'altitude': 125.0*u.km,
        'eclipse': [],
        'visibility': [],
        'stm': True,
    }
    gen = tell.prepare(initstate, proptime, events=events,
                       forceenv=tell.setgravity(20, 20),
                       propagator='numerical')
    _ = tell.propagate(gen, [proptime], output='pvt')

After propagation, the final STM is available in ``gen['final']``:

.. code-block:: python

    phi = gen['final']['stm']
    print(phi.shape)  # Typically (6, 6)

If drag-parameter partial derivatives are enabled in the force model, the
final parameter Jacobian is also stored in:

.. code-block:: python

    pjac = gen['final']['parameters jacobian']

Use the STM to estimate sensitivity of final Cartesian state to small changes
in the initial Cartesian state.

See :doc:`tutorials/stm1` for a step-by-step STM workflow.

Propagation Output
==================

The :func:`~tellurion.prepare` returns a generator dictionary with propagation
setup. The :func:`~tellurion.propagate` function then computes the position and velocity over the requested time steps as a :class:`~tellurion.PositionVelocityT` or an :ref:`ephemeris`.

Output Formats
--------------

Control output format with the ``output`` parameter::

    # Ephemeris table (default)
    ephem = tell.propagate(gen, times, output='et')

    # Position-Velocity-Time object
    pvt = tell.propagate(gen, times, output='pvt')

Ephemeris Table (ET)
~~~~~~~~~~~~~~~~~~~~

A :class:`~astropy.timeseries.TimeSeries` with columns for position, velocity, time, and any
requested auxiliary data::

    ephem = tell.propagate(gen, times, output='et')
    print(ephem)
    print(ephem['position'])  # Cartesian position
    print(ephem['velocity'])  # Cartesian velocity
    print(ephem['sunlight']) # Eclipse status (if detected)

Position-Velocity (PVT)
~~~~~~~~~~~~~~~~~~~~~~~

Native :class:`~tellurion.PositionVelocityT` object with methods for coordinate conversion::

    pvt = tell.propagate(gen, times, output='pvt')
    print(pvt.cartesian)
    print(pvt.spherical)
    kep = pvt.kepler()

Relative Frames
---------------

After generating states with :func:`tellurion.propagate`, you can express one
spacecraft relative to another in local orbital frames:

* :func:`tellurion.lvlh` (LVLH/RSW)
* :func:`tellurion.ntw` (NTW)

and map relative states back to inertial with:

* :func:`tellurion.inv_lvlh`
* :func:`tellurion.inv_ntw`

For frame definitions, conventions, inverse transforms, and complete examples,
see :ref:`posvel-relative-frames`.


Propagation examples
=======================================

Inertial motion of a single satellite
-------------------------------------

Here's a complete example comparing propagators on a LEO orbit with
eclipse detection:

.. jupyter-execute::

    import numpy as np
    import astropy.units as u
    import tellurion as tell

    # Define initial state
    initstate = tell.pvtcart([5740.1326835, 3314.06715, 0.,
                              -2.7508268, 4.7645718, 5.5016537],
                             tell.abstime('2025-01-01T00:00:00'))

    proptime = 8.0 * u.hour
    # Ephemeris steps every half hour for 6 hours
    times = np.linspace(0.5*u.hour, 6.0*u.hour, 12)

    # Propagation events
    events = {
        'altitude': 125.0*u.km,
        'eclipse': [True, True],  # [detect umbra, detect penumbra]
        'visibility': []
    }


Method 1: Keplerian (fast, no perturbations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    gen_kep = tell.prepare(initstate, proptime, events=events,
                          propagator='keplerian')
    result_kep = tell.propagate(gen_kep, times, output='et')
    result_kep.pprint_all(max_width=-1)


Method 2: Brouwer-Lyddane (medium speed, J2-J5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    forceenv_bl = tell.setgravity(5, 0)
    gen_bl = tell.prepare(initstate, proptime, events=events,
                          forceenv=forceenv_bl,
                          propagator='brouwer-lyddane')
    result_bl = tell.propagate(gen_bl, times, output='et')


Method 3: DSST (medium speed, full gravity, long-term)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    forceenv_dsst = tell.setgravity(20, 20)
    gen_dsst = tell.prepare(initstate, proptime, events=events,
                            forceenv=forceenv_dsst,
                            propagator='dsst')
    result_dsst = tell.propagate(gen_dsst, times, output='et')


Method 4: Numerical (slow, high accuracy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    forceenv_num = tell.setgravity(20, 20)
    gen_num = tell.prepare(initstate, proptime, events=events,
                           forceenv=forceenv_num,
                           propagator='numerical')
    result_num = tell.propagate(gen_num, times, output='et')


Compare results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    print("Keplerian position at t=4 hours:")
    print(result_kep[8]['position'])
    print("\nBrouwer-Lyddane position at t=4 hours:")
    print(result_bl[8]['position'])
    print("\nDSST position at t=4 hours:")
    print(result_dsst[8]['position'])
    print("\nNumerical position at t=4 hours:")
    print(result_num[8]['position'])

    # Check eclipse status
    print("\nEclipse transitions (Keplerian):")
    gen_kep['sun transition'].pprint_exclude_names =['position', 'velocity']
    gen_kep['sun transition'].pprint_all(max_width=-1)

Motion of a secondary satellite relative to a primary
-----------------------------------------------------

First, initialize the primary and secondary states and propagate them.

.. jupyter-execute::

    import numpy as np
    import astropy.units as u
    import tellurion as tell

    icref = np.asarray([5740.1326835, 3314.06715, 0.0, -2.7508268, 4.7645718, 5.5016537])
    icsec = icref + np.asarray([1.0, 0.0, 0.0, 0.0, -0.001, 0.0])

    init_p = tell.pvtcart(icref, tell.abstime('2025-01-01T00:00:00'))
    init_s = tell.pvtcart(icsec, tell.abstime('2025-01-01T00:00:00'))

    steps = np.linspace(5.0 * u.minute, 60.0 * u.minute, 12)
    primary_series = tell.propagate(tell.prepare(init_p, 1.0 * u.hour), steps, include_init=True, output='pvt')
    secondary_series = tell.propagate(tell.prepare(init_s, 1.0 * u.hour), steps, include_init=True, output='pvt')

Secondary ephemeris
~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

   secondary_series.ephemeris().pprint_all(max_width=-1)

Relative frame ephemeris
~~~~~~~~~~~~~~~~~~~~~~~~

.. jupyter-execute::

    rel_rsw = tell.lvlh(secondary_series, primary_series, tell.siunits)
    rel_rsw.ephemeris().pprint_all(max_width=-1)

    rel_ntw = tell.ntw(secondary_series, primary_series, tell.siunits)
    rel_ntw.ephemeris().pprint_all(max_width=-1)

.. _multibody:

Multibody (cislunar and interplanetary) propagation
===================================================

:func:`tellurion.mbprepare` and :func:`tellurion.mbpropagate` propagate a
spacecraft under the point-mass gravity of several named bodies at once. They
are separate from :func:`~tellurion.prepare` and :func:`~tellurion.propagate`
(which remain the Earth-orbit tools), and use the same two-step workflow.

How it works
------------

* The state is integrated in the inertial, solar-system-barycentric ICRF frame.
  This single frame is valid near Earth, near the Moon and in deep space, so
  the propagation never switches origin.
* Each body in the ``bodies`` list contributes one Orekit
  ``SingleBodyAbsoluteAttraction``, with positions from the Orekit ephemerides
  (``CelestialBodyFactory`` names such as ``"Sun"``, ``"Earth"``, ``"Moon"``,
  ``"Jupiter"``).  The default is Sun, Earth and Moon.
* The initial state must carry a ``frame`` descriptor
  (see :ref:`posvel-frames`). It is converted to barycentric ICRF internally,
  and results are tagged with their frame and may be requested in any
  supported frame with ``frame=``.  Supported frames are Earth-centered
  ``FrameDescriptor("Earth", "EME2000")`` and barycentric
  ``FrameDescriptor("Solar System Barycenter", "ICRF")``.
* Only point-mass gravity from the listed bodies is modeled: no Earth
  gravity field, drag, radiation pressure, or events.

Cislunar example
----------------

A highly elliptical orbit reaching lunar distance, started in Earth-centered
EME2000 and propagated for six days under Sun, Earth and Moon gravity::

    import numpy as np
    import astropy.units as u
    import tellurion as tell

    eci = tell.FrameDescriptor("Earth", "EME2000")
    t0 = tell.abstime("2026-04-02T02:00:00")

    kep = tell.kepler(
        tell.allplane({"altper": 200*u.km, "altapo": 400000*u.km,
                       "inc": 28.5*u.deg, "argper": 0*u.deg,
                       "raan": 0*u.deg, "ma": 0*u.deg}),
        t0)
    init = tell.pvt(kep)
    init.frame = eci                      # tag the initial state

    gen = tell.mbprepare(init, 6*u.day, ["Sun", "Earth", "Moon"])
    out = tell.mbpropagate(gen, np.linspace(1*u.day, 6*u.day, 6),
                           include_init=True, frame=eci)
    out.ephemeris().pprint_all(max_width=-1)

Omitting ``frame=`` returns the states in barycentric ICRF, which is the
natural frame for plotting transfers that leave the Earth-Moon system::

    bary = tell.mbpropagate(gen, [6*u.day])
    print(bary.frame)                     # Solar System Barycenter / ICRF

Interplanetary example
----------------------

For a departure or heliocentric cruise, list the planets that matter and use a
longer span and looser step limit.  Here a state at Earth, on a hyperbolic
escape trajectory, is propagated for 200
days under the Sun, Earth, Mars and Jupiter::

    pv_eci = np.array([7000.0, 0.0, 0.0, 0.0, 10.9, 0.0])   # km, km/s
    init = tell.pvtcart(pv_eci, tell.abstime("2026-10-01T00:00:00"),
                        frame=tell.FrameDescriptor("Earth", "EME2000"))
    gen = tell.mbprepare(init, 200*u.day,
                         ["Sun", "Earth", "Moon", "Mars", "Jupiter"],
                         maxstep=86400.0)
    cruise = tell.mbpropagate(gen, np.linspace(10*u.day, 200*u.day, 20))
    cruise.ephemeris().pprint_all(max_width=-1)

Position is in kilometres from the solar system barycenter.

Converting between frames
-------------------------

Any supported state can be moved between Earth-centered and barycentric frames
with :func:`tellurion.transform_pvt` (including the origin's velocity)::

    bary = tell.FrameDescriptor("Solar System Barycenter", "ICRF")
    earth_view = tell.transform_pvt(cruise, tell.FrameDescriptor("Earth", "EME2000"))
    back = tell.transform_pvt(earth_view, bary)

Using Horizons or OEM states
----------------------------

States from :func:`tellurion.horizons` and :func:`tellurion.oempvt` arrive
with frame descriptors, and can be used as initial states or as reference
trajectories to compare against.  A Horizons query with ``center="@399"`` and
``refplane="earth"`` is labelled ``"ICRF/Earth-equator"``, which is not yet
resolvable by Orekit; because that frame agrees with EME2000 to a few tens of
milliarcseconds, re-tag it before propagating::

    ref = tell.horizons("-1024", "@399",
                        "2026-04-02T02:00:00", "2026-04-05T02:00:00", "1h")
    init = ref[0]
    init.frame = tell.FrameDescriptor("Earth", "EME2000")

    gen = tell.mbprepare(init, 3*u.day, ["Sun", "Earth", "Moon"])
    mine = tell.mbpropagate(gen, np.linspace(1*u.hour, 3*u.day, 72))
    mine = tell.transform_pvt(mine, tell.FrameDescriptor("Earth", "EME2000"))

Comparing ``mine`` with the later Horizons or OEM rows then shows the effect of
unmodeled forces (thrust, solar radiation pressure, Earth's non-spherical
field), which is typically what limits cislunar accuracy.

Tuning and limitations
----------------------

* ``rtol``/``atol`` (default ``1e-12``/``1e-3``) control the
  Dormand-Prince 8(5,3) integrator; ``minstep`` and ``maxstep`` are in seconds.
* Frames that Orekit cannot resolve (Moon-centered axes, Horizons or OEM
  frames other than Earth/EME2000 and barycentric ICRF) raise ``ValueError``
  rather than being guessed.
* Orekit data (ephemerides, in particular DE) must be installed; see
  :doc:`installation`.

API Reference
=============

.. autofunction:: tellurion.mbprepare

.. autofunction:: tellurion.mbpropagate

.. autofunction:: tellurion.prepare

.. autofunction:: tellurion.propagate

.. autofunction:: tellurion.setgravity

.. autofunction:: tellurion.dragforce

.. seealso::

   | :ref:`element` — orbital element representations
   | :ref:`posvel` — position-velocity-time objects
   | :ref:`lambert` — Lambert transfer workflows
   | :doc:`tutorials/lambert1` — Lambert transfer tutorial notebook
   | :doc:`tutorials/stm1` — state transition matrix tutorial notebook
